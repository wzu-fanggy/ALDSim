#core/layer.py
# posA 状态定义（论文物理语义）
# 0 = 空
# 1 = Al–CH3
# 2 = –OH
# 3 = Al–O–Al（可选，后续扩展）
from __future__ import annotations
from typing import Optional, List
from .cell import C2DCell


class CLayer:
    """
    薄膜层结构（100×100 柱状结构，双向链表实现）。
    论文中每个 (i,j) 位置是一根“柱”，柱内由 C2DCell 双向链表构成。
    """

    def __init__(self, simusize: int = 100):
        self.SIMUSIZE = simusize

        # 100×100 grid，每个元素都是一根“柱”（链表）
        self.cellBottom: List[List[Optional[C2DCell]]] = [
            [None for _ in range(simusize)] for _ in range(simusize)
        ]
        self.cellTop: List[List[Optional[C2DCell]]] = [
            [None for _ in range(simusize)] for _ in range(simusize)
        ]

        self.cell_area = 1.0

        # 保存当前最薄的位置（用于优化）
        self.lowLayerNum = 0
        self.lowLayerX = 0
        self.lowLayerY = 0

    # ===============================================================
    #  初始化 Si 层（论文 4.4.1）
    # ===============================================================
    # def InitSiLayer(self):
    #     for i in range(self.SIMUSIZE):
    #         for j in range(self.SIMUSIZE):
    #             cell = C2DCell(currNum=0)
    #             cell.posA = 2  # 2 = OH
    #             cell.posB = 0
    #             cell.posC = 0
    #             self.cellBottom[i][j] = cell
    #             self.cellTop[i][j] = cell
    #
    #     self.lowLayerNum = 0
    def InitSiLayer(self):
        """
        初始化 Si(001)-OH 表面
        物理含义：
        - 第 0 层为 Si 基底表面
        - 所有表面位点被 OH 饱和
        """
        for i in range(self.SIMUSIZE):
            for j in range(self.SIMUSIZE):
                cell = C2DCell(currNum=0)

                # ---------- Si(001)-OH 初始表面 ----------
                cell.posA = 2  # 2 = OH（Si–OH）
                cell.posB = 0
                cell.posC = 0

                # 标记这是“基底层”（非常重要，后续可用）
                cell.is_substrate = True

                self.cellBottom[i][j] = cell
                self.cellTop[i][j] = cell

        self.lowLayerNum = 0

    # ===============================================================
    #  内部函数：周期性边界条件
    # ===============================================================
    def _wrap(self, i, j):
        si = self.SIMUSIZE
        if i >= si: i -= si
        if j >= si: j -= si
        if i < 0: i += si
        if j < 0: j += si
        return i, j

    # ===============================================================
    #  在(i,j)处创建新层（在链表顶部加一个 cell）
    # ===============================================================
    def CreateNewCell(self, i: int, j: int):
        i, j = self._wrap(i, j)

        top = self.cellTop[i][j]
        new_cell = C2DCell(currNum=top.currNum + 1)

        # 双向连接
        new_cell.previous = top
        top.next = new_cell

        # 更新顶部
        self.cellTop[i][j] = new_cell

        # 若增加了(i,j)处的高度，更新最薄层
        if new_cell.currNum > self.lowLayerNum:
            # 重新扫描最薄值（论文图 4-10 逻辑）
            self._update_lowest_layer()

    # ===============================================================
    #  填充某个 cell（论文 FillCell）
    # ===============================================================
    def FillCell(self, i: int, j: int, k: int, position: int, value: int, nCycle: int):
        """
        (i,j,k) 为柱(i,j)的第 k 层。
        position: 0=A, 1=B, 2=C
        value: 0/1/2（物理含义由化学反应逻辑决定）
        """

        if not (0 <= position <= 2 and 0 <= value <= 2):
            return -1

        i, j = self._wrap(i, j)

        # 若层数不足，自动创建
        while self.cellTop[i][j].currNum < k:
            self.CreateNewCell(i, j)

        # 在链表中查找第 k 层（自顶/自底选最短）
        target = self._find_cell(i, j, k)
        if target is None:
            return -1

        # 填充
        if position == 0:
            target.posA = value
        elif position == 1:
            target.posB = value
        elif position == 2:
            target.posC = value

        target.nCycle = nCycle
        return 0

    # ===============================================================
    #  内部：查找第 k 层 cell（论文中最关键的 O(h/4) 算法）
    # ===============================================================
    def _find_cell(self, i: int, j: int, k: int) -> Optional[C2DCell]:

        bottom = self.cellBottom[i][j]
        top = self.cellTop[i][j]

        if bottom is None or top is None:
            return None

        if k < bottom.currNum or k > top.currNum:
            return None

        # 从近的一端搜索
        mid = (bottom.currNum + top.currNum) // 2
        if k > mid:
            temp = top
            while temp and temp.currNum != k:
                temp = temp.previous
        else:
            temp = bottom
            while temp and temp.currNum != k:
                temp = temp.next

        return temp

    # ===============================================================
    # 删除某一层的 cell（论文 DeleteCell）
    # ===============================================================
    def DeleteCell(self, i: int, j: int, k: int):
        i, j = self._wrap(i, j)

        target = self._find_cell(i, j, k)
        if target is None:
            return False

        # 只有一个节点
        if target.previous is None and target.next is None:
            self.cellBottom[i][j] = None
            self.cellTop[i][j] = None
            return True

        # 删除顶部
        if target.next is None:
            self.cellTop[i][j] = target.previous
            target.previous.next = None
            return True

        # 删除底部
        if target.previous is None:
            self.cellBottom[i][j] = target.next
            target.next.previous = None
            return True

        # 删除中间
        target.previous.next = target.next
        target.next.previous = target.previous
        return True

    # ===============================================================
    # 删除整层（论文 DeleteLayer）
    # ===============================================================
    def DeleteLayer(self, k: int):
        for i in range(self.SIMUSIZE):
            for j in range(self.SIMUSIZE):
                self.DeleteCell(i, j, k)

        self._update_lowest_layer()

    # ===============================================================
    # 获取 cell 值（论文 GetCellValue）
    # ===============================================================
    def GetCellValue(self, i: int, j: int, k: int, index: int):
        """
        index = 0:A, 1:B, 2:C, 3:nCycle
        """
        if not (0 <= index <= 3):
            return -1

        i, j = self._wrap(i, j)

        cell = self._find_cell(i, j, k)
        if cell is None:
            return -1

        if index == 0:
            return cell.posA
        elif index == 1:
            return cell.posB
        elif index == 2:
            return cell.posC
        elif index == 3:
            return cell.nCycle

        return -1

    # ===============================================================
    # 查找最薄层（论文 FindLowLayer）
    # ===============================================================
    def FindLowLayer(self):
        return self.lowLayerNum

    # ===============================================================
    # 内部：更新最薄层位置（论文 4-10 逻辑）
    # ===============================================================
    def _update_lowest_layer(self):
        min_k = float("inf")
        min_x = 0
        min_y = 0

        for i in range(self.SIMUSIZE):
            for j in range(self.SIMUSIZE):
                top = self.cellTop[i][j]
                if top is None:
                    continue
                if top.currNum < min_k:
                    min_k = top.currNum
                    min_x = i
                    min_y = j

        self.lowLayerNum = min_k
        self.lowLayerX = min_x
        self.lowLayerY = min_y

    # ===============================================================
    # 获取 (i, j) 柱的当前高度（最高层层号）
    # ===============================================================
    def GetHeight(self, i: int, j: int) -> int:
        i, j = self._wrap(i, j)
        top = self.cellTop[i][j]
        if top is None:
            return 0
        return top.currNum

    # ===============================================================
    #  气体到达事件（对应论文到达事件）
    # ===============================================================
    def GasArrival(self, i: int, j: int, gas: str):
        """
        对应论文“到达事件”：仅在最顶层设置 posB 或 posC 标志，
        真正的化学反应由 React() 完成。
        """
        i, j = self._wrap(i, j)
        top = self.cellTop[i][j]
        if top is None:
            return

        # TMA 到达
        if gas == "TMA":
            # posB 作为“到达吸附”标志：1=吸附
            top.posB = 1

        # H2O 到达
        elif gas == "H2O":
            top.posC = 1

    # ===============================================================
    #  反应事件（对应论文反应事件）
    # ===============================================================
    def React(self, i: int, j: int, gas: str):
        """
        对应论文逻辑：
        TMA + OH -> AlCH3        (增长一层)
        H2O + AlCH3 -> OH        (不增长)
        """

        i, j = self._wrap(i, j)
        top = self.cellTop[i][j]
        if top is None:
            return

        # -----------------------------
        #  TMA 反应 (形成 Al–CH3，薄膜增长)
        # -----------------------------
        if gas == "TMA":
            # posA = 2 代表 –OH； posB = 1 代表吸附了 TMA
            if top.posA == 2 and top.posB == 1:
                # 反应：形成 AlCH3
                top.posA = 1       # 1 = AlCH3
                top.posB = 0       # 清除吸附标志

                # 增长薄膜：生成新层
                self.CreateNewCell(i, j)

        # -----------------------------
        #  H₂O 反应 (AlCH3 → OH，不增长)
        # -----------------------------
        elif gas == "H2O":
            if top.posA == 1 and top.posC == 1:
                # 反应：AlCH3 → OH
                top.posA = 2       # 2 = OH
                top.posC = 0       # 清除吸附标志

