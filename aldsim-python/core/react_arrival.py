# core/react_arrival.py
from __future__ import annotations
import math
import random
from typing import List


INFINITE_TIME = 1e30


class CReactantArrival:
    """
    前驱体到达事件时间表（对应论文 CReactantArrival）。
    每个 2D 单元 (i,j) 都有一个到达时间，用 100×100 浮点数组保存。

    论文中提到事件频繁更新，因此数组结构优于链表，访问成本低。
    """

    def __init__(self, simusize: int = 100):
        self.SIMUSIZE = simusize

        # 参数
        self.T = 0.0         # 温度 (K)
        self.P = 0.0         # 压力 (Pa)
        self.S = 0.0         # 2D 单元面积 (m^2)
        self.M = 0.0         # 分子量 (kg/mol)
        self.pTime = 0.0     # 脉冲时间 (s)

        self.prob = 0.0      # 分子到达概率（论文公式）

        # 时间列表（初始化为 +∞）
        self.arrivalList: List[List[float]] = [
            [INFINITE_TIME for _ in range(simusize)] for _ in range(simusize)
        ]

    # ===================================================================
    # 设置参数（论文接口 SetParam）
    # ===================================================================
    def SetParam(self, temp: float, pressure: float, cellArea: float,
                 moleWeight: float, pulseTime: float):
        """
        temp: 温度 (K)
        pressure: 压强 (Pa)
        cellArea: 单元面积 (m^2)
        moleWeight: 分子量 (kg/mol)
        pulseTime: 脉冲时间 (s)
        """
        self.T = temp
        self.P = pressure
        self.S = cellArea
        self.M = moleWeight
        self.pTime = pulseTime

        self.CalcuProb()

    # ===================================================================
    # 根据气体动理论计算分子到达概率（论文 CalcuProb）
    # ===================================================================
    def CalcuProb(self):
        """
        论文基于 Maxwell-Boltzmann 气体分子碰撞率得：

        分子到达通量： J = P / sqrt(2π M R T)

        单元面积 S，则单位时间到达期望数：
            λ = J * S

        KMC 采用指数分布生成到达时间：
            t = -ln(r) / λ
        """

        R = 8.314  # 气体常数 J/(mol·K)

        if self.P <= 0 or self.T <= 0 or self.M <= 0:
            self.prob = 0
            return

        # 分子到达通量 J
        J = self.P / math.sqrt(2 * math.pi * self.M * R * self.T)

        # 到达率 λ
        lam = J * self.S

        self.prob = lam  # 为方便，直接存 λ

    # ===================================================================
    # 初始化整个到达时间列表（论文 InitArrivalList）
    # ===================================================================
    """
    Arrival selection rules (according to ALD reaction mechanism):

    TMA arrival:
        allowed only on –OH sites (Si–OH or Al–OH)

    H2O arrival:
        allowed only on Al–CH3 sites

    Reaction steps (A2, B6) are handled in react_event.py
    """
    def InitArrivalList(self,layer=None,gas=None):
        """
        初始化每个格点的到达事件时间：
            t = -ln(1 - r) / λ
        """
        if self.prob <= 0:
            # 到达概率为 0 => 永不发生
            for i in range(self.SIMUSIZE):
                for j in range(self.SIMUSIZE):
                    self.arrivalList[i][j] = INFINITE_TIME
            return

        for i in range(self.SIMUSIZE):
            for j in range(self.SIMUSIZE):
                # ----------- 过滤不可吸附位置（关键） ----------
                if layer is not None and gas is not None:
                    top = layer.cellTop[i][j]

                    # TMA 只能吸附在 OH 位点 (posA == 2)
                    # if gas == "TMA":
                    #     if top.posA != 2:
                    #         self.arrivalList[i][j] = INFINITE_TIME
                    #         continue
                    if gas == "TMA":
                        # 只能吸附在 –OH 位点（Si–OH 或 Al–OH）
                        if top.posA != 2:
                            self.arrivalList[i][j] = INFINITE_TIME
                            continue

                        # 如果未来区分基底/生长层，可在此加活化能
                        # if top.is_substrate:
                        #     pass

                    # H2O 只能吸附在 AlCH3 位点 (posA == 1)
                    if gas == "H2O":
                        if top.posA != 1:
                            self.arrivalList[i][j] = INFINITE_TIME
                            continue

                # ----------- 生成到达时间 -----------
                r = random.random()
                t = -math.log(1 - r) / self.prob
                self.arrivalList[i][j] = t

    # ===================================================================
    # 更新某一格点的到达事件时间（论文 UpdateArrivalList）
    # ===================================================================
    def UpdateArrivalList(self, i: int, j: int) -> float:
        """
        当一个事件发生后，该位置的事件时间需要重新采样。
        返回新采样的时间（相对时间间隔）
        """
        if self.prob <= 0:
            self.arrivalList[i][j] = INFINITE_TIME
            return INFINITE_TIME

        r = random.random()
        t = -math.log(1 - r) / self.prob
        self.arrivalList[i][j] = t
        return t

    # ===================================================================
    # 获取到达时间
    # ===================================================================
    def GetTime(self, i: int, j: int) -> float:
        return self.arrivalList[i][j]

    # ===================================================================
    # 设置某位置为永不发生（用于 KMC 的过滤机制）
    # ===================================================================
    def SetInfinite(self, i: int, j: int):
        self.arrivalList[i][j] = INFINITE_TIME