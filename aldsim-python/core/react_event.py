# core/react_event.py
from __future__ import annotations
import math
import random
from typing import List

INFINITE_TIME = 1e30
# Reaction types (according to mechanism)
REACTION_NONE = 0
REACTION_A2 = 1   # TMA + OH → AlCH3
REACTION_B6 = 2   # H2O + AlCH3 → OH

class CReactEventList:
    """
    反应事件时间表（对应论文 CReactEventList）

    每个 (i,j) 位置记录一个反应事件时间：
        t = - ln(r) / k
    其中 k = A * exp(-Ea / (R * T))

    论文假设可用 Arrhenius 形式计算反应概率或速率，并用于 KMC。
    """

    def __init__(self, simusize: int = 100):
        self.SIMUSIZE = simusize

        # 参数
        self.M = 0.0        # 分子量 (kg/mol) —— 可用于扩展模型
        # self.Ea = 0.0       # 反应激活能 (J/mol)
        self.Ea_A2 = 0.0  # TMA 反应
        self.Ea_B6 = 0.0  # H2O 反应
        self.T = 0.0        # 温度 (K)
        self.pTime = 0.0    # 脉冲时间 (s)
        self.prob = 0.0     # 反应事件概率（反应速率常数 k）

        # 时间表数组
        self.eventList: List[List[float]] = [
            [INFINITE_TIME for _ in range(simusize)] for _ in range(simusize)
        ]

    # ===================================================================
    # 设置参数（论文接口）
    # ===================================================================
    # def SetParam(self, activeEnergy: float, temp: float,
    #              pulseTime: float, moleWeight: float):
    #     self.Ea = activeEnergy
    #     self.T = temp
    #     self.pTime = pulseTime
    #     self.M = moleWeight
    #
    #     self.CalcuProb()
    def SetParam(self, Ea_A2: float, Ea_B6: float,
                 temp: float, pulseTime: float):
        self.Ea_A2 = Ea_A2
        self.Ea_B6 = Ea_B6
        self.T = temp
        self.pTime = pulseTime

    # ===================================================================
    # 计算反应概率（论文 CalcuProb）
    # ===================================================================
    def CalcuProb(self):
        """
        论文采用 Arrhenius 形式：

            k = A * exp(-Ea / (R*T))

        注意：
        - A 为指前因子，论文未明确给出，通常取 10^13 s⁻¹（原子尺度反应的典型值）
        - 若需要可作为参数加入 config 系统中
        """
        R = 8.314
        A = 1e13   # 典型 ALD 原子反应的指前因子

        if self.T <= 0:
            self.prob = 0
            return

        # self.prob = A * math.exp(-self.Ea / (R * self.T))

    def _calc_rate(self, Ea):
        R = 8.314
        A = 1e13
        return A * math.exp(-Ea / (R * self.T))

    # ===================================================================
    # 初始化所有格点的反应事件时间（加入化学过滤）
    # ===================================================================
    def InitEventList(self, layer=None, gas=None):
        """
        根据当前 layer 状态与气体类型（TMA 或 H2O），初始化反应事件表。

        反应条件：
        - TMA：需要 OH + TMA 已吸附  → top.posA == 2 AND top.posB == 1
        - H2O：需要 AlCH3 + H2O 已吸附 → top.posA == 1 AND top.posC == 1
        """

        # 速率 k <= 0：永不反应
        # if self.prob <= 0:
        #     for i in range(self.SIMUSIZE):
        #         for j in range(self.SIMUSIZE):
        #             self.eventList[i][j] = INFINITE_TIME
        #     return
        for i in range(self.SIMUSIZE):
            for j in range(self.SIMUSIZE):
                self.eventList[i][j] = INFINITE_TIME

        for i in range(self.SIMUSIZE):
            for j in range(self.SIMUSIZE):

                # ----------- 过滤不可反应位置 ----------
                if layer is not None and gas is not None:
                    top = layer.cellTop[i][j]

                    reaction_type = REACTION_NONE

                    # if gas == "TMA":
                    #     # TMA 反应：OH + (TMA 已吸附)
                    #     if not (top.posA == 2 and top.posB == 1):
                    #         self.eventList[i][j] = INFINITE_TIME
                    #         continue
                    if gas == "TMA":
                        # A2: OH + adsorbed TMA
                        if top.posA == 2 and top.posB == 1:
                            reaction_type = REACTION_A2

                    # elif gas == "H2O":
                    #     # H2O 反应：AlCH3 + (H2O 已吸附)
                    #     if not (top.posA == 1 and top.posC == 1):
                    #         self.eventList[i][j] = INFINITE_TIME
                    #         continue
                    elif gas == "H2O":
                        # B6: AlCH3 + adsorbed H2O
                        if top.posA == 1 and top.posC == 1:
                            reaction_type = REACTION_B6

                    if reaction_type == REACTION_NONE:
                        self.eventList[i][j] = INFINITE_TIME
                        continue

                    if reaction_type == REACTION_A2:
                        k = self._calc_rate(self.Ea_A2)
                    elif reaction_type == REACTION_B6:
                        k = self._calc_rate(self.Ea_B6)

                    r = random.random()
                    t = -math.log(1 - r) / k
                    self.eventList[i][j] = t

                # # ----------- 生成反应时间 -----------
                # r = random.random()
                # t = -math.log(1 - r) / self.prob
                # self.eventList[i][j] = t

    # ===================================================================
    # 更新单个格点反应事件时间（根据新状态重新判断是否可反应）
    # ===================================================================
    # def UpdateEventList(self, i: int, j: int, layer=None, gas=None):
    #     """
    #     当 (i,j) 发生反应或吸附后，需要根据新的表面状态重新生成事件时间。
    #
    #     若不满足反应条件，则设为 ∞。
    #     """
    #
    #     # k <= 0：永不反应
    #     if self.prob <= 0:
    #         self.eventList[i][j] = INFINITE_TIME
    #         return INFINITE_TIME
    #
    #     # ----------- 根据化学状态过滤 ----------
    #     if layer is not None and gas is not None:
    #         top = layer.cellTop[i][j]
    #
    #         if gas == "TMA":
    #             if not (top.posA == 2 and top.posB == 1):
    #                 self.eventList[i][j] = INFINITE_TIME
    #                 return INFINITE_TIME
    #
    #         elif gas == "H2O":
    #             if not (top.posA == 1 and top.posC == 1):
    #                 self.eventList[i][j] = INFINITE_TIME
    #                 return INFINITE_TIME
    #
    #     # ----------- 状态满足 → 生成新反应时间 ----------
    #     r = random.random()
    #     t = -math.log(1 - r) / self.prob
    #     self.eventList[i][j] = t
    #     return t
    def UpdateEventList(self, i: int, j: int, layer=None, gas=None):

        # 默认：设为 ∞
        self.eventList[i][j] = INFINITE_TIME

        if layer is None or gas is None:
            return INFINITE_TIME

        top = layer.cellTop[i][j]

        # ----------- 判断是否满足反应条件 ----------
        if gas == "TMA":
            if not (top.posA == 2 and top.posB == 1):
                return INFINITE_TIME
            k = self._calc_rate(self.Ea_A2)

        elif gas == "H2O":
            if not (top.posA == 1 and top.posC == 1):
                return INFINITE_TIME
            k = self._calc_rate(self.Ea_B6)

        else:
            return INFINITE_TIME

        # ----------- 生成新反应时间 ----------
        r = random.random()
        t = -math.log(1 - r) / k
        self.eventList[i][j] = t
        return t

    # ===================================================================
    # 设置某事件为 ∞（用于反应条件不满足时过滤）
    # ===================================================================
    def SetInfinite(self, i: int, j: int):
        self.eventList[i][j] = INFINITE_TIME

    # ===================================================================
    # 获取反应时间
    # ===================================================================
    def GetTime(self, i: int, j: int) -> float:
        return self.eventList[i][j]

    # ===================================================================
    # EventHappen：用于进一步处理反应结果（由 SimulationProc 调用）
    # ===================================================================
    def EventHappen(self, i: int, j: int):
        """
        论文接口：反应发生时需要重新生成下一个事件时间。
        （具体反应如何影响薄膜结构由 SimulationProc 决定）
        """
        self.UpdateEventList(i, j)
