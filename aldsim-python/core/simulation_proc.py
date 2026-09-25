# core/simulation_proc.py
from __future__ import annotations
import time
import sqlite3
from typing import Tuple, Optional

from .layer import CLayer
from .react_arrival import CReactantArrival
from .react_event import CReactEventList
from .config import SimuConfig

INFINITE_TIME = 1e30


class CSimulationProc:
    """
    对应论文 CSimulationProc：
    - 控制 ALD 模拟总流程
    - 调用 Layer / Arrival / Event 各模块
    - 在每个气体脉冲阶段执行 KMC
    - 写入数据库
    """


    def __init__(self, config: SimuConfig, db_path="data/ald.db"):
        self.cfg = config
        self.db_path = db_path

        # 初始化薄膜结构
        self.layer = CLayer(simusize=config.SIMUSIZE)

        # TMA 到达 + 反应事件表
        self.arrival_TMA = CReactantArrival(config.SIMUSIZE)
        self.event_TMA = CReactEventList(config.SIMUSIZE)

        # H₂O 到达 + 反应事件表
        self.arrival_H2O = CReactantArrival(config.SIMUSIZE)
        self.event_H2O = CReactEventList(config.SIMUSIZE)

        # 数据库连接
        self._init_database()

        # ============================
        # 新增：统计数组
        # ============================
        self.stat_cycles = []
        self.stat_avg = []
        self.stat_rough = []
        self.stat_cover = []
        self.stat_rate = []

    # ==============================================================
    # 数据库初始化
    # ==============================================================
    def _init_database(self):
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()

        cur.execute("""
        CREATE TABLE IF NOT EXISTS simulation_raw (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sim_time REAL,
            cycle INTEGER,
            gas_phase TEXT,
            i INTEGER,
            j INTEGER,
            height INTEGER
        )
        """)

        conn.commit()
        conn.close()

    # ==============================================================
    # 将当前晶体结构写入数据库
    # ==============================================================
    def _save_layer_to_db(self, sim_time: float, cycle: int, gas: str):
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()

        for i in range(self.cfg.SIMUSIZE):
            for j in range(self.cfg.SIMUSIZE):
                h = self.layer.GetHeight(i, j)
                cur.execute("""
                INSERT INTO simulation_raw(sim_time, cycle, gas_phase, i, j, height)
                VALUES (?, ?, ?, ?, ?, ?)
                """, (sim_time, cycle, gas, i, j, h))

        conn.commit()
        conn.close()

    # ==============================================================
    # 初始化阶段
    # ==============================================================
    def InitSiLayer(self):
        """初始硅晶面构建"""
        self.layer.InitSiLayer()

    # ==============================================================
    # KMC 主控循环（论文图 4-12）
    # ==============================================================
    def run(self):
        sim_start = time.time()

        # 1. 初始化基底
        print("初始化基底层...")
        self.InitSiLayer()

        # 2. 初始化记录
        cycle = 0

        # ------------------------------------
        # 开始 ALD 循环
        # ------------------------------------
        for cycle in range(1, self.cfg.cycle_num + 1):
            print(f"=== 循环 {cycle} 开始 ===")

            # ===========================================
            # (1) TMA 气体脉冲
            # ===========================================
            print("TMA 脉冲阶段开始")
            t_end_TMA = self._run_single_gas(
                gas="TMA",
                arrival=self.arrival_TMA,
                event=self.event_TMA,
                temperature=self.cfg.T_TMA,
                pressure=self.cfg.P_TMA,
                pulse_time=self.cfg.t_TMA,
                mole_weight=0.072  # 可按模型需要补充
            )
            self._save_layer_to_db(t_end_TMA, cycle, "TMA")

            # ===========================================
            # (2) H2O 气体脉冲
            # ===========================================
            print("H₂O 脉冲阶段开始")
            t_end_H2O = self._run_single_gas(
                gas="H2O",
                arrival=self.arrival_H2O,
                event=self.event_H2O,
                temperature=self.cfg.T_H2O,
                pressure=self.cfg.P_H2O,
                pulse_time=self.cfg.t_H2O,
                mole_weight=0.018
            )
            self._save_layer_to_db(t_end_H2O, cycle, "H2O")

            # ===========================================
            # (3) 清洗阶段（无反应，仅等待）
            # ===========================================
            print("净化阶段 (Clean)")
            sim_time_clean = t_end_H2O + self.cfg.t_clean
            self._save_layer_to_db(sim_time_clean, cycle, "CLEAN")

            # -----------------------------------------------------
            # 新增：统计曲线数据
            # -----------------------------------------------------
            avg_h, rough, cover = self.calc_stats()

            self.stat_cycles.append(cycle)
            self.stat_avg.append(avg_h)
            self.stat_rough.append(rough)
            self.stat_cover.append(cover)

            # 沉积速率（平均厚度变化）
            if cycle == 1:
                self.stat_rate.append(avg_h)
            else:
                self.stat_rate.append(avg_h - self.stat_avg[-2])

        print("仿真结束！总用时 %.2f s" % (time.time() - sim_start))

    # ==============================================================
    # 应用配置参数
    # ==============================================================
    def apply_config(self, config: SimuConfig):
        """
        参数对话框更新后，必须重建所有依赖 SIMUSIZE 的结构，并重新初始化 Si 基底。
        """
        self.cfg = config

        # 重建所有与网格大小相关的对象
        self.layer = CLayer(simusize=config.SIMUSIZE)
        self.arrival_TMA = CReactantArrival(config.SIMUSIZE)
        self.event_TMA = CReactEventList(config.SIMUSIZE)
        self.arrival_H2O = CReactantArrival(config.SIMUSIZE)
        self.event_H2O = CReactEventList(config.SIMUSIZE)

        # 关键：必须初始化底层 Si 网格，否则 export_3d_points() 会收到空链表结构
        self.layer.InitSiLayer()


    # ==============================================================
    # 单个气体脉冲阶段 KMC
    # ==============================================================
    def _run_single_gas(
        self,
        gas: str,
        arrival: CReactantArrival,
        event: CReactEventList,
        temperature: float,
        pressure: float,
        pulse_time: float,
        mole_weight: float
    ) -> float:
        """
        KMC time evolution during a single precursor pulse.
        Arrival (A1/B1) and surface reaction (A2/B6) events
        are executed competitively based on their stochastic times.
        """

        # --------------------------------------------
        # 初始化参数（到达事件 + 反应事件）
        # --------------------------------------------
        arrival.SetParam(temp=temperature,
                         pressure=pressure,
                         cellArea=self.layer.cell_area,   # 必须提供单位面积
                         moleWeight=mole_weight,
                         pulseTime=pulse_time)
        arrival.InitArrivalList(layer=self.layer, gas=gas)

        # event.SetParam(activeEnergy=self.cfg.Ea,
        #                temp=temperature,
        #                pulseTime=pulse_time,
        #                moleWeight=mole_weight)
        event.SetParam(
            Ea_A2=self.cfg.Ea_TMA,  # TMA + OH → AlCH3
            Ea_B6=self.cfg.Ea_H2O,  # H2O + AlCH3 → OH
            temp=temperature,
            pulseTime=pulse_time
        )
        event.InitEventList(layer=self.layer, gas=gas)

        # --------------------------------------------
        # KMC 主循环：在 pulse_time 之前不断执行事件
        # --------------------------------------------
        current_time = 0.0
        loop_count = 0
        MAX_LOOP = 200000

        while True:
            loop_count += 1
            if loop_count > MAX_LOOP:
                print("⚠ KMC 死循环保护触发，强制退出")
                break
            # 找最小到达时间
            # t_arr, ia, ja = self._find_min(arrival.arrivalList)
            res_arr = self._find_min(arrival.arrivalList)
            if res_arr is None:
                print("ERROR: _find_min returned None for arrivalList -- aborting pulse")
                break
            t_arr, ia, ja = res_arr

            # 找最小反应时间
            # t_evt, ir, jr = self._find_min(event.eventList)
            res_evt = self._find_min(event.eventList)
            if res_evt is None:
                print("ERROR: _find_min returned None for eventList -- aborting pulse")
                break
            t_evt, ir, jr = res_evt

            # 如果 t_arr 或 t_evt 是 INFINITE_TIME，就直接退出循环
            if t_arr == INFINITE_TIME and t_evt == INFINITE_TIME:
                break

            # 在当前时间轴上加偏移
            t_arr += current_time
            t_evt += current_time

            # 若两个事件都超过脉冲时间 → 结束
            if t_arr > pulse_time and t_evt > pulse_time:
                break

        # ------------------------------------------------------
        # 执行最小事件
        # ------------------------------------------------------
        #     if t_arr <= t_evt and t_arr <= pulse_time:
        #         # 到达事件（吸附）
        #         current_time = t_arr
        #         self.layer.GasArrival(i=ia, j=ja, gas=gas)
        #         arrival.UpdateArrivalList(ia, ja)
        #
        #     elif t_evt <= pulse_time:
        #         # 反应事件（形成 Al-O 或 –OH 结构）
        #         current_time = t_evt
        #         # self.layer.React(i=ir, j=jr, gas=gas)
        #         # event.UpdateEventList(ir, jr,layer=self.layer, gas=gas)
        #         self.layer.React(i=ir, j=jr, gas=gas)
        #
        #         # 更新该点 arrival（是否还能吸附）
        #         arrival.UpdateArrivalList(ir, jr)
        #
        #         # arrival 之后，可能满足反应条件
        #         event.UpdateEventList(ia, ja, layer=self.layer, gas=gas)
        # return current_time
        if t_arr <= t_evt and t_arr <= pulse_time:
            # =========================
            # 到达事件 A1 / B1
            # =========================
            current_time = t_arr

            # 1️⃣ 执行吸附
            self.layer.GasArrival(i=ia, j=ja, gas=gas)

            # 2️⃣ 更新该点 arrival（下一次到达）
            arrival.UpdateArrivalList(ia, ja)

            # 3️⃣ arrival 之后，可能触发反应（A2 / B6）
            event.UpdateEventList(ia, ja, layer=self.layer, gas=gas)

        elif t_evt <= pulse_time:
            # =========================
            # 反应事件 A2 / B6
            # =========================
            current_time = t_evt

            # 1️⃣ 执行反应（可能生长一层）
            self.layer.React(i=ir, j=jr, gas=gas)

            # 2️⃣ 反应后，该点的 arrival 条件发生变化
            arrival.UpdateArrivalList(ir, jr)

            # 3️⃣ 反应后，该点的反应条件也发生变化
            event.UpdateEventList(ir, jr, layer=self.layer, gas=gas)

    # ==============================================================
    # 找最小事件时间（更健壮）
    # ==============================================================
    def _find_min(self, arr):
        try:
            min_v = INFINITE_TIME
            min_i = -1
            min_j = -1
            N = self.cfg.SIMUSIZE

            # 基本验证：arr 必须为二维列表且维度匹配
            if arr is None:
                print("ERROR: _find_min received arr is None")
                return INFINITE_TIME, -1, -1

            if not hasattr(arr, '__len__'):
                print("ERROR: _find_min received non-iterable arr:", type(arr))
                return INFINITE_TIME, -1, -1

            # 防护：若 arr 行数与 N 不一致，使用实际行数以避免索引错误
            rows = len(arr)
            cols = len(arr[0]) if rows > 0 else 0

            if rows < N or cols < N:
                # 报警但仍尝试使用最小维度，避免抛出索引异常
                # 这有助于甄别配置/初始化不一致的问题
                # 打印诊断信息
                print(f"WARNING: arr shape ({rows},{cols}) smaller than cfg.SIMUSIZE={N}. Using min dims.")
                N_r = min(rows, N)
                N_c = min(cols, N)
            else:
                N_r = N
                N_c = N

            for i in range(N_r):
                row = arr[i]
                for j in range(N_c):
                    v = row[j]
                    # 防护：若元素为 None 或 非数值，跳过并打印
                    if v is None:
                        # 不直接抛异常，记录并继续
                        # 这通常是导致之前异常的罪魁
                        print(f"DEBUG: _find_min encountered None at ({i},{j})")
                        continue
                    try:
                        if v < min_v:
                            min_v = v
                            min_i = i
                            min_j = j
                    except TypeError:
                        print(f"DEBUG: _find_min type error comparing value at ({i},{j}): {v} (type {type(v)})")
                        continue

            # return min_v, min_i, min_j
            # ========= 修复 1：强制合法返回 =========
            if min_i == -1 or min_j == -1:
                return INFINITE_TIME, -1, -1

            return float(min_v), min_i, min_j
        except Exception as e:
            # 捕获不可预期异常，提供诊断输出并保证返回值类型一致
            print("_find_min() unexpected exception:", repr(e))
            return INFINITE_TIME, -1, -1

    # # ==============================================================
    # # 导出三维点云（给 View3D 使用）
    # # ==============================================================
    # def export_3d_points(self):
    #     # 如果 layer 未初始化，直接返回空
    #     if self.layer is None:
    #         return [], []
    #
    #     pts = []
    #     colors = []
    #
    #     N = self.cfg.SIMUSIZE
    #
    #     try:
    #         for i in range(N):
    #             for j in range(N):
    #                 h = self.layer.GetHeight(i, j)
    #                 pts.append((i, j, h))
    #
    #                 if h <= 1:
    #                     colors.append("blue")
    #                 elif h <= 3:
    #                     colors.append("green")
    #                 else:
    #                     colors.append("red")
    #
    #         return pts, colors
    #
    #     except Exception as e:
    #         print("export_3d_points() failed:", e)
    #         return [], []

    def calc_stats(self):
        N = self.cfg.SIMUSIZE
        h_list = []

        for i in range(N):
            for j in range(N):
                h_list.append(self.layer.GetHeight(i, j))

        # 平均厚度
        avg_h = sum(h_list) / len(h_list)

        # 粗糙度（厚度标准差）
        rough = (sum((h - avg_h) ** 2 for h in h_list) / len(h_list)) ** 0.5

        # 覆盖率（OH 位置占比）
        cover_cnt = 0
        total = N * N
        for i in range(N):
            for j in range(N):
                if self.layer.cellTop[i][j].posA == 2:  # OH
                    cover_cnt += 1
        coverage = cover_cnt / total

        return avg_h, rough, coverage

    def get_2d_data(self):
        return (
            self.stat_cycles,
            self.stat_avg,
            self.stat_rough,
            self.stat_cover,
            self.stat_rate,
        )

    # ==============================================================
    # 导出三维原子结构（论文级，用于球模型）
    # ==============================================================
    def export_3d_atoms(self):
        """
        返回原子列表:
        (x, y, z, element)
        element ∈ {"Si", "Al", "O", "H"}
        """
        atoms = []
        N = self.cfg.SIMUSIZE

        for i in range(N):
            for j in range(N):
                cell = self.layer.cellTop[i][j]

                # ---- Si 基底 ----
                atoms.append((i, j, 0.0, "Si"))

                z = 0.0

                # ---- posA: 表面官能团 ----
                if cell.posA == 1:  # *OH
                    z += 0.8
                    atoms.append((i, j, z, "O"))
                    atoms.append((i, j, z + 0.4, "H"))

                elif cell.posA == 2:  # =O
                    z += 0.8
                    atoms.append((i, j, z, "O"))

                # ---- posB: Al ----
                if cell.posB == 1:
                    z += 1.2
                    atoms.append((i, j, z, "Al"))

                # ---- posC: 上层 O ----
                if cell.posC == 1:
                    z += 1.0
                    atoms.append((i, j, z, "O"))

        return atoms