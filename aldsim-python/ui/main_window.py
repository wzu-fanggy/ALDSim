# ui/main_window.py
import sys
from PyQt5.QtWidgets import (
    QMainWindow, QApplication, QAction, QMessageBox, QWidget, QStackedLayout
)
from PyQt5.QtCore import Qt

# 引入 UI 模块
from ui.dialog_params import ParamDialog
from ui.dialog_data import DataDialog
from ui.view_2d import View2D
from ui.view_3d import View3D

from core.config import SimuConfig

# 仿真处理模块
from core.simulation_proc import CSimulationProc

class MainWindow(QMainWindow):
    """ALD 模拟系统主窗口（论文图 5-1 总界面）"""

    def __init__(self):
        super().__init__()

        self.config = SimuConfig()

        self.setWindowTitle("ALD 模拟系统（基于 KMC 的原子层沉积仿真）")
        self.resize(900, 650)

        # ----------------------------------------------------------
        # 仿真模块
        # ----------------------------------------------------------
        self.sim = CSimulationProc(self.config)

        # ----------------------------------------------------------
        # Central widget: stacked views
        # ----------------------------------------------------------
        # central = QWidget()
        # self.stack = QStackedLayout()
        self.view2d = View2D()
        self.view3d = View3D(self.sim)

        self.stack = QStackedLayout()
        self.stack.addWidget(self.view2d)
        self.stack.addWidget(self.view3d)

        central = QWidget()
        central.setLayout(self.stack)
        self.setCentralWidget(central)

        # ----------------------------------------------------------
        # 菜单栏
        # ----------------------------------------------------------
        self.create_menu()

    # ==============================================================
    # 菜单栏结构（与论文一致）
    # ==============================================================
    def create_menu(self):
        menubar = self.menuBar()

        # ------------- 仿真菜单 -------------
        menu_sim = menubar.addMenu("仿真")

        act_params = QAction("参数设置", self)
        act_params.triggered.connect(self.open_params)
        menu_sim.addAction(act_params)

        act_run = QAction("开始仿真", self)
        act_run.triggered.connect(self.run_simulation)
        menu_sim.addAction(act_run)

        act_data = QAction("原始数据选择", self)
        act_data.triggered.connect(self.open_data_dialog)
        menu_sim.addAction(act_data)

        # ------------- 视图菜单 -------------
        menu_view = menubar.addMenu("视图")

        act_view2d = QAction("二维曲线图", self)
        act_view2d.triggered.connect(self.show_view_2d)
        menu_view.addAction(act_view2d)

        act_view3d = QAction("三维显示", self)
        act_view3d.triggered.connect(self.show_view_3d)
        menu_view.addAction(act_view3d)

    # ==============================================================
    # 打开参数设置窗口 (图 5-1)
    # ==============================================================
    def open_params(self):
        dlg = ParamDialog(self.config)
        if dlg.exec_():
            params = dlg.get_params()
            # 配置 simulation_proc 参数
            self.sim.apply_config(params)
            QMessageBox.information(self, "参数设置", "参数已成功保存。")

    # ==============================================================
    # 打开原始数据管理窗口 (图 5-2)
    # ==============================================================
    def open_data_dialog(self):
        dlg = DataDialog(self)
        dlg.exec_()

    # ==============================================================
    # 开始仿真
    # ==============================================================
    def run_simulation(self):
        reply = QMessageBox.question(
            self,
            "开始仿真",
            "确定要开始仿真吗？注意：仿真可能耗时较长。",
            QMessageBox.Yes | QMessageBox.No,
        )

        if reply == QMessageBox.No:
            return

        try:
            # 运行仿真，不需要返回值
            self.sim.run()
        except Exception as e:
            QMessageBox.critical(self, "仿真错误", f"仿真过程中出现错误：\n{e}")
            return

        # 正确地获取 2D 数据
        cycles, avg, rough, cover, rate = self.sim.get_2d_data()

        # 更新曲线
        self.view2d.set_data(cycles, avg, rough, cover, rate)

        QMessageBox.information(self, "仿真完成", "仿真已完成")# ==============================================================
    # 切换到二维视图
    # ==============================================================
    def show_view_2d(self):
        self.stack.setCurrentWidget(self.view2d)

    # ==============================================================
    # 切换到三维视图
    # ==============================================================
    # def show_view_3d(self):
    #     self.stack.setCurrentWidget(self.view3d)
    #
    #     self.view3d.update_from_sim()
    #
    #     # --- 安全解包 ---
    #     res = self.sim.export_3d_points()
    #     if res is None or len(res) < 2:
    #         QMessageBox.warning(self, "错误", "仿真数据为空或未初始化")
    #         return
    #
    #     pts, colors = res[0], res[1]
    #
    #     print("pts length =", len(pts))
    #     print("colors length =", len(colors))
    #     print("first few =", pts[:10])
    #
    #     self.view3d.set_points(pts, colors=colors)
    def show_view_3d(self):
        self.stack.setCurrentWidget(self.view3d)

        try:
            atoms = self.sim.export_3d_atoms()
            if not atoms:
                QMessageBox.warning(self, "错误", "三维原子数据为空")
                return

            self.view3d.set_atoms(atoms)

        except Exception as e:
            QMessageBox.critical(self, "三维显示错误", str(e))


# ==============================================================
# 单独调试 main_window
# ==============================================================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())