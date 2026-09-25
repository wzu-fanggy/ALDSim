# ui/view_2d.py
import sys
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QCheckBox, QDialog,
    QPushButton, QLabel, QFormLayout
)
from PyQt5.QtCore import Qt
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class CurveSelectionDialog(QDialog):
    """曲线选择对话框（双击主图弹出）"""

    def __init__(self, flags: dict):
        super().__init__()
        self.setWindowTitle("选择显示曲线")
        self.flags = flags

        self.chk_avg = QCheckBox("平均厚度")
        self.chk_rough = QCheckBox("表面粗糙度")
        self.chk_cover = QCheckBox("覆盖率")
        self.chk_rate = QCheckBox("沉积速率")

        self.chk_avg.setChecked(flags["avg"])
        self.chk_rough.setChecked(flags["rough"])
        self.chk_cover.setChecked(flags["cover"])
        self.chk_rate.setChecked(flags["rate"])

        layout = QFormLayout()
        layout.addRow(self.chk_avg)
        layout.addRow(self.chk_rough)
        layout.addRow(self.chk_cover)
        layout.addRow(self.chk_rate)

        btn_ok = QPushButton("确定")
        btn_ok.clicked.connect(self.accept)

        root = QVBoxLayout()
        root.addLayout(layout)
        root.addWidget(btn_ok)

        self.setLayout(root)
        self.resize(250, 200)

    def get_flags(self):
        return {
            "avg": self.chk_avg.isChecked(),
            "rough": self.chk_rough.isChecked(),
            "cover": self.chk_cover.isChecked(),
            "rate": self.chk_rate.isChecked(),
        }


class View2D(QWidget):
    """
    二维图显示界面（论文 图 5-3）
    显示四条曲线：
        - 平均厚度
        - 表面粗糙度
        - 覆盖率
        - 沉积速率
    """

    def __init__(self):
        super().__init__()

        self.setWindowTitle("二维曲线图")

        # 曲线显示开关
        self.show_flags = {
            "avg": True,
            "rough": True,
            "cover": True,
            "rate": True
        }

        # Matplotlib 图表
        self.fig = Figure(figsize=(6, 5))
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.ax = self.fig.add_subplot(111)

        self.status = QLabel("X=0, Y=0")
        self.status.setAlignment(Qt.AlignLeft)

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.status)
        self.setLayout(layout)

        # 鼠标事件：显示坐标
        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.canvas.mpl_connect("button_press_event", self.on_double_click)

        # 空数据初始化
        self.cycles = []
        self.avg_thick = []
        self.roughness = []
        self.coverage = []
        self.rate = []

        self.draw_plot()

    # ======================================================================
    # 设置数据接口 —— simulation_proc 调用
    # ======================================================================
    def set_data(self, cycles, avg, rough, cover, rate):
        self.cycles = cycles
        self.avg_thick = avg
        self.roughness = rough
        self.coverage = cover
        self.rate = rate
        self.draw_plot()

    # ======================================================================
    # 鼠标移动显示坐标
    # ======================================================================
    def on_mouse_move(self, event):
        if event.xdata is None or event.ydata is None:
            return
        self.status.setText(f"X={event.xdata:.2f},  Y={event.ydata:.2f}")

    # ======================================================================
    # 双击弹出曲线选择窗口
    # ======================================================================
    def on_double_click(self, event):
        if event.dblclick:
            dlg = CurveSelectionDialog(self.show_flags)
            if dlg.exec_():
                self.show_flags = dlg.get_flags()
                self.draw_plot()

    # ======================================================================
    # 绘图
    # ======================================================================
    def draw_plot(self):
        self.ax.clear()
        self.ax.set_title("ALD 仿真二维曲线")
        self.ax.set_xlabel("反应循环次数")
        self.ax.set_ylabel("值")

        if len(self.cycles) == 0:
            self.canvas.draw()
            return

        # 绘制各曲线
        if self.show_flags["avg"]:
            self.ax.plot(self.cycles, self.avg_thick, label="平均厚度")

        if self.show_flags["rough"]:
            self.ax.plot(self.cycles, self.roughness, label="表面粗糙度")

        if self.show_flags["cover"]:
            self.ax.plot(self.cycles, self.coverage, label="覆盖率")

        if self.show_flags["rate"]:
            self.ax.plot(self.cycles, self.rate, label="沉积速率")

        # 显示图例
        self.ax.legend()
        self.ax.grid(True)

        self.canvas.draw()


# ============================================================
# 单独测试
# ============================================================
if __name__ == "__main__":
    from PyQt5.QtWidgets import QApplication
    import numpy as np

    app = QApplication(sys.argv)
    w = View2D()

    # 测试数据
    x = list(range(1, 51))
    w.set_data(
        x,
        [i * 0.12 for i in x],
        [0.5 + 0.1 * (i % 5) for i in x],
        [1 - 0.8 * (1 / (i + 1)) for i in x],
        [0.1 + 0.01 * i for i in x],
    )

    w.resize(700, 550)
    w.show()
    sys.exit(app.exec_())