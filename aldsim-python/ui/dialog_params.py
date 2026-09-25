# ui/dialog_params.py
import sys
from PyQt5.QtWidgets import (
    QDialog, QLabel, QLineEdit, QPushButton, QGridLayout,
    QHBoxLayout, QVBoxLayout, QApplication, QMessageBox
)
from core.config import SimuConfig


class ParamDialog(QDialog):
    """
    ALD 参数设置对话框
    对应论文 图 5-1 “参数设置”窗口
    """

    def __init__(self, config: SimuConfig):
        super().__init__()

        self.setWindowTitle("仿真参数设置")
        self.config = config

        # =============== 创建控件 ===============
        # TMA 参数
        self.edit_TMA_T = QLineEdit(str(config.T_TMA))
        self.edit_TMA_P = QLineEdit(str(config.P_TMA))
        self.edit_TMA_t = QLineEdit(str(config.t_TMA))

        # H2O 参数
        self.edit_H2O_T = QLineEdit(str(config.T_H2O))
        self.edit_H2O_P = QLineEdit(str(config.P_H2O))
        self.edit_H2O_t = QLineEdit(str(config.t_H2O))

        # 清洗参数
        self.edit_clean_t = QLineEdit(str(config.t_clean))

        # 公共参数
        self.edit_cycle = QLineEdit(str(config.cycle_num))
        self.edit_SIMU = QLineEdit(str(config.SIMUSIZE))

        # =============== 布局 ===============
        layout = QVBoxLayout()

        grid = QGridLayout()
        row = 0

        # ---------- TMA ----------
        grid.addWidget(QLabel("TMA 脉冲温度 (K)"), row, 0)
        grid.addWidget(self.edit_TMA_T, row, 1)
        row += 1

        grid.addWidget(QLabel("TMA 脉冲压力 (Pa)"), row, 0)
        grid.addWidget(self.edit_TMA_P, row, 1)
        row += 1

        grid.addWidget(QLabel("TMA 脉冲时间 (s)"), row, 0)
        grid.addWidget(self.edit_TMA_t, row, 1)
        row += 1

        # ---------- H2O ----------
        grid.addWidget(QLabel("H₂O 脉冲温度 (K)"), row, 0)
        grid.addWidget(self.edit_H2O_T, row, 1)
        row += 1

        grid.addWidget(QLabel("H₂O 脉冲压力 (Pa)"), row, 0)
        grid.addWidget(self.edit_H2O_P, row, 1)
        row += 1

        grid.addWidget(QLabel("H₂O 脉冲时间 (s)"), row, 0)
        grid.addWidget(self.edit_H2O_t, row, 1)
        row += 1

        # ---------- Clean ----------
        grid.addWidget(QLabel("清洗气体脉冲时间 (s)"), row, 0)
        grid.addWidget(self.edit_clean_t, row, 1)
        row += 1

        # ---------- 公共 ----------
        grid.addWidget(QLabel("模拟循环次数"), row, 0)
        grid.addWidget(self.edit_cycle, row, 1)
        row += 1

        grid.addWidget(QLabel("二维模拟尺寸 (N×N)"), row, 0)
        grid.addWidget(self.edit_SIMU, row, 1)
        row += 1

        layout.addLayout(grid)

        # =============== 按钮区域 ===============
        btn_ok = QPushButton("确定")
        btn_cancel = QPushButton("取消")

        btn_ok.clicked.connect(self.save_params)
        btn_cancel.clicked.connect(self.reject)

        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(btn_ok)
        btn_layout.addWidget(btn_cancel)

        layout.addLayout(btn_layout)

        self.setLayout(layout)

    # ==============================================================
    # 保存参数
    # ==============================================================
    def save_params(self):
        try:
            self.config.T_TMA = float(self.edit_TMA_T.text())
            self.config.P_TMA = float(self.edit_TMA_P.text())
            self.config.t_TMA = float(self.edit_TMA_t.text())

            self.config.T_H2O = float(self.edit_H2O_T.text())
            self.config.P_H2O = float(self.edit_H2O_P.text())
            self.config.t_H2O = float(self.edit_H2O_t.text())

            self.config.t_clean = float(self.edit_clean_t.text())

            self.config.cycle_num = int(self.edit_cycle.text())
            self.config.SIMUSIZE = int(self.edit_SIMU.text())

        except ValueError:
            QMessageBox.warning(self, "输入错误", "请确保所有参数均为数字！")
            return

        self.accept()

    # ==============================================================
    # 获取参数
    # ==============================================================
    def get_params(self):
        """
        返回当前配置对象，如果输入有误则弹出提示并保持原配置
        """
        self.save_params()  # 尝试保存参数
        return self.config


# ======================================================================
# 单独测试窗口（不影响主程序）
# ======================================================================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    cfg = SimuConfig()
    dlg = ParamDialog(cfg)
    if dlg.exec_():
        print(cfg)
    sys.exit(0)