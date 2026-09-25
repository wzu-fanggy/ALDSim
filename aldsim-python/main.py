# main.py
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['SimHei']
matplotlib.rcParams['axes.unicode_minus'] = False
import os
import sys
from PyQt5.QtWidgets import QApplication

from ui.main_window import MainWindow


def ensure_directories():
    """
    创建项目所需的目录结构：
        data/
        data/export/
        data/models/
    """
    dirs = [
        "data",
        "data/export",
        "data/models",
    ]

    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
            print(f"[INIT] 创建目录: {d}")


def main():
    # --------------------------------------------------
    # 初始化目录
    # --------------------------------------------------
    ensure_directories()

    # --------------------------------------------------
    # 启动 Qt 主应用
    # --------------------------------------------------
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()