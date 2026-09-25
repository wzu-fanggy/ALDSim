# ui/dialog_data.py
import sys
import sqlite3
import csv
import os

from PyQt5.QtWidgets import (
    QDialog, QTableWidget, QTableWidgetItem, QVBoxLayout, QHBoxLayout,
    QPushButton, QMessageBox, QFileDialog, QApplication, QTextEdit, QLabel
)
from PyQt5.QtCore import Qt


DB_PATH = os.path.join("data", "ald.db")


class DataDialog(QDialog):
    """
    仿真数据管理对话框（论文 图 5-2）
    功能：
        - 浏览所有仿真记录
        - 查看某条记录详情
        - 删除记录
        - 导出数据
    """

    def __init__(self,parent=None):
        super().__init__()
        self.setWindowTitle("仿真数据管理")

        self.conn = sqlite3.connect(DB_PATH)
        self.cursor = self.conn.cursor()

        # 表格
        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(["ID", "时间", "循环数", "SIMUSIZE", "备注"])
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)

        # 按钮
        btn_view = QPushButton("查看详情")
        btn_delete = QPushButton("删除记录")
        btn_export = QPushButton("导出数据")
        btn_close = QPushButton("关闭")

        btn_view.clicked.connect(self.view_detail)
        btn_delete.clicked.connect(self.delete_record)
        btn_export.clicked.connect(self.export_data)
        btn_close.clicked.connect(self.close)

        btn_layout = QHBoxLayout()
        btn_layout.addWidget(btn_view)
        btn_layout.addWidget(btn_delete)
        btn_layout.addWidget(btn_export)
        btn_layout.addStretch()
        btn_layout.addWidget(btn_close)

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        layout.addLayout(btn_layout)
        self.setLayout(layout)

        self.load_data()

    # ===============================================================
    # 载入数据库记录
    # ===============================================================
    def load_data(self):
        self.table.setRowCount(0)

        try:
            self.cursor.execute(
                "SELECT id, datetime, cycle, simusize, note FROM simulation ORDER BY id DESC"
            )
            rows = self.cursor.fetchall()
        except Exception as e:
            QMessageBox.critical(self, "数据库错误", str(e))
            return

        for row_idx, row in enumerate(rows):
            self.table.insertRow(row_idx)
            for col_idx, val in enumerate(row):
                item = QTableWidgetItem(str(val))
                item.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(row_idx, col_idx, item)

    # ===============================================================
    # 查看记录详情
    # ===============================================================
    def view_detail(self):
        row = self.table.currentRow()
        if row < 0:
            QMessageBox.warning(self, "提示", "请先选择记录。")
            return

        rec_id = int(self.table.item(row, 0).text())

        # 查询详情
        try:
            self.cursor.execute(
                "SELECT detail FROM simulation WHERE id=?", (rec_id,)
            )
            result = self.cursor.fetchone()
        except Exception as e:
            QMessageBox.critical(self, "数据库错误", str(e))
            return

        detail_text = result[0] if result else ""

        dlg = QDialog(self)
        dlg.setWindowTitle(f"记录 {rec_id} 详情")
        layout = QVBoxLayout()

        text = QTextEdit()
        text.setReadOnly(True)
        text.setText(detail_text)

        btn = QPushButton("关闭")
        btn.clicked.connect(dlg.close)

        layout.addWidget(QLabel("数据详情："))
        layout.addWidget(text)
        layout.addWidget(btn)
        dlg.setLayout(layout)
        dlg.resize(600, 500)
        dlg.exec_()

    # ===============================================================
    # 删除记录
    # ===============================================================
    def delete_record(self):
        row = self.table.currentRow()
        if row < 0:
            QMessageBox.warning(self, "提示", "请选择要删除的记录。")
            return

        rec_id = int(self.table.item(row, 0).text())

        if QMessageBox.question(
            self, "删除确认", f"确定删除记录 {rec_id}？删除后不可恢复。",
            QMessageBox.Yes | QMessageBox.No
        ) == QMessageBox.No:
            return

        try:
            self.cursor.execute("DELETE FROM simulation WHERE id=?", (rec_id,))
            self.conn.commit()
        except Exception as e:
            QMessageBox.critical(self, "数据库错误", str(e))
            return

        QMessageBox.information(self, "成功", "记录已删除。")
        self.load_data()

    # ===============================================================
    # 导出数据 CSV
    # ===============================================================
    def export_data(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "导出 CSV", "simulation.csv", "CSV Files (*.csv)"
        )
        if not path:
            return

        try:
            self.cursor.execute(
                "SELECT id, datetime, cycle, simusize, note, detail FROM simulation"
            )
            rows = self.cursor.fetchall()
        except Exception as e:
            QMessageBox.critical(self, "数据库错误", str(e))
            return

        with open(path, "w", newline='', encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["ID", "时间", "循环数", "SIMUSIZE", "备注", "详情"])
            writer.writerows(rows)

        QMessageBox.information(self, "成功", f"已导出到 {path}")


# =================================================================
# 测试该对话框（不影响主程序）
# =================================================================
if __name__ == "__main__":
    app = QApplication(sys.argv)

    if not os.path.exists("data"):
        os.mkdir("data")
    if not os.path.exists(DB_PATH):  # 创建测试数据库
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("""
            CREATE TABLE simulation (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                datetime TEXT,
                cycle INTEGER,
                simusize INTEGER,
                note TEXT,
                detail TEXT
            )
        """)
        conn.commit()
        conn.close()

    dlg = DataDialog()
    dlg.resize(800, 500)
    dlg.exec_()

    sys.exit(0)