# ui/view_3d.py
# ================== 原子可视化常量（论文级） ==================
ATOM_RADIUS = {   # 共价半径，单位 Å
    "H": 0.31,
    "C": 0.76,
    "O": 0.66,
    "Al": 1.21,
    "Si": 1.11,
}

ATOM_COLOR = {
    "H": (255, 255, 255, 255),     # white
    "C": (0, 0, 0, 255),           # black
    "O": (220, 50, 50, 255),        # red
    "Al": (180, 180, 180, 255),     # silver
    "Si": (40, 60, 160, 255),       # dark blue
}

RADIUS_SCALE = 0.35  # 视觉缩放因子（Å → OpenGL 空间单位）
import numpy as np
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel
from PyQt5.QtCore import Qt
import pyqtgraph as pg
import pyqtgraph.opengl as gl


def _color_to_rgba(color_name):
    """简单颜色名映射到 RGBA (0-255)。可根据需要扩展。"""
    cmap = {
        "blue": (50, 130, 255, 255),
        "green": (50, 200, 50, 255),
        "red": (220, 50, 50, 255),
        "yellow": (230, 220, 50, 255),
        "gray": (150, 150, 150, 255),
        "white": (255, 255, 255, 255),
        "black": (0, 0, 0, 255),
    }
    if isinstance(color_name, (tuple, list)) and len(color_name) in (3, 4):
        # already rgb or rgba (0-255)
        if len(color_name) == 3:
            return (int(color_name[0]), int(color_name[1]), int(color_name[2]), 255)
        return tuple(int(x) for x in color_name)
    return cmap.get(color_name, cmap["blue"])


class View3D(QWidget):
    """
    PyQtGraph + OpenGL 版本的 3D 点云视图（稳定、性能好）
    使用方法：
        v = View3D(sim)                # sim 可为 None
        v.set_points(pts, colors, sizes)
    pts: iterable of (x,y,z)
    colors: iterable of color names or rgb tuples (len == len(pts)) 或 None
    sizes: iterable of numbers or None
    """

    def __init__(self, sim=None):
        super().__init__()
        self.sim = sim
        self.points = np.empty((0, 3), dtype=float)
        self.colors = None
        self.sizes = None
        self._init_ui()
        self.atom_items = []  # 保存所有球，方便清空

    # def _init_ui(self):
    #     self.setWindowTitle("三维结构显示（PyQtGraph）")
    #     layout = QVBoxLayout()
    #     self.setLayout(layout)
    #
    #     # GL 视图
    #     self.view = gl.GLViewWidget()
    #     self.view.opts['distance'] = 40  # 默认摄像机距离，可调整
    #     layout.addWidget(self.view)
    #
    #     # 状态栏
    #     self.status = QLabel("Ready")
    #     self.status.setAlignment(Qt.AlignLeft)
    #     layout.addWidget(self.status)
    #
    #     # 添加网格参考
    #     grid = gl.GLGridItem()
    #     grid.scale(1, 1, 1)
    #     self.view.addItem(grid)
    #
    #     # Scatter item 占位（后面 set_points 会替换）
    #     self.scatter = None
    #
    #     # 允许鼠标交互（GLViewWidget 默认支持）
    #     # 绑定鼠标移动事件（可选）
    #     self.view.mouseMoveEvent = self._on_mouse_move_gl
    def _init_ui(self):
        self.setWindowTitle("三维结构显示（PyQtGraph）")

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)  # 去掉四周空白
        layout.setSpacing(0)  # 去掉控件间距
        self.setLayout(layout)

        # ===== OpenGL 视图 =====
        self.view = gl.GLViewWidget()
        self.view.opts['distance'] = 40
        self.view.setBackgroundColor((235, 235, 235))  # ★新增：背景色★
        self.view.setSizePolicy(
            pg.QtWidgets.QSizePolicy.Expanding,
            pg.QtWidgets.QSizePolicy.Expanding
        )
        layout.addWidget(self.view, stretch=1)  # ★关键：拉满★

        # ===== 状态栏 =====
        self.status = QLabel("Ready")
        self.status.setFixedHeight(22)  # ★固定高度★
        self.status.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        layout.addWidget(self.status, stretch=0)

        # ===== 网格 =====
        grid = gl.GLGridItem()
        grid.scale(1, 1, 1)
        self.view.addItem(grid)

        self.scatter = None
        self.view.mouseMoveEvent = self._on_mouse_move_gl

    def _on_mouse_move_gl(self, ev):
        # 显示鼠标在窗口中的位置（屏幕坐标）
        pos = ev.pos()
        self.status.setText(f"x={pos.x()} y={pos.y()}")
        # 保持默认行为
        super(gl.GLViewWidget, self.view).mouseMoveEvent(ev)

    def set_atoms(self, atoms):
        """
        atoms: [(x, y, z, element), ...]
        """
        self._clear_atoms()

        for (x, y, z, elem) in atoms:
            item = self._add_atom_sphere(x, y, z, elem)
            self.atom_items.append(item)

    # 兼容旧接口：直接从 simulation_proc 导出并显示
    def update_from_sim(self):
        if self.sim is None:
            return
        # 要求 simulation_proc 提供 export_3d_points() 返回 pts, colors, sizes(opt)
        try:
            res = self.sim.export_3d_points()
            if res is None:
                return
            if len(res) == 2:
                pts, colors = res
                sizes = None
            else:
                pts, colors, sizes = res
            self.set_points(pts, colors=colors, sizes=sizes)
        except Exception as e:
            # 不抛出致命错误，显示在状态栏
            self.status.setText(f"更新失败: {e}")

    def _add_atom_sphere(self, x, y, z, element):
        """
        在 (x,y,z) 位置添加一个原子球
        """
        radius = ATOM_RADIUS[element] * RADIUS_SCALE
        color = np.array(ATOM_COLOR[element], dtype=np.ubyte) / 255.0

        md = gl.MeshData.sphere(rows=16, cols=16, radius=radius)
        mesh = gl.GLMeshItem(
            meshdata=md,
            smooth=True,
            color=color,
            shader="shaded",
            drawEdges=False
        )
        mesh.translate(x, y, z)
        self.view.addItem(mesh)

        return mesh

    def _clear_atoms(self):
        for item in self.atom_items:
            self.view.removeItem(item)
        self.atom_items.clear()

    def set_points(self, pts, colors=None, sizes=None):
        """
        兼容旧接口：将 (x,y,z) 点云临时映射为 Si 原子
        """
        atoms = []
        for (x, y, z) in pts:
            atoms.append((x, y, z, "Si"))  # 临时全部当作 Si

        self.set_atoms(atoms)