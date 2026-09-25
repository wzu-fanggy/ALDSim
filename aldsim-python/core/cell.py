# core/cell.py
from __future__ import annotations
from typing import Optional


class C2DCell:
    """
    单个二维单元（2D Cell），对应论文中 C2DCell 结构。
    采用双向链表结构，用于构建垂直方向的“柱”状薄膜。
    """

    def __init__(
        self,
        posA: int = 0,
        posB: int = 0,
        posC: int = 0,
        nCycle: int = 0,
        currNum: int = 0
    ):
        # 三个占位（代表 A/B/C 三原子位置）
        self.posA: int = posA
        self.posB: int = posB
        self.posC: int = posC

        # 本层的生成循环次数（用于记录反应周期）
        self.nCycle: int = nCycle

        # 当前层的层号（0 为最底层）
        self.currNum: int = currNum

        # 上一层（更高层）的指针
        self.next: Optional["C2DCell"] = None

        # 下一层（更低层）的指针
        self.previous: Optional["C2DCell"] = None

    def __repr__(self):
        return (
            f"C2DCell(A={self.posA}, B={self.posB}, C={self.posC}, "
            f"k={self.currNum}, cyc={self.nCycle})"
        )