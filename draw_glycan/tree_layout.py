"""Tree layout — clean-room Reingold–Tilford "tidy tree" implementation.

Strategy (equivalent to d3.tree with default separation for typical glycans):
    1. Assign each leaf an integer x slot in left-to-right DFS order
       (each leaf takes one slot, siblings are adjacent).
    2. Each internal node is centered over its children.
    3. After placing children, check that this subtree does not overlap with
       the previous sibling subtree; if it does, shift it right.

This produces the same left-to-right topology as d3.tree() for the
shallow, multi-branch trees typical of glycans.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional


class Node:
    __slots__ = ('data', 'depth', 'height', 'parent', 'children', 'x', 'y')

    def __init__(self, data: Dict[str, Any], depth: int = 0, parent: Optional['Node'] = None):
        self.data = data
        self.depth = depth
        self.parent = parent
        self.children: List[Node] = []
        self.x: float = 0.0
        self.y: float = 0.0
        self.height: int = 0

    def each(self, fn):
        fn(self)
        for c in self.children:
            c.each(fn)

    def descendants(self) -> List['Node']:
        out = [self]
        for c in self.children:
            out.extend(c.descendants())
        return out

    def leaves(self) -> List['Node']:
        if not self.children:
            return [self]
        out: List[Node] = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def hierarchy(data: Dict[str, Any]) -> Node:
    root = Node(data, depth=0, parent=None)
    _build(root)
    _assign_height(root)
    return root


def _build(node: Node):
    for child_data in node.data.get('children', []) or []:
        child = Node(child_data, depth=node.depth + 1, parent=node)
        node.children.append(child)
        _build(child)


def _assign_height(node: Node) -> int:
    if not node.children:
        node.height = 0
        return 0
    h = 0
    for c in node.children:
        ch = _assign_height(c) + 1
        if ch > h:
            h = ch
    node.height = h
    return h


def tree_layout(root: Node, size: tuple):
    """Tidy tree layout matching d3.tree().size([width, height]).

    After the call each node has:
        x in [0, width]  (breadth)
        y in [0, height] (depth * layer_height)
    """
    width, height = size
    if root is None:
        return

    if not root.leaves():
        return

    # 1) assign preliminary x: each leaf gets sequential integer slot in
    #    left-to-right DFS order; internal node = midpoint of children.
    counter = [0]
    _next_leaf_slot(root, counter)
    _assign_preliminary_x(root)

    # 2) resolve overlaps: walk through siblings and shift right subtrees
    #    that would overlap the right contour of the left neighbour.
    _resolve_overlaps(root)

    # 3) normalize to [0, width] over leaves, [0, height] over depth
    leaves = root.leaves()
    x0 = min(n.x for n in leaves)
    x1 = max(n.x for n in leaves)
    kx = (width / (x1 - x0)) if x1 > x0 else 0.0
    ky = (height / root.height) if root.height > 0 else 0.0

    for n in root.descendants():
        n.x = (n.x - x0) * kx
        n.y = n.depth * ky


def _assign_preliminary_x(node: Node):
    if not node.children:
        return
    for c in node.children:
        _assign_preliminary_x(c)
    node.x = (node.children[0].x + node.children[-1].x) / 2.0


def _next_leaf_slot(node: Node, counter: List[int]):
    """Assign leaf x positions via DFS (left to right)."""
    if not node.children:
        node.x = float(counter[0])
        counter[0] += 1
        return
    for c in node.children:
        _next_leaf_slot(c, counter)


def _resolve_overlaps(node: Node):
    if not node.children:
        return
    for c in node.children:
        _resolve_overlaps(c)

    # After children are placed, check sibling overlap from left to right
    siblings = node.children
    for i in range(1, len(siblings)):
        left = siblings[i - 1]
        right = siblings[i]
        shift = _overlap_shift(left, right)
        if shift > 0:
            _shift_subtree(right, shift)


def _overlap_shift(left: Node, right: Node) -> float:
    """How far right `right` subtree must shift so that its left contour is
    at least `separation` to the right of `left`'s right contour.
    Sibling separation is 1 (matching d3.tree default).
    """
    def right_contour(node):
        out = []
        cur = [node]
        while cur:
            out.append(max(n.x for n in cur))
            nxt = []
            for n in cur:
                nxt.extend(n.children)
            cur = nxt
        return out

    def left_contour(node):
        out = []
        cur = [node]
        while cur:
            out.append(min(n.x for n in cur))
            nxt = []
            for n in cur:
                nxt.extend(n.children)
            cur = nxt
        return out

    rc = right_contour(left)
    lc = left_contour(right)
    max_shift = 0.0
    for d in range(min(len(rc), len(lc))):
        gap = lc[d] - rc[d]
        needed = 1.0 - gap
        if needed > max_shift:
            max_shift = needed
    return max_shift


def _shift_subtree(node: Node, shift: float):
    node.x += shift
    for c in node.children:
        _shift_subtree(c, shift)
