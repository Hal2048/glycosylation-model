"""Render a glycan tree to a self-contained SVG string.

Ported from GlycoGlyph `src/modules/d3glycanstruc.js` and `savesvg.js`
(inline-monos behaviour). The output SVG embeds every used monosaccharide
symbol inline, so it can be saved to a file or rendered without external
resources.
"""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .tree_layout import Node, hierarchy, tree_layout


# Default drawing settings (mirror GlycoGlyph's `drawingSettings`).
DEFAULTS: Dict[str, Any] = {
    'width': 320,
    'height': 465,
    'orientation': 'bottom-to-top',   # 'bottom-to-top' or 'right-to-left'
    'symbsize': 25,
    'drawingareachoice': 'fixglycansize',
    'fucopt': 'fucdown',              # 'fucdown' | 'fucout' | 'fucleft' | 'fucoriginal'
    'linkageVisible': True,
    'linkRotate': True,
    'linkAbbr': False,
    'linkFontSize': 12,
    'margin': {'top': 40, 'right': 40, 'bottom': 40, 'left': 40},
    'inline_symbols': True,            # embed monos.svg symbols inline
}


# ----- monos.svg loading / caching -----

_SYMBOLS_CACHE: Optional[Dict[str, Tuple[str, str]]] = None  # id -> (viewBox, innerSVG)


def _load_symbols() -> Dict[str, Tuple[str, str]]:
    """Parse monos.svg once and return {id: (viewBox, inner_svg)}."""
    global _SYMBOLS_CACHE
    if _SYMBOLS_CACHE is not None:
        return _SYMBOLS_CACHE
    path = Path(__file__).parent / 'monos.svg'
    text = path.read_text(encoding='utf-8')
    # strip xml namespace to simplify ElementTree parsing
    text_clean = re.sub(r'xmlns(:\w+)?="[^"]+"', '', text)
    root = ET.fromstring(text_clean)
    symbols: Dict[str, Tuple[str, str]] = {}
    for sym in root.iter('symbol'):
        sid = sym.get('id')
        viewbox = sym.get('viewBox') or '0 0 50 50'
        # serialize inner children
        inner = ''.join(ET.tostring(child, encoding='unicode') for child in list(sym))
        # remove namespace prefixes ET may add
        inner = re.sub(r'ns\d+:', '', inner)
        symbols[sid] = (viewbox, inner)
    _SYMBOLS_CACHE = symbols
    return symbols


# ----- main drawing entry -----

def render_svg(glycan_name: str, options: Optional[Dict[str, Any]] = None) -> str:
    """Render a CFG glycan name to a self-contained SVG string."""
    # parse name to tree data, then sort children BEFORE layout
    # (mirrors GlycoGlyph's order: parse -> hierarchy -> sortchildren -> treemap)
    from .parser import glycan_to_json  # local import to avoid circular

    cfg = dict(DEFAULTS)
    cfg['margin'] = dict(DEFAULTS['margin'])
    if options:
        for k, v in options.items():
            if k == 'margin' and isinstance(v, dict):
                cfg['margin'].update(v)
            else:
                cfg[k] = v

    data = glycan_to_json(glycan_name)
    root = hierarchy(data)
    depth = _get_depth(data)

    margin = cfg['margin']
    symb_size = cfg['symbsize']

    # In fixglycansize mode, give the layout enough horizontal space
    # based on the number of leaves, so multi-branch glycans are not
    # squeezed into the default width
    n_leaves = max(1, len(root.leaves()))
    leaf_spacing = int(symb_size * 2.4)

    if cfg['orientation'] == 'right-to-left':
        # swap width/height for layout purposes (as original does)
        width = cfg['height'] - margin['top'] - margin['bottom']
        height = cfg['width'] - margin['left'] - margin['right']
        # then "fixdrawingarea" branch sets the other dim
        if cfg['drawingareachoice'] == 'fixdrawingarea':
            height = cfg['height'] - margin['top'] - margin['bottom']
            width = cfg['width'] - margin['left'] - margin['right']
        elif cfg['drawingareachoice'] == 'fixglycansize':
            if depth <= 2:
                width = 50 * depth
            elif depth < 9:
                width = 50 * depth
            else:
                width = cfg['width'] - margin['left'] - margin['right']
            height = max(height, n_leaves * leaf_spacing)
        # orientations mapping (size is [height, width] for r2l)
        size = (height, width)
    else:
        # bottom-to-top (default)
        width = cfg['width'] - margin['left'] - margin['right']
        height = cfg['height'] - margin['top'] - margin['bottom']
        if cfg['drawingareachoice'] == 'fixdrawingarea':
            height = cfg['height'] - margin['top'] - margin['bottom']
        elif cfg['drawingareachoice'] == 'fixglycansize':
            if depth <= 2:
                height = 50 * depth
            elif depth < 9:
                height = 50 * depth
            else:
                height = cfg['height'] - margin['top'] - margin['bottom']
            width = max(width, n_leaves * leaf_spacing)
        size = (width, height)

    # sort children before layout (as GlycoGlyph does)
    if cfg['orientation'] == 'right-to-left':
        _sort_children(root, reverse=True)
    else:
        _sort_children(root, reverse=False)

    # run tidy tree layout
    tree_layout(root, size)

    # fix fucose positions (after layout, before drawing)
    _fix_fucers(root, cfg)

    # Coordinate transform functions (defined early so we can compute the
    # content bounding box)
    if cfg['orientation'] == 'right-to-left':
        def ox(n): return size[1] - n.y
        def oy(n): return n.x
    else:
        def ox(n): return n.x
        def oy(n): return size[1] - n.y

    # Compute the bounding box of the actual drawn content (in layout
    # coordinates), including the symbol radius and the margin for linkage
    # text labels
    half = symb_size / 2
    link_pad = (cfg['linkFontSize'] + 4) if cfg['linkageVisible'] else 0
    desc = list(root.descendants())
    cx0 = min(ox(n) for n in desc) - half - link_pad
    cx1 = max(ox(n) for n in desc) + half + link_pad
    cy0 = min(oy(n) for n in desc) - half - link_pad
    cy1 = max(oy(n) for n in desc) + half + link_pad

    # ----- build SVG -----
    # group translation (horizontal centering)
    if cfg['orientation'] == 'right-to-left':
        dxg = root.x - size[0] / 2
        if dxg != 0:
            if depth < 10:
                tx = size[1] + margin['left'] - 50 * depth
                ty = margin['top'] - dxg
            else:
                tx = margin['left']
                ty = margin['top'] - dxg
        else:
            tx = margin['left']
            ty = margin['top']
    else:
        dxg = root.x - size[0] / 2
        if dxg != 0:
            tx = margin['left'] - dxg
            ty = margin['top']
        else:
            tx = margin['left']
            ty = margin['top']

    # Translate the content bounding box into canvas coordinates; if it goes
    # beyond the left/top edge, adjust tx/ty to push it back inside
    abs_x0 = tx + cx0
    abs_x1 = tx + cx1
    abs_y0 = ty + cy0
    abs_y1 = ty + cy1
    pad = 5
    if abs_x0 < pad:
        shift = pad - abs_x0
        tx += shift
        abs_x1 += shift
        abs_x0 = pad
    if abs_y0 < pad:
        shift = pad - abs_y0
        ty += shift
        abs_y1 += shift
        abs_y0 = pad

    # Canvas size = content bounding box + padding (snug fit).
    # If the user explicitly set width/height in options, take the larger
    # value to respect their setting
    svg_w = int(math.ceil(abs_x1 + pad))
    svg_h = int(math.ceil(abs_y1 + pad))
    if options:
        if 'width' in options:
            svg_w = max(svg_w, int(options['width']))
        if 'height' in options:
            svg_h = max(svg_h, int(options['height']))

    parts: List[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{_fmt(svg_w)}" height="{_fmt(svg_h)}" '
        f'viewBox="0 0 {_fmt(svg_w)} {_fmt(svg_h)}">'
    )

    # background
    parts.append(
        f'<rect width="{_fmt(svg_w)}" height="{_fmt(svg_h)}" fill="white"/>'
    )

    parts.append(f'<g transform="translate({_fmt(tx)},{_fmt(ty)})">')

    # ----- draw links -----
    for n in root.descendants():
        if n is root:
            continue
        p = n.parent
        x1, y1 = ox(n), oy(n)
        x2, y2 = ox(p), oy(p)
        parts.append(
            f'<path d="M{_fmt(x1)},{_fmt(y1)}L{_fmt(x2)},{_fmt(y2)}" '
            f'stroke-width="2" stroke="black" fill="none"/>'
        )

    # ----- draw nodes -----
    symbols = _load_symbols() if cfg['inline_symbols'] else {}
    symb_size = cfg['symbsize']
    for idx, n in enumerate(root.descendants()):
        transform = _rotate_fuc(n, ox, oy, cfg)
        parts.append(f'<g transform="{transform}">')

        # mono symbol
        mono = n.data.get('monosaccharide') or n.data.get('name') or ''
        sid = mono.lower() if mono else 'unknown'
        if cfg['inline_symbols'] and sid in symbols:
            viewbox, inner = symbols[sid]
            half = symb_size / 2
            parts.append(
                f'<svg x="{_fmt(-half)}" y="{_fmt(-half)}" '
                f'width="{_fmt(symb_size)}" height="{_fmt(symb_size)}" '
                f'viewBox="{viewbox}">{inner}</svg>'
            )
        else:
            # unknown / missing — draw a placeholder
            half = symb_size / 2
            parts.append(
                f'<rect x="{_fmt(-half)}" y="{_fmt(-half)}" '
                f'width="{_fmt(symb_size)}" height="{_fmt(symb_size)}" '
                f'fill="white" stroke="black" stroke-width="1"/>'
            )
            parts.append(
                f'<text x="0" y="0" font-size="8" text-anchor="middle" '
                f'dominant-baseline="central" font-family="Arial">?</text>'
            )

        # linkage text
        if cfg['linkageVisible']:
            link_txt = _linkage_text(n, idx, cfg)
            if link_txt:
                transform_link = _transform_link_text(n, idx, cfg, ox, oy)
                parts.append(
                    f'<text transform="{transform_link}" '
                    f'font-size="{cfg["linkFontSize"]}" text-anchor="middle" '
                    f'font-family="Arial,Sans">{_xml_escape(link_txt)}</text>'
                )

        # substituent text
        sub_txt = _sub_text(n)
        if sub_txt:
            if cfg['orientation'] == 'bottom-to-top':
                sx = -cfg['linkFontSize'] - 3
                sy = 3
                anchor = 'end'
            else:
                sx = 3
                sy = cfg['linkFontSize'] + symb_size / 2
                anchor = 'middle'
            parts.append(
                f'<text x="{_fmt(sx)}" y="{_fmt(sy)}" font-size="{cfg["linkFontSize"]}" '
                f'text-anchor="{anchor}" font-family="Arial,Sans">'
                f'{_xml_escape(sub_txt)}</text>'
            )

        parts.append('</g>')

    parts.append('</g>')
    parts.append('</svg>')
    return ''.join(parts)


# ----- helpers -----

def _fmt(v):
    if isinstance(v, float):
        s = f'{v:.4f}'.rstrip('0').rstrip('.')
        return s if s else '0'
    return str(v)


def _xml_escape(s: str) -> str:
    return (s.replace('&', '&amp;')
             .replace('<', '&lt;')
             .replace('>', '&gt;'))


def _get_depth(obj: Dict[str, Any]) -> int:
    depth = 0
    for c in obj.get('children', []) or []:
        d = _get_depth(c)
        if d > depth:
            depth = d
    return 1 + depth


def _sort_children(node: Node, reverse: bool = False):
    """Mirror sortchildren.js — sort by linknum; special core-fucose case."""
    if node.children and len(node.children) > 0:
        if reverse:
            node.children.sort(key=lambda c: c.data.get('linknum', 0), reverse=True)
        else:
            node.children.sort(key=lambda c: c.data.get('linknum', 0))

        # core fucose detection (only when right-to-left / reverse=True in
        # the original code — original checks regardless of reverse flag but
        # the sort it does is the custom one; we follow the same logic)
        has_fuc = any(c.data.get('monosaccharide') == 'Fuc' for c in node.children)
        has_glcnac_sib = any(c.data.get('monosaccharide') == 'GlcNAc' for c in node.children)
        if (has_fuc
                and node.depth <= 1
                and node.data.get('monosaccharide') == 'GlcNAc'
                and has_glcnac_sib):
            node.children.sort(
                key=lambda c: 0 if c.data.get('monosaccharide') == 'Fuc' else 1
            )

        for c in node.children:
            _sort_children(c, reverse)


def _fix_fucers(node: Node, cfg: Dict[str, Any]):
    """Mirror fixfucers() in d3glycanstruc.js."""
    if cfg['fucopt'] == 'fucoriginal':
        return

    # width used for fucose shift
    if cfg['orientation'] == 'right-to-left':
        width = cfg['height'] - cfg['margin']['top'] - cfg['margin']['bottom']
    else:
        width = cfg['width'] - cfg['margin']['left'] - cfg['margin']['right']

    centerx = node.x

    for n in node.descendants():
        if n.data.get('monosaccharide') != 'Fuc':
            continue
        if n.depth == 0:
            continue
        if not n.parent or not n.parent.children:
            continue

        siblings = n.parent.children
        xmin = siblings[0].x
        xmax = siblings[-1].x
        nonfuc = sum(1 for s in siblings if 'Fuc' not in (s.data.get('name') or ''))

        # drop fucose to parent y
        n.y = n.parent.y
        parentx = n.parent.x

        if cfg['fucopt'] == 'fucout':
            if n.x < centerx:
                n.x = n.parent.x - width / 5
            elif n.x > centerx:
                n.x = n.parent.x + width / 5
            else:
                if n.parent.x < centerx:
                    n.x = n.parent.x - width / 5
                elif n.parent.x > centerx:
                    n.x = n.parent.x + width / 5
                else:
                    n.x = n.parent.x - width / 5
        elif cfg['fucopt'] == 'fucdown':
            if n.x < n.parent.x:
                n.x = n.parent.x - width / 5
            else:
                n.x = n.parent.x + width / 5
        elif cfg['fucopt'] == 'fucleft':
            n.x = n.parent.x - width / 5

        # find first non-fucose sibling diff
        diffx = None
        for s in siblings:
            if s.data.get('monosaccharide') != 'Fuc' and diffx is None:
                diffx = s.x - parentx
        if diffx is not None:
            for s in siblings:
                if s.data.get('monosaccharide') != 'Fuc' and s.data.get('name') != n.parent.data.get('name'):
                    for child in s.descendants():
                        child.x = child.x - diffx


def _rotate_fuc(n: Node, ox, oy, cfg: Dict[str, Any]) -> str:
    if n.depth == 0:
        return f'translate({_fmt(ox(n))},{_fmt(oy(n))})'
    if cfg['orientation'] == 'right-to-left':
        return f'translate({_fmt(ox(n))},{_fmt(oy(n))})'
    # default behavior: non-original => rotate fucose
    if cfg['fucopt'] == 'fucoriginal':
        return f'translate({_fmt(ox(n))},{_fmt(oy(n))})'
    if n.data.get('monosaccharide') != 'Fuc':
        return f'translate({_fmt(ox(n))},{_fmt(oy(n))})'
    # rotate ±30° around its position; bump y slightly
    n.y = n.y + 4
    if n.x < n.parent.x:
        return f'translate({_fmt(ox(n))},{_fmt(oy(n))}) rotate(-30)'
    else:
        return f'translate({_fmt(ox(n))},{_fmt(oy(n))}) rotate(30)'


def _linkage_text(n: Node, idx: int, cfg: Dict[str, Any]) -> str:
    s = n.data.get('linkage') or ''
    if s == '':
        return ''
    if idx == 0:
        return s.replace('a', 'α').replace('b', 'β')

    abbr = cfg['linkAbbr']
    if cfg['orientation'] == 'bottom-to-top' \
            and n.data.get('monosaccharide') == 'Fuc' \
            and n.x > n.parent.x:
        if abbr:
            s = re.sub(r'[\d?]\-', '', s)
        s = s[::-1]
        s = re.sub(r'a', 'α', s)
        s = re.sub(r'b', 'β', s)
        return s

    if abbr:
        s = re.sub(r'a(?=[\d?]?\-[\d?])', 'α', s)
        s = re.sub(r'b(?=[\d?]?\-[\d?])', 'β', s, count=1)
        s = re.sub(r'[\d?]\-', '', s)
        return s

    if re.search(r'[a,b,?][\d,?]\-[\d,?]', s):
        if re.search(r'a[\d,?]\-[\d,?]', s):
            s = re.sub(r'a(?=[\d?]\-[\d?])', 'α', s, count=1)
        elif re.search(r'b[\d?]\-[\d?]', s):
            s = re.sub(r'b(?=[\d?]\-[\d?])', 'β', s, count=1)
        return s[max(s.index('-') - 2, 0):]
    elif re.search(r'[a,b,?]\-', s):
        return s[max(s.index('-') - 1, 0):]
    elif re.search(r'[a,b,?][\d,?]', s):
        if len(s) >= 2 and s[-2] == 'a':
            return 'α' + s[-1]
        elif len(s) >= 2 and s[-2] == 'b':
            return 'β' + s[-1]
        else:
            return s
    elif s == 'a':
        return 'α'
    elif s == 'b':
        return 'β'
    else:
        return s


def _transform_link_text(n: Node, idx: int, cfg: Dict[str, Any], ox, oy) -> str:
    if idx == 0:
        return f'translate({_fmt(cfg["symbsize"] / 4)},{_fmt(cfg["symbsize"])})'

    if cfg['orientation'] == 'right-to-left':
        angle = math.degrees(math.atan2(n.y - n.parent.y, n.parent.x - n.x))
        if n.data.get('monosaccharide') == 'Fuc' and cfg.get('fuctype', '') != 'fucoriginal':
            if angle <= 45:
                angle = 90
                return (f'translate({_fmt((n.y - n.parent.y) / 2 + cfg["linkFontSize"] / 2)},'
                        f'{_fmt((n.parent.x - n.x) / 2)}) rotate({_fmt(angle)})')
            else:
                angle = -90
                return (f'translate({_fmt((n.y - n.parent.y) / 2 - cfg["linkFontSize"] / 2)},'
                        f'{_fmt((n.parent.x - n.x) / 2)}) rotate({_fmt(angle)})')
        if not cfg['linkRotate']:
            if n.x < n.parent.x:
                return (f'translate({_fmt((n.y - n.parent.y) / 2 + cfg["linkFontSize"])},'
                        f'{_fmt((n.parent.x - n.x) / 2 - cfg["linkFontSize"])})')
            elif n.x == n.parent.x:
                return (f'translate({_fmt((n.y - n.parent.y) / 2)},'
                        f'{_fmt((n.parent.x - n.x) / 2 - cfg["linkFontSize"] / 2)})')
            else:
                return (f'translate({_fmt((n.y - n.parent.y) / 2 + cfg["linkFontSize"])},'
                        f'{_fmt((n.parent.x - n.x) / 2 + cfg["linkFontSize"])})')
        angle = math.degrees(math.atan2(n.parent.x - n.x, n.y - n.parent.y))
        return (f'translate({_fmt((n.y - n.parent.y) / 2)},'
                f'{_fmt((n.parent.x - n.x) / 2 - cfg["linkFontSize"] / 2)}) rotate({_fmt(angle)})')

    # bottom-to-top
    angle = math.degrees(math.atan2(n.y - n.parent.y, n.parent.x - n.x))
    if n.data.get('monosaccharide') == 'Fuc' and cfg.get('fuctype', '') != 'fucoriginal' and n.depth > 0:
        angle = 30 if n.x < n.parent.x else 150
    if n.data.get('monosaccharide') == 'Fuc' and cfg.get('fuctype', '') != 'fucoriginal':
        if angle <= 45:
            return (f'translate({_fmt((n.parent.x - n.x - 10) / 2)},'
                    f'{_fmt((n.y - n.parent.y + cfg["symbsize"] - 6) / 2)}) rotate({_fmt(angle)})')
        else:
            angle = -30
            return (f'translate({_fmt((n.parent.x - n.x + 10) / 2)},'
                    f'{_fmt((n.y - n.parent.y + cfg["symbsize"] - 6) / 2)}) rotate({_fmt(angle)})')
    if not cfg['linkRotate']:
        n_val = 2 if cfg['linkAbbr'] else 4
        if n.x < n.parent.x:
            return (f'translate({_fmt((n.parent.x - n.x - 6 * cfg["linkFontSize"] + cfg["linkFontSize"] * n_val) / (cfg["linkFontSize"] / 2))},'
                    f'{_fmt((n.y - n.parent.y + 5) / 2)})')
        elif n.x == n.parent.x:
            return (f'translate({_fmt((n.parent.x - n.x + 6 * cfg["linkFontSize"] + cfg["linkFontSize"] * n_val) / (cfg["linkFontSize"] / 2))},'
                    f'{_fmt((n.y - n.parent.y + 5) / 2)})')
        else:
            return (f'translate({_fmt((n.parent.x - n.x + 6 * cfg["linkFontSize"] - cfg["linkFontSize"] * n_val) / (cfg["linkFontSize"] / 2))},'
                    f'{_fmt((n.y - n.parent.y + 5) / 2)})')
    return (f'translate({_fmt((n.parent.x + 10 - n.x) / 2)},'
            f'{_fmt((n.y - n.parent.y) / 2)}) rotate({_fmt(angle)})')


def _sub_text(n: Node) -> str:
    s = n.data.get('substituent')
    if not s:
        return ''
    s = re.sub(r'\[(?=[0-9a-zA-Z])', '', s)
    s = re.sub(r'\](?=[0-9a-zA-Z])', ',', s)
    s = s.replace(']', '')
    return s
