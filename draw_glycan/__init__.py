"""draw_glycan — Python interface to render CFG glycan names as SVG.

Usage:
    from draw_glycan import draw_glycan

    svg = draw_glycan("Galb1-4GlcNAcb1-2Mana1-3(Mana1-6)Manb1-4GlcNAcb1-4GlcNAcb1")
    with open("glycan.svg", "w", encoding="utf-8") as f:
        f.write(svg)
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from .render import render_svg


def draw_glycan(iupac_name: str, options: Optional[Dict[str, Any]] = None) -> str:
    """Render a CFG/IUPAC glycan name to a self-contained SVG string.

    Args:
        iupac_name: Glycan sequence using CFG symbols.
            Examples:
                "Galb1-4GlcNAcb1-2Mana1-3(Mana1-6)Manb1-4GlcNAcb1-4GlcNAcb1"
                "Mana1-6(Mana1-3)Manb1-4GlcNAcb1-4GlcNAcb1"
                "Galb1-3GalNAca1"
        options: Optional overrides. Supported keys:
            width (int): canvas width in px (default 320)
            height (int): canvas height in px (default 465)
            orientation (str): 'bottom-to-top' (default) or 'right-to-left'
            symbsize (int): monosaccharide symbol size in px (default 25)
            linkageVisible (bool): show linkage text (default True)
            linkFontSize (int): font size for linkage text (default 12)
            margin (dict): {'top':40,'right':40,'bottom':40,'left':40}
            inline_symbols (bool): embed symbol geometry inline (default True)

    Returns:
        SVG string. Save with `open(..., 'w', encoding='utf-8')`.
    """
    return render_svg(iupac_name, options)


__all__ = ['draw_glycan']
