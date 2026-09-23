"""
glycan_drawer — minimal Python interface for drawing glycan structures.

Quick start:
    from glycan_drawer import draw
    draw("Gal(b1-4)Glc", "out.svg")
    draw("Mana1-6(Mana1-3)Manb1-4GlcNAcb1-4GlcNAcb1", "out.png")

Setup (once):
    cd glycan_drawer/node && npm install

How it works:
    The drawing core (d3 + GlycoGlyph) runs headlessly under Node.js + jsdom — no real
    browser is needed. Emitted SVGs are self-contained (every <use> is inlined). PNG
    conversion falls back from cairosvg to a headless Edge/Chrome screenshot.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_NODE_DIR = _HERE / "node"
_DRAW_JS = _NODE_DIR / "draw.js"

# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

def _ensure_node_deps() -> None:
    """Auto-run `npm install` if node_modules is missing (first-time convenience)."""
    if (_NODE_DIR / "node_modules").exists():
        return
    if not shutil.which("node") or not shutil.which("npm"):
        raise RuntimeError("Node.js + npm are required (https://nodejs.org).")
    subprocess.run(["npm", "install"], cwd=str(_NODE_DIR), check=True)


def _render_svg(iupac: str, width: int, height: int, orientation: str) -> str:
    """Call `node draw.js --json` and parse the resulting SVG string."""
    _ensure_node_deps()
    proc = subprocess.run(
        ["node", str(_DRAW_JS), iupac, "--json", str(width), str(height), orientation],
        cwd=str(_NODE_DIR), capture_output=True, text=True, timeout=60,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"node draw.js failed (code {proc.returncode}):\n{proc.stderr}")
    try:
        result = json.loads(proc.stdout)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Bad Node output: {e}\n{proc.stdout[:300]}") from e
    if result.get("error"):
        raise RuntimeError(f"Rendering failed: {result['error']}")
    return result["svg"]


def _find_browser() -> str | None:
    """Locate Edge / Chrome in common install paths, then PATH."""
    for p in (
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    ):
        if Path(p).exists():
            return p
    for name in ("msedge", "chrome", "chromium"):
        if found := shutil.which(name):
            return found
    return None


def _svg_to_png(svg: str, out_path: Path, scale: float) -> None:
    """Try cairosvg first, then fall back to a headless-browser screenshot."""
    try:
        import cairosvg
        cairosvg.svg2png(
            bytestring=svg.encode("utf-8"),
            write_to=str(out_path),
            output_width=int(scale) if scale else None,
            background_color="white",
        )
        return
    except Exception:
        pass

    browser = _find_browser()
    if not browser:
        raise RuntimeError("PNG needs cairosvg or a headless Edge/Chrome installation.")

    # Use SVG's declared size as the base for the screenshot window.
    w = float((re.search(r'\bwidth="([\d.]+)"', svg) or [None, "320"])[1])
    h = float((re.search(r'\bheight="([\d.]+)"', svg) or [None, "500"])[1])
    s = scale if scale > 0 else 2.0

    with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False, encoding="utf-8") as f:
        f.write(svg)
        tmp_svg = f.name
    try:
        cmd = [
            browser, "--headless=new", "--disable-gpu", "--no-sandbox",
            f"--window-size={int(w * s)},{int(h * s)}",
            f"--screenshot={out_path.resolve()}",
            "file:///" + tmp_svg.replace("\\", "/"),
        ]
        if subprocess.run(cmd, capture_output=True, timeout=60).returncode != 0 or not out_path.exists():
            raise RuntimeError("Headless browser screenshot failed.")
    finally:
        try: os.unlink(tmp_svg)
        except OSError: pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def to_svg(iupac: str, width: int = 320, height: int = 500, orientation: str = "bottom-to-top") -> str:
    """Render and return the SVG as a self-contained string."""
    return _render_svg(iupac, width, height, orientation)


def to_png(iupac: str, width: int = 320, height: int = 500, orientation: str = "bottom-to-top", scale: float = 2.0) -> bytes:
    """Render and return PNG bytes. Tries cairosvg, then headless Edge/Chrome."""
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tmp_png = Path(tmp.name)
    try:
        _svg_to_png(_render_svg(iupac, width, height, orientation), tmp_png, scale)
        return tmp_png.read_bytes()
    finally:
        try: tmp_png.unlink()
        except OSError: pass


def draw(
    iupac: str,
    output: str | os.PathLike,
    width: int = 320,
    height: int = 500,
    orientation: str = "bottom-to-top",
    scale: float = 2.0,
) -> Path:
    """Render to a file. Format is auto-detected from the extension (.png -> PNG, else SVG)."""
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    svg = _render_svg(iupac, width, height, orientation)
    if out_path.suffix.lower() == ".png":
        _svg_to_png(svg, out_path, scale)
    else:
        out_path.write_text(svg, encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description="IUPAC-condensed -> glycan structure image (SVG/PNG).")
    ap.add_argument("iupac", help="IUPAC-condensed glycan, e.g. 'Gal(b1-4)GlcNAc'.")
    ap.add_argument("-o", "--output", default="out.svg")
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=500)
    ap.add_argument("--orientation", choices=["bottom-to-top", "right-to-left"], default="bottom-to-top")
    ap.add_argument("--scale", type=float, default=2.0)
    args = ap.parse_args()
    p = draw(args.iupac, args.output, args.width, args.height, args.orientation, args.scale)
    print(f"Wrote {p} ({p.stat().st_size} bytes)")


if __name__ == "__main__":
    main()