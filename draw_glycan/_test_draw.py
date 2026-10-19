"""Quick smoke test for draw_glycan — generates SVGs for several structures."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from draw_glycan import draw_glycan

TEST_STRUCTURES = [
    "Galb1-4GlcNAcb1",                              # disaccharide
    "Mana1-6(Mana1-3)Manb1-4GlcNAcb1-4GlcNAcb1",    # N-core
    "Galb1-3GalNAca1",                              # Core-1 O-glycan
    "Fuca1-2Galb1-3(Fuca1-4)GlcNAcb1",              # Lewis B
    "Neu5Aca2-3Galb1-4GlcNAcb1",                    # sialyl LacNAc
    "Mana1-6(GlcNAcb1-4)(Mana1-3)Manb1-4GlcNAcb1-4(Fuca1-6)GlcNAcb1",  # bisecting + core fuc
]

out_dir = Path(__file__).parent / "_test_output"
out_dir.mkdir(exist_ok=True)

for i, name in enumerate(TEST_STRUCTURES):
    print(f"[{i+1}/{len(TEST_STRUCTURES)}] {name}")
    try:
        svg = draw_glycan(name)
        out_file = out_dir / f"test_{i+1}.svg"
        out_file.write_text(svg, encoding='utf-8')
        print(f"    -> {out_file.name} ({len(svg)} bytes)")
    except Exception as e:
        print(f"    ERROR: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
