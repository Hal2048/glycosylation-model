import ast 
from pathlib import Path  
import matplotlib  
matplotlib.use("Agg")  
import matplotlib.pyplot as plt  
from glycoform import Glycoform
from draw_glycan import draw_glycan


INPUT_FILE = Path(__file__).parent / "distribution_data.txt"
OUTPUT_FILE = Path(__file__).parent / "top10.png"
OUTPUT_SVG = Path(__file__).parent / "top10_structures.svg"
TOP_N = 10

GLYCOFORM_FIELDS = (
    "man1", "man2", "man3", "fuc", "gnb",
    "br1", "br2", "br3", "br4", "gal", "sia", "br_o",
)


def load_distribution(path):
    """Read distribution_data.txt; each line is one glycoform record."""
    entries = []  # holds all glycoforms
    with open(path, "r", encoding="utf-8") as f:  # open the file
        for line in f:  # read line by line
            line = line.rstrip("\n")  # strip the trailing newline
            if not line or line.startswith("#"):  # skip blank and comment lines
                continue
            parts = line.split("\t")  # split on tabs
            if len(parts) < 5:  # skip lines with too few fields
                continue
            try:
                tup = ast.literal_eval(parts[2])  # field 3 is a string like "(1,0,2,...)"; parse it back into a tuple
            except (ValueError, SyntaxError):
                continue
            entries.append({  # pack the needed info into a dict
                "rank": int(parts[0]),  # rank
                "id": int(parts[1]),  # glycoform ID
                "tuple": tup,  # glycoform composition (12-tuple)
                "raw": float(parts[3]),  # raw intensity
                "abundance": float(parts[4]),  # normalized abundance
            })
    return entries  # return all records to the caller


def _coerce_code(glycoform_or_tuple):
    """
    Accept either a Glycoform (namedtuple with 12 fields) or a 12-tuple and
    return the canonical (man1, man2, man3, fuc, gnb, br1, br2, br3, br4,
    gal, sia, br_o) tuple.
    """
    if hasattr(glycoform_or_tuple, "_fields") and len(glycoform_or_tuple._fields) == 12:
        return tuple(getattr(glycoform_or_tuple, name) for name in GLYCOFORM_FIELDS)
    tup = tuple(glycoform_or_tuple)
    if len(tup) != 12:
        raise ValueError(
            f"Expected 12-tuple Glycoform, got {len(tup)} values: {tup}"
        )
    return tup


def krambeck_to_iupac(glycan_form: Glycoform):
    """
    Convert the 12-digit Krambeck code to IUPAC notation
    """

    glycan_code = _coerce_code(glycan_form)
    man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o = _coerce_code(glycan_code)
    
    frag = {
        0:                                          "",
        1:                                   "(GlcNAc",
        2:                            "(Galb1-4GlcNAc",
        3:                  "(Neu5Aca2-3Galb1-4GlcNAc",
        4:                  "(GlcNAcb1-3Galb1-4GlcNAc",
        5:           "(Galb1-4GlcNAcb1-3Galb1-4GlcNAc",
        6: "(Neu5Aca2-3Galb1-4GlcNAcb1-3Galb1-4GlcNAc",
    }

    def expand(br, link):
        return f"{frag[br]}{link}" + ")" if br > 0 else ""

    # Upper arm: Br1 (b1-6) + Br2 (b1-2) -> Man(a1-6)
    if br1 > 0 or br2 > 0:
        upper = expand(br1, "b1-6") + expand(br2, "b1-2") + "Mana1-6"
    else:
        # man1/man2 are the Man counts on the a1-6 main arm / a1-3 branch
        main = "Mana1-2" * (man1 - 1) + "Mana1-6" if man1 > 0 else ""
        branch = ""
        if man2 == 2:
            branch = "(Mana1-2Mana1-3)" 
        elif man2 == 1:
            branch = "(Mana1-3)"
        upper = main + branch + "Mana1-6"

    lower = "Mana1-3"
    # Lower arm: Br3 (b1-4) + Br4 (b1-2) -> Man(a1-3)
    if br_o > 0:
        lower = "(" + "Mana1-6" * br_o + ")" + lower
    
    if br3 > 0 or br4 > 0:
        lower = expand(br3, "b1-4") + expand(br4, "b1-2")+ lower
    else:
        lower = "Mana1-2" * man3 + lower

    # Core
    core = "Manb1-4GlcNAcb1-4"
    if fuc:
        core += "(Fuca1-6)"
    core += "GlcNAcb1"

    # Bisecting GlcNAc hangs off the central Man
    if gnb:
        core = "(GlcNAcb1-4)" + core


    return upper + (f"({lower})" if lower else "") + core


def iupac_to_structure(iupac_name, output_path=None):
    """Take an IUPAC glycan name and return the SVG string; if output_path is
    given, also write it to that file."""
    svg = draw_glycan(iupac_name)
    if output_path is not None:
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(svg)
    return svg


def top10_structures(input_file=INPUT_FILE, output_file=OUTPUT_SVG, top_n=TOP_N):
    """Read the top-N glycans, convert to IUPAC and draw structures, compose
    them into one SVG with 5 columns, and label rank/ID below each subplot.
    Each sub-SVG is scaled proportionally to its intrinsic size to fit its
    cell."""
    import re
    top = sorted(load_distribution(input_file), key=lambda e: e["rank"])[:top_n]

    cols = 5
    rows = (len(top) + cols - 1) // cols
    cell_w, cell_h = 320, 495
    label_h = 25  # height reserved at the bottom of each cell for the Rank/ID text

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" '
           f'width="{cols * cell_w}" height="{rows * cell_h}">']
    for i, e in enumerate(top):
        svg = iupac_to_structure(krambeck_to_iupac(e["tuple"]))
        # the sub-SVG's actual size varies with glycan complexity; extract its intrinsic width/height
        m = re.search(r'<svg[^>]*?width="([\d.]+)"[^>]*?height="([\d.]+)"', svg)
        sub_w, sub_h = float(m.group(1)), float(m.group(2))
        # scale proportionally to fit the cell (leaving label_h at the bottom in both directions)
        s = min(cell_w / sub_w, (cell_h - label_h) / sub_h)
        drawn_w, drawn_h = sub_w * s, sub_h * s
        x, y = (i % cols) * cell_w, (i // cols) * cell_h
        x_off = x + (cell_w - drawn_w) / 2
        y_off = y + (cell_h - label_h - drawn_h) / 2
        # rewrite the sub-SVG's width/height (viewBox is kept, so content scales automatically)
        svg = re.sub(r'(<svg[^>]*?)width="[\d.]+" height="[\d.]+"',
                     rf'\1width="{drawn_w:.1f}" height="{drawn_h:.1f}"',
                     svg, count=1)
        out.append(svg.replace("<svg ", f'<svg x="{x_off:.1f}" y="{y_off:.1f}" ', 1))
        out.append(f'<text x="{x + cell_w / 2:.1f}" y="{y + cell_h - 8}" '
                   f'text-anchor="middle" font-size="14">'
                   f'Rank {e["rank"]}: ID {e["id"]}</text>')
    out.append("</svg>")

    Path(output_file).write_text("\n".join(out), encoding="utf-8")
    print(f"Image saved to: {output_file}")


def plot_top_n(entries, top_n, save_path):
    """Plot a bar chart of the top-N glycoform abundances and save it as PNG."""
    if not entries:  # exit if no data was read
        print("No glycoform data read.", file=__import__("sys").stderr)
        return

    # take the first N by rank
    top = sorted(entries, key=lambda e: e["rank"])[:top_n]
    ranks = [e["rank"] for e in top]  # bars evenly spaced by rank
    ids = [e["id"] for e in top]  # glycoform IDs shown on the X axis
    abundances = [e["abundance"] for e in top]  # Y axis: abundance

    # set up the figure
    fig, ax = plt.subplots(figsize=(max(8, top_n * 0.9), 6))
    # draw the bars: positions by rank, IDs as the X tick labels
    ax.bar(ranks, abundances, color="#3a78b5", edgecolor="#1f3a5f")

    # title and axis labels
    ax.set_xlabel("ID")
    ax.set_ylabel("Normalized Abundance")
    ax.set_title(f"Top {len(top)} Glycoforms")

    # X axis: tick positions use ranks, tick labels show glycoform IDs
    ax.set_xticks(ranks)
    ax.set_xticklabels([str(i) for i in ids], fontsize=10)

    # leave some headroom at the top of the Y axis for the value labels
    ax.set_ylim(0, max(abundances) * 1.25 if abundances else 1.0)

    # write the abundance value on top of each bar
    for x, v in zip(ranks, abundances):
        ax.text(x, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)

    # save as a PNG file, then close the figure
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Image saved to: {save_path}")


def draw_single_structure(glycoform, name):
    """
    Draw a single glycoform structure and save it as an SVG file
    """
    example_glyco1 = glycoform
    iupac_name = krambeck_to_iupac(example_glyco1)
    print(f"The IUPAC name of glycoform {name} is: {iupac_name}")
    iupac_to_structure(iupac_name, Path(f"{name}.svg"))
    print(f"The structure of glycoform {name} has been saved as {name}.svg")


if __name__ == "__main__":
    if not INPUT_FILE.exists():  # report an error and exit if the input file is missing
        print(f"Data file not found: {INPUT_FILE}", file=__import__("sys").stderr)
        raise SystemExit(1)

    entries = load_distribution(INPUT_FILE)  # step 1: read the data
    print(f"Read {len(entries)} glycoform records in total")
    plot_top_n(entries, TOP_N, OUTPUT_FILE)  # step 2: plot + save
    top10_structures()                       # step 3: combined top-10 glycan structure figure

    draw_single_structure(Glycoform(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0), "example1")


