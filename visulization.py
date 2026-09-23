GLYCOFORM_FIELDS = (
    "man1", "man2", "man3", "fuc", "gnb",
    "br1", "br2", "br3", "br4", "gal", "sia", "br_o",
)


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


def krambeck_to_iupac(code):
    """Transform the current Glycoform encoding to IUPAC-condensed notation."""
    man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o = _coerce_code(code)
    _ = (man1, man2, man3, gal, sia)  # currently used only to keep the variable alive

    # Branch completion level -> terminal fragment, mirroring Table III of
    # Krambeck & Betenbaugh (2005). Levels 5 and 6 carry an extra
    # Gal-GlcNAc repeat, hence the duplicated entry.
    frag = {
        0: "",
        1: "[GlcNAc",
        2: "[Gal(b1-4)GlcNAc",
        3: "[Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
        4: "[GlcNAc(b1-3)Gal(b1-4)GlcNAc",
        5: "[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc",
        6: "[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc",
    }

    def expand(br, link):
        return f"{frag[br]}({link})" + "]" if br > 0 else ""

    # Upper arm: Br1 (b1-6) + Br2 (b1-2) -> Man(a1-6)
    upper = (
        expand(br1, "b1-6") + expand(br2, "b1-2") + "Man(a1-6)"
        if (br1 or br2)
        else ""
    )

    # Lower arm: Br3 (b1-4) + Br4 (b1-2) -> Man(a1-3)
    lower = (
        expand(br3, "b1-4") + expand(br4, "b1-2") + "Man(a1-3)"
        if (br3 or br4)
        else ""
    )

    # Core
    core = "Man(b1-4)GlcNAc(b1-4)"
    if fuc:
        core += "[Fuc(a1-6)]"
    core += "GlcNAc(b1-N)Asn(-CHAR)"

    # Bisecting GlcNAc hangs off the central Man
    if gnb:
        core = "[GlcNAc(b1-4)]" + core

    # Yeast-like outer-chain mannose (Och1 / Mnn9). We render ``br_o``
    # alpha-1,6-linked mannose residues hanging off the central trimannosyl
    # Man of the core.
    outer = ""
    if br_o and br_o > 0:
        outer = "[" + "Man(a1-6)" * br_o + "]"

    return outer + upper + (f"[{lower}]" if lower else "") + core


if __name__ == "__main__":
    # Examples drawn from Krambeck 2005 (using the 12-tuple equivalent)
    examples = [
        # (man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o)
        (2, 0, 1, 1, 1, 6, 5, 3, 1, 5, 2, 0),  # paper Figure 1c
        (1, 1, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0),  # biantennary
        (1, 1, 0, 1, 0, 0, 2, 0, 2, 2, 0, 0),  # Table XIII 310020220
        (0, 1, 0, 1, 0, 2, 2, 2, 2, 4, 0, 0),  # tetra-antennary
        (2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0),  # M9 (initial substrate)
        (1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1),  # M8 with one outer Man
    ]

    for code in examples:
        print(f"{code} -> {krambeck_to_iupac(code)}")



