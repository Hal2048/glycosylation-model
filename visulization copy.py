from glycoform import Glycoform

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

def krambeck_to_iupac(glycan_form: Glycoform):
    """
    将12位Krambeck编码转换为IUPAC缩略表示法
    
    编码格式：[man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o]
    
    基础结构: Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-N)Asn
    man1: a1-6分支上的a1-2延伸Man数
    man2: a1-3分支上的a1-2延伸Man数  
    man3: a1-3分支上额外a1-2延伸Man数
    br_o: a1-3分支上a1-6延伸的Man数
    """

    glycan_code = _coerce_code(glycan_form)
    man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o = _coerce_code(glycan_code)
    
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
    if br1 > 0 or br2 > 0:
        upper = expand(br1, "b1-6") + expand(br2, "b1-2") + "Man(a1-6)"
    elif man1 == 2 and man2 == 2:
        upper = "Man(a1-2)Man(a1-6)[Man(a1-2)Man(a1-3)]" + "Man(a1-6)"
    elif man1 == 2 and man2 == 1:
        upper = "Man(a1-2)Man(a1-6)[Man(a1-3)]" + "Man(a1-6)"
    elif man1 == 1 and man2 == 2:
        upper = "Man(a1-6)[Man(a1-2)Man(a1-3)]" + "Man(a1-6)"
    elif man1 == 1 and man2 == 1:
        upper = "Man(a1-6)[Man(a1-3)]" + "Man(a1-6)"

    # Lower arm: Br3 (b1-4) + Br4 (b1-2) -> Man(a1-3)
    lower = "Man(a1-3)"
    if br_o > 0:
        lower = "[" + "Man(a1-6)" * br_o + "]" + "Man(a1-3)"
    
    if br3 > 0 or br4 > 0:
        lower = expand(br3, "b1-4") + expand(br4, "b1-2")+ lower
    else:
        lower = "Man(a1-2)" * man3 + lower

    # Core
    core = "Man(b1-4)GlcNAc(b1-4)"
    if fuc:
        core += "[Fuc(a1-6)]"
    core += "GlcNAc(b1-N)Asn(-CHAR)"

    # Bisecting GlcNAc hangs off the central Man
    if gnb:
        core = "[GlcNAc(b1-4)]" + core


    return upper + (f"[{lower}]" if lower else "") + core


if __name__ == "__main__":
    examples = [
        Glycoform(man1=0, man2=0, man3=0, fuc=0, gnb=0, br1=2, br2=5, br3=1, br4=5, gal=5, sia=0, br_o=1),
    ]


    for code in examples:
        print(f"{code} -> {krambeck_to_iupac(code)}")