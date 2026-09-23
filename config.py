import numpy as np


# Residence time for each compartment (CGC, MGC, TGC, TGN)
TAU = np.array([5.56, 5.56, 5.56, 5.56])

# Donor concentrations (uM)
DONOR_CONC = {
    "UDP-GlcNAc": 9200.0,
    "UDP-Gal": 3800.0,
    "CMP-NeuAc": 2400.0,
    "GDP-Fuc": 5000.0,
    "GDP-Man": 2000.0,
    "H20": 1
}

# Total enzyme concentrations (uM), Compartment distribution ratio (CGC, MGC, TGC, TGN)
ENZYME_CONC = {
    "ManI_9":   [1.78, [0.05, 0.15, 0.40, 0.40]],
    "ManI_8":   [1.78, [0.05, 0.15, 0.40, 0.40]],
    "ManI_7":   [1.78, [0.05, 0.15, 0.40, 0.40]],
    "ManI_6":   [1.78, [0.05, 0.15, 0.40, 0.40]],
    "ManII_5":  [1.32, [0.20, 0.40, 0.25, 0.15]],
    "ManII_4":  [1.32, [0.20, 0.40, 0.25, 0.15]],
    "FucT":     [2.50, [0.05, 0.15, 0.40, 0.40]],
    "GnTI":     [3.05, [0.15, 0.30, 0.35, 0.20]],
    "GnTII":    [1.29, [0.10, 0.25, 0.40, 0.25]],
    "GnTIV":    [3.62, [0.10, 0.20, 0.35, 0.35]],
    "GnTV":     [0.4,  [0.10, 0.25, 0.30, 0.35]],
    "GnTE_Br1": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "GnTE_Br2": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "GnTE_Br3": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "GnTE_Br4": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "GalT_Br1": [0.66, [0.05, 0.20, 0.35, 0.40]],
    "GalT_Br2": [0.66, [0.05, 0.20, 0.35, 0.40]],
    "GalT_Br3": [0.66, [0.05, 0.20, 0.35, 0.40]],
    "GalT_Br4": [0.66, [0.05, 0.20, 0.35, 0.40]],
    "SiaT_Br1": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "SiaT_Br2": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "SiaT_Br3": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "SiaT_Br4": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "Och1":     [0,    [0.25, 0.25, 0.25, 0.25]],
    "Mnn9":     [0,    [0.25, 0.25, 0.25, 0.25]],
    "GnTIII":   [0,    [0.25, 0.25, 0.25, 0.25]],
}




ENZYME_RULES = [
    # ManI, 1
    {
        "name": "ManI_M9",
        "condition": "g.man1 == 2 and g.man2 == 2 and g.man3 == 2", # M9
        "product_struc": "Glycoform(g.man1, g.man2, g.man3-1, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManI, 2
    {
        "name": "ManI_M8",
        "condition": "g.man1 == 2 and g.man2 == 2 and g.man3 == 1", # M8
        "product_struc": "Glycoform(g.man1-1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManI, 3
    {
        "name": "ManI_M7",
        "condition": "g.man1 == 1 and g.man2 == 2 and g.man3 == 1", # M7
        "product_struc": "Glycoform(g.man1, g.man2, g.man3-1, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManI, 3
    {
        "name": "ManI_M6",
        "condition": "g.man1 == 1 and g.man2 == 2 and g.man3 == 0", # M6
        "product_struc": "Glycoform(g.man1, g.man2-1, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManII, 1
    {
        "name": "ManII_M5",
        "condition": "g.man1 ==1 and g.man2 == 1 and g.man3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1, g.man2-1, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",    
        "Km": 200, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManII, 2
    {
        "name": "ManII_M4",
        "condition": "g.man1 ==1 and g.man2 == 0 and g.man3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1-1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",   
        "Km": 200, "Kmd": 0, "kf": 1924,
        "adjustments": [("g.man1 ==1 and g.man2 == 0", 0.50)]
    },

    # FucT
    {
        "name": "FucT",
        "condition": "g.fuc == 0 and g.br4 > 0 and g.gnb == 0 and g.gal == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc+1, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "GDP-Fuc",
        "Km": 25, "Kmd": 46, "kf": 253,
        "adjustments": []
    },

    # GnTI
    {
        "name": "GnTI",
        "condition": "g.br4 == 0 and g.man1 + g.man2 + g.man3 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 260, "Kmd": 170, "kf": 990,
        "adjustments": []
    },

    # GnTII
    {
        "name": "GnTII",
        "condition": "g.br2 == 0 and g.man1 + g.man2 + g.man3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 190, "Kmd": 3100, "kf": 607,
        "adjustments": []
    },

    # # GnTIII
    # {
    #     "name": "GnTIII",
    #     "condition": "g.gnb == 0 and g.br4 > 0 and g.gal == 0",
    #     "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb+1, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
    #     "cosubstrate": "UDP-GlcNAc",
    #     "Km": 4000, "Kmd": 3100, "kf": 607,
    #     "adjustments": [("g.br2 > 0", 0.048)]
    # },

    # GnTIV
    {
        "name": "GnTIV",
        "condition": "g.br3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 3400, "Kmd": 8300, "kf": 187,
        "adjustments": [("g.br2 == 0", 5), ("g.br2 > 1 or g.br1 > 1", 1.5), ("g.br1 > 0", 0.178)]
    },

    # GnTV
    {
        "name": "GnTV",
        "condition": "g.br1 == 0 and g.br2 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 130, "Kmd": 3500, "kf": 1410,
        "adjustments": [("g.br3 > 0", 0.692)]
    },

    # GnTE, 1
    {
        "name": "GnTE_Br1",
        "condition": "g.br1 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1+2, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 700, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br2 == 0", 4)]
    },

    # GnTE, 2
    {
        "name": "GnTE_Br2",
        "condition": "g.br2 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2+2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 700, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br1 == 0", 4)]
    },

    # GnTE, 3
    {
        "name": "GnTE_Br3",
        "condition": "g.br3 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3+2, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 7000, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br4 == 0", 4)]
    },

    # GnTE, 4
    {
        "name": "GnTE_Br4",
        "condition": "g.br4 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+2, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 7000, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br3 == 0", 4)]
    },

    # GalT, 1
    {
        "name": "GalT_Br1",
        "condition": "g.br1 == 1 or g.br1 == 4",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal+1, g.sia, g.br_o)",
        "cosubstrate": "UDP-Gal",
        "Km": 150, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },

    # GalT, 2
    {
        "name": "GalT_Br2",
        "condition": "g.br2 == 1 or g.br2 == 4",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal+1, g.sia, g.br_o)",
        "cosubstrate": "UDP-Gal",
        "Km": 135, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },

    # GalT, 3
    {
        "name": "GalT_Br3",
        "condition": "g.br3 == 1 or g.br3 == 4",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal+1, g.sia, g.br_o)",
        "cosubstrate": "UDP-Gal",
        "Km": 80, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },

    # GalT, 4
    {
        "name": "GalT_Br4",
        "condition": "g.br4 == 1 or g.br4 == 4",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal+1, g.sia, g.br_o)",
        "cosubstrate": "UDP-Gal",
        "Km": 4000, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },

    # SiaT, 1
    {
        "name": "SiaT_Br1",
        "condition": "g.br1 == 2 or g.br1 == 5",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal, g.sia+1, g.br_o)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },

    # SiaT, 2
    {
        "name": "SiaT_Br2",
        "condition": "g.br2 == 2 or g.br2 == 5",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal, g.sia+1, g.br_o)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },

    # SiaT, 3
    {
        "name": "SiaT_Br3",
        "condition": "g.br3 == 2 or g.br3 == 5",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal, g.sia+1, g.br_o)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },

    # SiaT, 4
    {
        "name": "SiaT_Br4",
        "condition": "g.br4 == 2 or g.br4 == 5",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal, g.sia+1, g.br_o)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },

    # Och1
    {
        "name": "Och1",
        "condition": "g.man1 > 0 and g.man2 > 0 and g.man3 > 0 and g.fuc == 0 and g.gnb == 0 and g.br1 == 0 and g.br2 == 0 and g.br3 == 0 and g.br4 == 0 and g.gal == 0 and g.sia == 0 and g.br_o == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o+1)",
        "cosubstrate": "GDP-Man",
        "Km": 390, "Kmd": 10, "kf": 5.2,
        "adjustments": [("g.man1 + g.man2 + g.man3 == 6", 0.74), ("g.man1 == 2 and g.man2 == 2 and g.man3 == 1", 0.55), ("g.man1 == 1 and g.man2 == 2 and g.man3 == 2", 0.27), ("g.man1 == 2 and g.man2 == 1 and g.man3 == 1", 0.6), ("g.man1 == 1 and g.man2 == 1 and g.man3 == 2", 0.25), ("g.man1 == 1 and g.man2 == 2 and g.man3 == 1", 0.14), ("g.man1 == 1 and g.man2 == 1 and g.man3 == 1", 0.16)]
    },

    # Mnn9
    {
        "name": "Mnn9",
        "condition": "g.br_o == 1",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o+1)",
        "cosubstrate": "GDP-Man",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },
]


n_compartments = 4

def build_enzyme_dist():
    """
    Computes the concentration of each enzyme in every compartment.
    Returns [{enzyme_name: concentration}, ...], length equals number of compartments.
    """
    enzyme_conc = []
    enzyme_names = list(ENZYME_CONC.keys())
    for j in range(n_compartments):
        comp_dict = {}
        for i, enzyme in enumerate(enzyme_names):
            # 隔室浓度 = 总浓度 × 分布比例
            distr = ENZYME_CONC[enzyme][1][j]
            comp_dict[enzyme] = ENZYME_CONC[enzyme][0] * distr
        enzyme_conc.append(comp_dict)
    return enzyme_conc


if __name__ == "__main__": 
    enzymez_dist = build_enzyme_dist()
    for i, comp in enumerate(enzymez_dist):
        print(f"Compartment {i+1}:")
        for enzyme, conc in comp.items():
            print(f"  {enzyme}: {conc:.2f}")