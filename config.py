import numpy as np


# Residence time for each compartment (CGC, MGC, TGC, TGN)
TAU = np.array([5.56, 5.56, 5.56, 5.56])

# Donor concentrations (uM)
DONOR_CONC = {
    "UDP-GlcNAc": 9143.0,
    "UDP-Gal": 3810.0,
    "CMP-NeuAc": 2286.0,
    "GDP-Fuc": 5000.0,
    "UDP-GalNAc": 3000.0,
}

# Total enzyme concentrations (uM)
TOTAL_ENZYME_CONC = {
    "ManI": 24.8,
    "ManII": 36.4,
    "FucT": 14.3,
    "GnTI": 29.7,
    "GnTII": 18.2,
    "GnTIII": 111.0,
    "GnTIV": 18.7,
    "GnTV": 21.4,
    "GalT": 87.1,
    "SiaT": 48.4,
}

# Compartment distribution ratio (4 compartments x 10 total enzyme groups)
DIST_MATRIX = np.array([
    [0.05, 0.20, 0.05, 0.15, 0.10, 0.10, 0.10, 0.05, 0.05, 0.05],  # CGC
    [0.15, 0.40, 0.15, 0.30, 0.25, 0.20, 0.25, 0.15, 0.20, 0.15],  # MGC
    [0.40, 0.25, 0.40, 0.35, 0.40, 0.35, 0.30, 0.40, 0.35, 0.40],  # TGC
    [0.40, 0.15, 0.40, 0.20, 0.25, 0.35, 0.35, 0.40, 0.40, 0.40],  # TGN
])


n_compartments = 4

def build_enzyme_dist():
    """
    Computes the concentration of each enzyme in every compartment.
    Returns [{enzyme_name: concentration}, ...], length equals number of compartments.
    """
    enzyme_conc = []
    enzyme_names = list(TOTAL_ENZYME_CONC.keys())
    for j in range(n_compartments):
        comp_dict = {}
        for i, enzyme in enumerate(enzyme_names):
            # 隔室浓度 = 总浓度 × 分布比例
            comp_dict[enzyme] = TOTAL_ENZYME_CONC[enzyme] * DIST_MATRIX[j, i]
        enzyme_conc.append(comp_dict)
    return enzyme_conc