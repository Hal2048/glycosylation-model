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

# Total enzyme names corresponding to DIST_MATRIX columns
_BASE_ENZYME_NAMES = [
    "ManI", "ManII", "FucT", "GnTI", "GnTII",
    "GnTIII", "GnTIV", "GnTV", "GalT", "SiaT",
]

# Maps 20 specific enzyme names in the network to 10 total enzyme groups
_ENZYME_GROUP_MAP = {
    "ManI": "ManI", "ManII": "ManII", "FucT": "FucT",
    "GnTI": "GnTI", "GnTII": "GnTII", "GnTIII": "GnTIII",
    "GnTIV": "GnTIV", "GnTV": "GnTV",
    "GnTE_Br1": "GnTIII", "GnTE_Br2": "GnTIII",
    "GnTE_Br3": "GnTIII", "GnTE_Br4": "GnTIII",
    "GalT_Br1": "GalT", "GalT_Br2": "GalT",
    "GalT_Br3": "GalT", "GalT_Br4": "GalT",
    "SiaT_Br1": "SiaT", "SiaT_Br2": "SiaT",
    "SiaT_Br3": "SiaT", "SiaT_Br4": "SiaT",
}


def build_enzyme_dist():
    """
    Computes the concentration of each enzyme in every compartment.
    Returns [{enzyme_name: concentration}, ...], length equals number of compartments.
    """
    n_compartments = len(TAU)
    enzyme_dist = []
    for j in range(n_compartments):
        # First compute compartment concentrations for 10 total enzyme groups
        base_dict = {}
        for idx, base_name in enumerate(_BASE_ENZYME_NAMES):
            base_dict[base_name] = TOTAL_ENZYME_CONC[base_name] * DIST_MATRIX[j, idx]

        # Expand to 20 specific enzyme names
        full_dict = {}
        for full_name, base_name in _ENZYME_GROUP_MAP.items():
            full_dict[full_name] = base_dict.get(base_name, 0.0)

        enzyme_dist.append(full_dict)
    return enzyme_dist
