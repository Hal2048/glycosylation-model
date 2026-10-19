import sys
from pathlib import Path

# Add the directory containing this file to the module search path
# to ensure it can run from any location
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

from glycoform import Glycoform
from enzyme import build_enzymes
from reaction import ReactionNetwork
from kinetics import KineticsCalculator
from golgi_model import GolgiModel
from config import TAU, DONOR_CONC, ENZYME_CONC, build_enzyme_dist, ENZYME_RULES, TOT_GLYCAN_CONC


GLYCOFORM_FIELDS = ("man1", "man2", "man3", "fuc", "gnb",
                    "br1", "br2", "br3", "br4", "gal", "sia", "br_o")


def export_network(network: ReactionNetwork, out_dir: Path = None):
    """Exports structure and reaction information of the network to txt files."""
    if out_dir is None:
        out_dir = Path(__file__).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Export reaction list
    reactions_file = out_dir / "network_reactions.txt"
    with open(reactions_file, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("Reaction Network - All Reactions\n")
        f.write(f"Total Reactions: {len(network.reactions)}\n")
        f.write("=" * 80 + "\n\n")
        for index, reaction in network.reactions.items():
            line = (
                f"Reaction ID: {reaction.id}, "
                f"Enzyme: {reaction.enzyme.enzyme_name}, "
                f"Substrate ID: {reaction.substrate_id}, "
                f"Product ID: {reaction.product_id}, "
                f"Km: {reaction.km_value}\n"
            )
            f.write(line)
    print(f"Reaction list exported: {reactions_file}")

    # Export structure list
    structures_file = out_dir / "network_structures.txt"
    with open(structures_file, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("Reaction Network - All Structures\n")
        f.write(f"Total Structures: {len(network.structures)}\n")
        f.write("=" * 80 + "\n\n")
        for index, structure in network.structures.items():
            line = (
                f"Structure ID: {index}, "
                f"Glycoform: {structure}, "
                f"Consuming Reactions: {structure._consuming_reactions}, "
                f"Producing Reactions: {structure._producing_reactions}\n"
            )
            f.write(line)
    print(f"Structure list exported: {structures_file}")


def export_distribution(
    network: ReactionNetwork,
    final_concs: np.ndarray,
    out_dir: Path = None,
    compartment_idx: int = -1,
):
    """
    Exports the steady-state glycoform distribution to a tab-separated
    ``distribution_data.txt`` file.
    """
    if out_dir is None:
        out_dir = Path(__file__).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    col_conc = final_concs[:, compartment_idx]
    total = col_conc.sum()
    abundances = col_conc / total if total > 0 else col_conc

    order = np.argsort(-col_conc)  # descending by raw concentration
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(order) + 1)

    dist_file = out_dir / "distribution_data.txt"
    with open(dist_file, "w", encoding="utf-8") as f:
        f.write("# Glycoform Distribution\n")
        f.write(f"# Compartment Index: {compartment_idx}\n")
        f.write(f"# Total Structures: {len(network.structures)}\n")
        f.write(f"# Total Concentration: {total:.6e}\n")
        f.write(
            "# Columns: Rank\tID\tGlycoform Tuple ("
            + ", ".join(GLYCOFORM_FIELDS)
            + ")\tRaw Concentration\tNormalized Abundance\n"
        )
        for idx in order:
            glyco = network.structures[int(idx)]
            tup = tuple(getattr(glyco, name) for name in GLYCOFORM_FIELDS)
            f.write(
                f"{ranks[idx]}\t{int(idx)}\t{tup}\t{col_conc[idx]:.6e}\t"
                f"{abundances[idx]:.6e}\n"
            )
    print(f"Distribution exported: {dist_file}")



def run_model(initial_structures, enzyme_conc=None, donor_conc=None,
              enzyme_rules=None, verbose=False):
    """
    Build and solve the model, returning (network, final_concs)
    """
    eff_enzyme_conc = ENZYME_CONC if enzyme_conc is None else {**ENZYME_CONC, **enzyme_conc}
    eff_donor = DONOR_CONC if donor_conc is None else {**DONOR_CONC, **donor_conc}
    rules = ENZYME_RULES if enzyme_rules is None else enzyme_rules

    # ---------- Generate enzyme list ----------
    enzymes = build_enzymes(rules)

    # ---------- Generate reaction network ----------
    network = ReactionNetwork()

    network.generate_network(initial_structures, enzymes)

    initial_feed = np.zeros(len(network.structures))
    # initial_feed[0] =  TOT_GLYCAN_CONC * 0.91
    # initial_feed[1] =  TOT_GLYCAN_CONC * 0.09
    initial_feed[0] =  TOT_GLYCAN_CONC * 0.5
    initial_feed[1] =  TOT_GLYCAN_CONC * 0.5

    # ---------- Initialize KineticsCalculator ----------
    kinetics = KineticsCalculator(network)

    # ---------- Build GolgiModel ----------
    enzyme_dist = build_enzyme_dist(eff_enzyme_conc)
    model = GolgiModel(
        network=network,
        kinetics=kinetics,
        tau=TAU,
        enzyme_dist=enzyme_dist,
        donor_conc=eff_donor,
        initial_feed=initial_feed,
    )

    # ---------- Solve ----------
    # Damped iteration for small-scale debugging;
    # recommend method='krylov' for large-scale production runs
    final_concs = model.solve_sparse(tol=1e-6, max_iter=50, verbose=verbose)
    return network, final_concs


def main():
    initial_structures = [
        Glycoform(2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        Glycoform(2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    ]

    network, final_concs = run_model(initial_structures=initial_structures, verbose=True)
    n_structs = len(network.structures)
    print(f"Structures: {n_structs}, Reactions: {len(network.reactions)}")

    # ---------- 7. Export network & result output ----------
    export_network(network)

    print("\nSolution complete.")
    print(f"Final concentration range: [{final_concs.min():.6f}, {final_concs.max():.4f}]")

    print("\nTotal concentration per compartment:")
    for j in range(len(TAU)):
        total = final_concs[:, j].sum()
        print(f"  Compartment {j}: {total:.2f}")

    last_comp = final_concs[:, -1]
    last_comp_norm = last_comp / last_comp.sum()
    top_idx = np.argsort(last_comp_norm)[-10:][::-1]

    print("\nTop 10 glycoforms by concentration in TGN:")
    for idx in top_idx:
        g = network.structures[idx]
        print(f"  ID {idx:4d}: {g} = {last_comp_norm[idx]:.6f} (raw={last_comp[idx]:.4f})")

    # ---------- 8. Export glycoform distribution ----------
    export_distribution(network, final_concs, compartment_idx=-1)

    return final_concs


if __name__ == "__main__":
    final_concs = main()
