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
from config import TAU, DONOR_CONC, build_enzyme_dist


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


def main():
    # ---------- 1. Generate enzyme list ----------
    enzymes = build_enzymes()

    # ---------- 2. Generate reaction network ----------
    network = ReactionNetwork()
    initial_structures = [
        Glycoform(9, 0, 0, 0, 0, 0, 0, 0, 0),
        Glycoform(8, 0, 0, 0, 0, 0, 0, 0, 0),
    ]
    network.generate_network(initial_structures, enzymes)
    n_structs = len(network.structures)
    print(f"Structures: {n_structs}, Reactions: {len(network.reactions)}")

    # Export network to txt
    export_network(network)

    # ---------- 3. Initial feed ----------
    initial_feed = np.zeros(n_structs)
    initial_feed[:2] = 0.5  # M9 and M8 each account for 50%

    # ---------- 4. Initialize KineticsCalculator ----------
    kinetics = KineticsCalculator(network)

    # ---------- 5. Build GolgiModel ----------
    enzyme_dist = build_enzyme_dist()
    model = GolgiModel(
        network=network,
        kinetics=kinetics,
        tau=TAU,
        enzyme_dist=enzyme_dist,
        donor_conc=DONOR_CONC,
        initial_feed=initial_feed,
    )
    print(f"Total variables: {n_structs} x {len(TAU)} = {n_structs * len(TAU)}")

    # ---------- 6. Solve ----------
    # Damped iteration for small-scale debugging;
    # recommend method='krylov' for large-scale production runs
    final_concs = model.solve_damped(tol=1e-6, max_iter=50, damping=0.3)

    # ---------- 7. Result output ----------
    print("\nSolution complete.")
    print(f"Final concentration range: [{final_concs.min():.6f}, {final_concs.max():.4f}]")

    print("\nTotal concentration per compartment:")
    for j in range(len(TAU)):
        total = final_concs[:, j].sum()
        print(f"  Compartment {j}: {total:.2f}")

    last_comp = final_concs[:, -1]
    top_idx = np.argsort(last_comp)[-10:][::-1]
    print("\nTop 10 glycoforms by concentration in TGN:")
    for idx in top_idx:
        if last_comp[idx] > 0.001:
            g = network.structures[idx]
            print(f"  ID {idx:4d}: {g} = {last_comp[idx]:.4f}")

    return final_concs


if __name__ == "__main__":
    final_concs = main()
