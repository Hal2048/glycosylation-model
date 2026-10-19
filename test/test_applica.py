import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend: save only, no window popup
import matplotlib.pyplot as plt

# Add the project root (parent of this file's directory) to the module search path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config as config
from glycoform import Glycoform
from enzyme import build_enzymes
from reaction import ReactionNetwork
from kinetics import KineticsCalculator
from golgi_model import GolgiModel

# Protein production rates to scan
RATES = [i for i in range(50000,200000,100)]

# Scaling factors to scan
S_LIST = [0.2, 0.4, 0.6, 0.8, 1.0]


# ENZYME_CONC stores the unscaled baseline concentrations; back them up first
BASE_CONC = {k: config.ENZYME_CONC[k][0] for k in config.S_SCALED_ENZYMES}


def set_s(new_s):
    """Set s in this script: rewrite the concentrations of the s-scaled
    enzymes to baseline concentration x new_s."""
    for k in config.S_SCALED_ENZYMES:
        config.ENZYME_CONC[k][0] = BASE_CONC[k] * new_s

GRAPH_PATH = Path(__file__).parent / "test3_graph_2.png"  # save path for the line plots
DATA_PATH = GRAPH_PATH.with_name(f"{GRAPH_PATH.stem}_data.txt")          # scan data from main
DATA_PATH_2 = GRAPH_PATH.with_name(f"{GRAPH_PATH.stem}_rate05_vs_s.txt")  # (s, rate_0.5) data from main_2
DATA_PATH_3 = GRAPH_PATH.with_name(f"{GRAPH_PATH.stem}_entropy.txt")      # (s, rate, Shannon entropy) data from main_3


def compute_br_abundance(network, final_concs):
    """Normalized abundances in the last compartment, aggregated by branch
    count (number of nonzero values among br1-br4)."""
    last = final_concs[:, -1]
    total = last.sum()
    norm = last / total if total > 0 else last
    abund = {}
    for idx, g in network.structures.items():
        br_count = sum(1 for b in (g.br1, g.br2, g.br3, g.br4) if b > 0)
        abund[br_count] = abund.get(br_count, 0.0) + float(norm[idx])
    return abund


def compute_shannon_entropy(final_concs):
    """Take the last-compartment normalized abundance p_i per glycoform
    (each structure in the reaction network) and compute the Shannon entropy
    H = -sum p_i*log2(p_i) (base 2, unit: bit).
    Glycoforms with zero abundance do not contribute. Returns 0 if the total
    abundance is 0."""
    last = final_concs[:, -1]
    total = last.sum()
    if total <= 0:
        return 0.0
    p = last / total
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def solve_model(rate):
    """Run the full model once with the given protein production rate,
    returning (network, final_concs)."""
    # compute the total glycan concentration locally for the current rate,
    # without modifying the config module
    tot_glycan_conc = rate * config.tau / config.COMPARTMENT_VOLUME_SINGLE

    # same pipeline as main.py: enzymes -> reaction network -> kinetics ->
    # Golgi model -> solve
    enzymes = build_enzymes(config.ENZYME_RULES)
    network = ReactionNetwork()
    network.generate_network(
        [Glycoform(2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0),  # M9
         Glycoform(2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)],  # M8
        enzymes,
    )

    initial_feed = np.zeros(len(network.structures))
    initial_feed[0] =  tot_glycan_conc * 0.91
    initial_feed[1] =  tot_glycan_conc * 0.09

    model = GolgiModel(
        network=network,
        kinetics=KineticsCalculator(network),
        tau=config.TAU,
        enzyme_dist=config.build_enzyme_dist(),
        donor_conc=config.DONOR_CONC,
        initial_feed=initial_feed,
    )
    final_concs = model.solve_sparse(tol=1e-6, max_iter=50, verbose=False)
    return network, final_concs


def run_one(rate):
    """Run the full model once with the given protein production rate,
    directly returning {branch count: abundance}."""
    network, final_concs = solve_model(rate)
    return compute_br_abundance(network, final_concs)


def run_one_entropy(rate):
    """Run the full model once with the given protein production rate,
    returning the Shannon entropy (bit) computed per glycoform."""
    _, final_concs = solve_model(rate)
    return compute_shannon_entropy(final_concs)


def plot(rates, results):
    """results: {s: [{branch count: abundance} for each rate]}. Draw one
    figure per s; curves with the same branch count share a color, and the s
    value is shown in the legend (format br=X, s=Y).
    X-axis values are divided by 1e4 and the factor is noted in the unit."""
    rates_plot = np.asarray(rates, dtype=float) / 1e4
    br_values = sorted({b for ab in results.values() for a in ab for b in a})
    colors = {b: c for b, c in zip(br_values, plt.cm.tab10.colors)}
    for s, all_abund in results.items():
        fig, ax = plt.subplots(figsize=(8, 6))
        for b in br_values:
            ys = [a.get(b, 0.0) for a in all_abund]
            ax.plot(rates_plot, ys, color=colors[b], label=f"br={b}")
        ax.set_xlabel("Glycoprotein production rate (10⁴ mol/min)", fontsize=12)
        ax.set_ylabel("Abundance", fontsize=12)
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))  # minor X ticks every 0.5 (i.e. 5000)
        # ax.set_title("Glycoform abundance vs protein yield")
        ax.legend(title="Branch count")
        path = GRAPH_PATH.with_name(f"{GRAPH_PATH.stem}_s{s}{GRAPH_PATH.suffix}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Figure saved: {path}")


def save_results(rates, results, path=DATA_PATH):
    """Write the scan results to a tab-separated txt: each line is
    s, rate, abundances per branch count.
    Lines starting with # are comments (including the branch-count column
    names); load_results reads the file directly when replotting."""
    br_values = sorted({b for ab in results.values() for a in ab for b in a})
    with open(path, "w", encoding="utf-8") as f:
        f.write("# main scan data: branch abundances at different protein production rates for each s\n")
        f.write("# column 1: s (enzyme concentration scaling factor)\n")
        f.write("# column 2: protein production rate (mol/cell/min)\n")
        f.write("# columns 3..n: abundances per branch count, column names:\n")
        f.write("# br:\t" + "\t".join(str(b) for b in br_values) + "\n")
        for s, all_abund in results.items():
            for rate, a in zip(rates, all_abund):
                row = [f"{s:.6f}", f"{rate:.6f}"]
                row += [f"{a.get(b, 0.0):.6f}" for b in br_values]
                f.write("\t".join(row) + "\n")
    print(f"Data saved: {path}")


def load_results(path=DATA_PATH):
    """Read (rates, results) back from txt, in the format written by
    save_results."""
    br_values = []
    rows_by_s = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# br:"):
                br_values = [int(x) for x in line[5:].split()]
                continue
            if line.startswith("#"):
                continue
            s, rate, *vals = line.split("\t")
            a = {b: float(v) for b, v in zip(br_values, vals)}
            rows_by_s.setdefault(float(s), []).append((float(rate), a))
    results, rates = {}, None
    for s, rows in rows_by_s.items():
        results[s] = [a for _, a in rows]
        rs = [r for r, _ in rows]
        rates = rs if rates is None else rates
    return rates, results


def replot():
    """Replot main's figure directly from DATA_PATH without rerunning the model."""
    rates, results = load_results()
    plot(rates, results)


def main():
    results = {}
    for s in S_LIST:
        set_s(s)
        results[s] = [run_one(rate) for rate in RATES]
        print(f"s={s} computation done")
    save_results(RATES, results)
    # read back from the just-saved file and plot, which also verifies
    # write/read consistency
    rates, results = load_results()
    plot(rates, results)

def find_rate_0_5(rates, all_abund):
    """Protein production rate at which the br=1 abundance first drops to 0.5
    (linear interpolation between adjacent points).
    Returns None if 0.5 is never crossed."""
    ys = [a.get(1, 0.0) for a in all_abund]
    for i in range(1, len(rates)):
        if (ys[i - 1] - 0.5) * (ys[i] - 0.5) <= 0 and ys[i] != ys[i - 1]:
            return rates[i - 1] + (0.5 - ys[i - 1]) * (rates[i] - rates[i - 1]) / (ys[i] - ys[i - 1])
    return None


def plot_rate05_vs_s(ss, rate_05s):
    """Plot the rate_0.5 vs s curve. Y-axis values are divided by 1e4 and
    the factor is noted in the unit."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(ss, np.asarray(rate_05s, dtype=float) / 1e4, marker="o")
    ax.set_xlabel("s", fontsize=12)
    ax.set_ylabel("rate at glycan diversity (10⁴ mol/cell/min)", fontsize=12)
    path = GRAPH_PATH.with_name(f"{GRAPH_PATH.stem}_rate05_vs_s{GRAPH_PATH.suffix}")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {path}")


def save_rate05_data(ss, rate_05s, path=DATA_PATH_2):
    """Write (s, rate_0.5) to a tab-separated txt; lines starting with # are
    comments."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("# data for the first drop of br=1 abundance to 0.5\n")
        f.write("# column 1: s (enzyme concentration scaling factor)\n")
        f.write("# column 2: protein production rate at br=1 abundance = 0.5 (mol/cell/min)\n")
        for s, r in zip(ss, rate_05s):
            f.write(f"{s:.6f}\t{r:.6f}\n")
    print(f"Data saved: {path}")


def load_rate05_data(path=DATA_PATH_2):
    """Read (s, rate_0.5) back from txt, skipping blank lines and #
    comment lines."""
    ss, rate_05s = [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            a, b = line.split("\t")
            ss.append(float(a))
            rate_05s.append(float(b))
    return ss, rate_05s


def replot_2():
    """Replot main_2's figure directly from DATA_PATH_2 without rerunning the
    model, and run a linear regression, annotating the fit line, the
    function, and R^2 on the figure."""
    ss, rate_05s = load_rate05_data()
    x = np.asarray(ss, dtype=float)
    y = np.asarray(rate_05s, dtype=float) / 1e4  # same unit as the Y axis (10^4 mol/cell/min)

    # linear regression: y = k*x + b
    k, b = np.polyfit(x, y, 1)
    y_fit = k * x + b
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    print(f"Linear fit: y = {k:.4f}*s {'+' if b >= 0 else '-'} {abs(b):.4f}, R^2 = {r_squared:.6f}")

    # plot: raw data + fit line + annotation
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x, y, marker="o", linestyle="", color="#1f77b4")
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, k * x_line + b, "-",color = "#1f77b4")
    ax.set_xlabel("s", fontsize=12)
    ax.set_ylabel("rate at glycan diversity (10⁴ mol/cell/min)", fontsize=12)
    ax.text(
        0.05, 0.95,
        f"y = {k:.4f}·x {'+' if b >= 0 else '−'} {abs(b):.4f}\n$R^2$ = {r_squared:.4f}",
        transform=ax.transAxes, va="top", fontsize=12,
    )
    path = GRAPH_PATH.with_name(f"{GRAPH_PATH.stem}_rate05_vs_s{GRAPH_PATH.suffix}")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {path}")


def main_2():
    ss, rate_05s = [], []
    for s in np.arange(0.001, 0.6, 0.02):
        set_s(s)
        all_abund = [run_one(rate) for rate in RATES]
        r = find_rate_0_5(RATES, all_abund)
        print(f"s={s} computation done, rate_0.5={r}")
        if r is not None:
            ss.append(s)
            rate_05s.append(r)
    save_rate05_data(ss, rate_05s)
    # read back from the just-saved file and plot, which also verifies
    # write/read consistency
    ss, rate_05s = load_rate05_data()
    plot_rate05_vs_s(ss, rate_05s)


def plot_entropy(rates, entropies):
    """entropies: {s: [Shannon entropy for each rate]}. All s values are
    drawn in one figure; curve color varies with s and a colorbar is shown.
    X-axis values are divided by 1e4 and the factor is noted in the unit."""
    rates_plot = np.asarray(rates, dtype=float) / 1e4
    ss = sorted(entropies)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=min(ss), vmax=max(ss)))
    fig, ax = plt.subplots(figsize=(8, 6))
    for s in ss:
        ax.plot(rates_plot, entropies[s], color=sm.to_rgba(s))
    cb = fig.colorbar(sm, ax=ax)
    cb.set_label("s", fontsize=12)
    cb.ax.tick_params(labelsize=12)
    ax.set_xlabel("Glycoprotein production rate (10⁴ mol/min)", fontsize=12)
    ax.set_ylabel("Shannon entropy (bit)", fontsize=12)
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))  # minor X ticks every 0.5 (i.e. 5000)
    path = GRAPH_PATH.with_name(f"{GRAPH_PATH.stem}_entropy{GRAPH_PATH.suffix}")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {path}")


def save_entropy_data(rates, entropies, path=DATA_PATH_3):
    """Write the scan results to a tab-separated txt: each line is
    s, rate, Shannon entropy.
    Lines starting with # are comments; load_entropy_data reads the file
    directly when replotting."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("# main_3 scan data: glycoform Shannon entropy at different protein production rates for each s\n")
        f.write("# column 1: s (enzyme concentration scaling factor)\n")
        f.write("# column 2: protein production rate (mol/cell/min)\n")
        f.write("# column 3: Shannon entropy H = -sum p_i*log2(p_i), computed from the last-compartment normalized abundance per glycoform (bit)\n")
        for s, hs in entropies.items():
            for rate, h in zip(rates, hs):
                f.write(f"{s:.6f}\t{rate:.6f}\t{h:.6f}\n")
    print(f"Data saved: {path}")


def load_entropy_data(path=DATA_PATH_3):
    """Read (rates, entropies) back from txt, in the format written by
    save_entropy_data."""
    rows_by_s = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            s, rate, h = line.split("\t")
            rows_by_s.setdefault(float(s), []).append((float(rate), float(h)))
    entropies = {s: [h for _, h in rows] for s, rows in rows_by_s.items()}
    # the rate column is identical for every s, so take the first s's
    first_s = next(iter(rows_by_s))
    rates = [r for r, _ in rows_by_s[first_s]]
    return rates, entropies


def replot_3():
    """Replot main_3's figure directly from DATA_PATH_3 without rerunning the
    model."""
    rates, entropies = load_entropy_data()
    plot_entropy(rates, entropies)


def main_3():
    entropies = {}
    for s in np.arange(0.001, 0.6, 0.02):
        set_s(s)
        entropies[s] = [run_one_entropy(rate) for rate in RATES]
        print(f"s={s} computation done")
    save_entropy_data(RATES, entropies)
    # read back from the just-saved file and plot, which also verifies
    # write/read consistency
    rates, entropies = load_entropy_data()
    plot_entropy(rates, entropies)


if __name__ == "__main__":
    main()
