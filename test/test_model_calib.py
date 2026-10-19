"""
Model parameter calibration.

Adjustable parameters:
  - total concentration and 4-compartment distribution fraction of each
    enzyme group in FIT_ENZYME_GROUPS (all enzymes in the same group share
    one set of parameters, e.g. multiple ManI variants)
  - concentration of each donor substrate in FIT_DONORS

Fitting objective:
  make the model's glycoform-class abundances in the last compartment
  (aggregated after classify() classification) approach the user-given
  REFERENCE_ABUNDANCE.

Optimization method:
  numerical gradient descent (central-difference gradient + backtracking
  line search for the step size); parameters are log-transformed to keep
  them positive, and distribution parameters are normalized via softmax.
"""

import contextlib
import sys
from pathlib import Path

# Add the project root (parent of this file's directory) to the module search path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from glycoform import Glycoform
import numpy as np
from main import run_model
from config import ENZYME_CONC, DONOR_CONC


# Enzymes to fit, organized in "parameter-sharing groups": all enzymes in
# the same sublist share one set of total concentration and distribution
# fractions. For _pack/_unpack, one group = 5 parameters
# (1 log concentration + 4 log distribution weights); every enzyme in the
# group gets the same values. For example, ManI_9 / ManI_8 are essentially
# the same ManI acting on different substrates; putting them in one group
# lets the optimizer adjust them together.
FIT_ENZYME_GROUPS = [
    ["Och1"],
    ["Mnn9"]
]

# Donor substrate concentrations to fit

FIT_DONORS = ["GDP-Man"]

# Reference abundances: {class name: relative abundance}
REFERENCE_ABUNDANCE = {
    "M9": 0.0476,
    "M10": 0.2857,
    "M11": 0.6667,
}

# Fitting hyperparameters
MAX_ITERS = 1000        # maximum number of outer iterations
LEARNING_RATE = 0.5     # initial step size (decayed on backtracking failure)
GRAD_EPS = 1e-3         # relative step size for the numerical gradient
ATOL = 1e-6             # absolute loss threshold: stop once loss < ATOL ("stop when good enough")


def classify(g):
    """Classify by total mannose count
    """
    n_man = 3 + g.man1 + g.man2 + g.man3 + g.br_o
    return f"M{n_man}"



def _pack(raw):
    """Parameter vector -> (enzyme_conc override dict, donor_conc override dict).

    Each group's shared parameters occupy 5 consecutive components of the
    vector; every enzyme in the group is written the same [conc, dist], so
    run_model assigns them identical compartment concentrations.
    """
    enzymes, donors, i = {}, {}, 0
    for group in FIT_ENZYME_GROUPS:
        conc = float(np.exp(raw[i])); i += 1
        log_w = raw[i:i + 4]; i += 4
        w = np.exp(log_w - log_w.max())
        dist = (w / w.sum()).tolist()
        for name in group:
            enzymes[name] = [conc, dist]
    for name in FIT_DONORS:
        donors[name] = float(np.exp(raw[i])); i += 1
    return enzymes, donors


def _unpack():
    """Build the initial parameter vector (log space) from the current
    defaults in config.

    Each group only reads the config value of its first enzyme as the
    initial value of the group's shared parameters; if the enzymes in a
    group have different config defaults, they are silently unified to the
    first one's values.
    """
    raw = []
    for group in FIT_ENZYME_GROUPS:
        conc, dist = ENZYME_CONC[group[0]]
        raw.append(np.log(conc))
        d = np.clip(np.asarray(dist, dtype=float), 1e-8, None)
        raw.extend(np.log(d / d.sum()))
    for name in FIT_DONORS:
        raw.append(np.log(DONOR_CONC[name]))
    return np.asarray(raw)


# --------------------------------------------------------------------
# Loss: squared error between predicted class abundances and reference
# --------------------------------------------------------------------

def _class_abundances(network, final_concs):
    """Normalized abundances in the last compartment, aggregated by the
    classes of classify()."""
    last = final_concs[:, -1]
    total = last.sum()
    norm = last / total if total > 0 else last
    out = {}
    for idx, g in network.structures.items():
        cls = classify(g)
        out[cls] = out.get(cls, 0.0) + float(norm[idx])
    return out


def _loss(pred, ref):
    """Squared error; ref is first normalized to a probability distribution,
    and classes missing from pred are counted as 0."""
    keys = sorted(ref.keys())
    p = np.array([pred.get(k, 0.0) for k in keys])
    r = np.asarray([ref[k] for k in keys], dtype=float)
    r = r / r.sum()
    return float(np.sum((p - r) ** 2))


def _format_distribution(abundances):
    """Sort {class: abundance} by the number after M and format as
    'M8:  5.00%  M9: 16.67% ...'."""
    def _key(item):
        name = item[0]
        if isinstance(name, str) and name.startswith("M") and name[1:].isdigit():
            return (0, int(name[1:]))
        return (1, str(name))
    return "  ".join(
        f"{k}: {v * 100:6.2f}%"
        for k, v in sorted(abundances.items(), key=_key)
    )


# --------------------------------------------------------------------
# Numerical gradient descent
# --------------------------------------------------------------------

def _numerical_gradient(loss_fn, x):
    """Central-difference gradient."""
    g = np.empty_like(x)
    for i in range(len(x)):
        h = GRAD_EPS * max(abs(x[i]), 1.0)
        xp = x.copy(); xp[i] += h
        xm = x.copy(); xm[i] -= h
        g[i] = (loss_fn(xp) - loss_fn(xm)) / (2 * h)
    return g


def _gradient_descent(loss_fn, x0, state):
    """Simple gradient descent + backtracking line search, returning
    (best parameters, best loss).

    state is a dict shared with loss_fn: after each solve, loss_fn writes
    the class abundances into state['abundances'], which is printed here
    at iter 0 and after every accepted new parameter set.
    """
    x, f = x0.copy(), loss_fn(x0)
    print(f"[fit] iter 0  loss = {f:.6e}")
    print(f"       pred: {_format_distribution(state['abundances'])}")
    if f < ATOL:
        print("[fit] converged (loss at iter 0 already below ATOL)")
        return x, f
    for it in range(1, MAX_ITERS + 1):
        g = _numerical_gradient(loss_fn, x)
        if not np.all(np.isfinite(g)) or np.linalg.norm(g) < 1e-14:
            print("[fit] gradient degenerated, stopping")
            break
        # backtracking line search: until a point with lower loss is found
        step, x_new, f_new = LEARNING_RATE, None, None
        for _ in range(25):
            cand = x - step * g
            f_cand = loss_fn(cand)
            if np.isfinite(f_cand) and f_cand < f:
                x_new, f_new = cand, f_cand
                break
            step *= 0.5
        if x_new is None:
            print("[fit] line search failed, stopping")
            break
        rel = (f - f_new) / max(f, 1e-12)
        x, f = x_new, f_new
        print(f"[fit] iter {it:3d}  loss = {f:.6e}  step = {step:.2e}  rel_improve = {rel:.3e}")
        print(f"       pred: {_format_distribution(state['abundances'])}")
        if f < ATOL:
            print("[fit] converged (loss below ATOL)")
            break
    return x, f


def main():
    if not REFERENCE_ABUNDANCE:
        raise ValueError("Please fill in the reference abundances for each class in REFERENCE_ABUNDANCE first")

    initial_structures = [
        Glycoform(2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        Glycoform(2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    ]

    state = {"abundances": None}

    def loss_fn(x):
        enz_over, don_over = _pack(x)
        network, final_concs = run_model(
            initial_structures=initial_structures, enzyme_conc=enz_over, donor_conc=don_over, verbose=False
        )
        ab = _class_abundances(network, final_concs)
        state["abundances"] = ab
        return _loss(ab, REFERENCE_ABUNDANCE)

    x0 = _unpack()
    ref_norm = {k: v / sum(REFERENCE_ABUNDANCE.values())
                for k, v in REFERENCE_ABUNDANCE.items()}
    print(f"[fit] number of parameters to fit: {len(x0)}")
    print(f"[fit] reference distribution: {_format_distribution(ref_norm)}")
    x_best, f_best = _gradient_descent(loss_fn, x0, state)

    # Re-evaluate the best parameters once, to make sure state holds the
    # distribution at x_best (if we exited due to backtracking failure,
    # state may still hold the distribution of some rejected candidate)
    loss_fn(x_best)

    # output the best parameters
    enz_best, don_best = _pack(x_best)
    print("\n[fit] best parameters:")
    for group in FIT_ENZYME_GROUPS:
        conc, dist = enz_best[group[0]]
        names = "/".join(group)
        print(f"  {names:18s}  conc = {conc:.4f}   dist = {[round(d, 4) for d in dist]}")
    for name in FIT_DONORS:
        print(f"  {name:12s} conc = {don_best[name]:.4f}")
    print(f"[fit] best distribution: {_format_distribution(state['abundances'])}")
    print(f"[fit] final loss = {f_best:.6e}")


class _Tee:
    """Output stream that writes to both stdout and a file."""
    def __init__(self, fp):
        self._fp = fp
        self._out = sys.stdout

    def write(self, s):
        self._out.write(s)
        self._fp.write(s)

    def flush(self):
        self._out.flush()
        self._fp.flush()


if __name__ == "__main__":
    log_path = Path(__file__).parent / "fit_log_test1.txt"
    with open(log_path, "w", encoding="utf-8") as fp, contextlib.redirect_stdout(_Tee(fp)):
        main()
