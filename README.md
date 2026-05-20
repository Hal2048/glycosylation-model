# Glycoform Reaction Network Model

Steady-state solver for the Golgi N-glycosylation reaction network based on the modified Michaelis-Menten equation.

## Model Overview

This model describes the enzymatic modification of N-glycans across 4 cascaded CSTR compartments (CGC, MGC, TGC, TGN). It solves the steady-state mass balance equations to obtain the concentration distribution of all glycoforms.

Steady-state equation:
```
c_{i,j} = c_{i,j-1} + tau_j * r_{i,j}
```

Reaction rate (modified Michaelis-Menten):
```
r = kf * [Et] * [donor] * [Pi] / (Kmi * (Kmd + [donor]) * (1 + sum(Pj)/Kmj))
```

## File Structure

```
model/
├── __init__.py        # Package init
├── glycoform.py       # Glycoform class (9-digit encoded glycan structure)
├── enzyme.py          # Enzyme class + 20 enzyme rules + build_enzymes()
├── reaction.py        # Reaction class + ReactionNetwork class
├── kinetics.py        # KineticsCalculator (competitive sum 1+sum(c/Km))
├── golgi_model.py     # GolgiModel (residual function + steady-state solver)
├── config.py          # Parameter configuration (tau, enzyme dist, donor conc, etc.)
├── main.py            # Main program entry point
└── README.md          # This file
```

## How to Run

### Option 1: Direct execution (recommended)

```bash
cd Model_versions
python model/main.py
```

### Option 2: As a module

```bash
cd Model_versions
python -m model.main
```

## Parameter Configuration

All tunable parameters are centralized in `config.py`:

| Parameter | Description | Current Value |
|-----------|-------------|---------------|
| `TAU` | Residence time per compartment (s) | [5.56, 5.56, 5.56, 5.56] |
| `DONOR_CONC` | Donor concentrations (uM) | See config.py |
| `TOTAL_ENZYME_CONC` | Total enzyme concentrations (uM) | See config.py |
| `DIST_MATRIX` | Enzyme compartment distribution ratio | 4x10 matrix |

## Output

After running, the following files are generated in the `model/` directory:

- `network_reactions.txt` — List of all reactions
- `network_structures.txt` — All glycoform structures and their participating reactions

Console output includes:
- Number of structures / reactions / total variables
- Iteration convergence info
- Total concentration per compartment
- Top 10 glycoforms by concentration in TGN

## Solver Selection

`GolgiModel` provides two solving methods:

1. **`solve_damped()`** — Damped Picard iteration
   - No Jacobian computation, memory-friendly
   - Suitable for debugging and small-scale testing

2. **`solve(method='krylov')`** — Krylov subspace method
   - Suitable for large-scale sparse problems (~30K variables)
   - Recommended for production calculations

Switch in `main.py`:
```python
# Damped iteration (debugging)
final_concs = model.solve_damped(tol=1e-6, max_iter=50, damping=0.3)

# Krylov method (large-scale)
final_concs = model.solve(method='krylov', tol=1e-6, max_iter=100)
```

## Dependencies

- Python >= 3.8
- numpy
- scipy
