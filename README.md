# Glycosylation Model

A model capable of quantitative predictions of glycosylation results, embedded with a steady-state solver for the Golgi N-glycosylation reaction network and visualization methods.

## Model Overview

The model simulates the enzymatic modification process of N-glycans across 4 cascaded Golgi compartments (CGC, MGC, TGC, TGN). It solves the steady-state mass balance equations to obtain the concentration distribution of all glycoforms.

Steady-state equation:
$$
c_{i,j} = c_{i,j-1} + \tau * r_{i,j}
$$

Reaction rate (modified Michaelis-Menten):
$$
r = \frac{k_f [E_t] [\text{donor}] [G_i] }{ K_{m,i} (K_{m,d} + [\text{donor}])  (1 + \sum{\frac{(G_j)}{K_{m,j}}})}
$$

## File Structure

```
glycosylation-model
├── __init__.py        # Package init
├── glycoform.py       # Defines the glycoform class
├── enzyme.py          # Defines the Enzyme class 
├── reaction.py        # Defines the Reaction class and ReactionNetwork class
├── kinetics.py        # Defines the KineticsCalculator class
├── golgi_model.py     # Defines the GolgiModel class
├── config.py          # Parameter configuration 
├── main.py            # Main program 
├── visualization.py   # Load distribution data and visualize
├── draw_glycan/       # Glycan strucutre renderer
├── test/              # Testing scripts used in [wiki](https://2026.igem.wiki/fudan/model#results)
└── README.md          # This file
```

## How to Run

```bash
cd glycosylation-model
python glycosylation-model/main.py
```

For visualization:

```bash
python visualization.py
```


## Parameter Configuration

All tunable parameters are in `config.py`:


## Output

After running, the following files are generated in the `glycosylation-model/` directory:
- `network_reactions.txt` — List of all reactions
- `network_structures.txt` — All glycoform structures and their participating reactions
- `distribution_data.txt` - Abundance of glycoform structures

For visualization, the following files are additionally generated:

- `top10.png` — Bar chart of the top-10 glycoform normalized abundances
- `top10_structures.svg` — Top-10 glycan structures



## Dependencies

- Python >= 3.10
- numpy
- scipy

For visualization:

- matplotlib


## Reference
[1] Krambeck, F. J., & Betenbaugh, M. J. (2005). A mathematical model of N-linked glycosylation. Biotechnology and bioengineering, 92(6), 711–728. https://doi.org/10.1002/bit.20645    
[2] Mehta, A. Y., & Cummings, R. D. (2020). GlycoGlyph: a glycan visualizing, drawing and naming application. Bioinformatics (Oxford, England), 36(11), 3613–3614. https://doi.org/10.1093/bioinformatics/btaa190

---