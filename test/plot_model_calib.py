"""
Plot the model calibration results
"""

import numpy as np
import matplotlib.pyplot as plt

# The three data series (percentages)
categories = ["M9", "M10", "M11"]
initial    = [72.15, 27.81, 0.04]  
optimized  = [7.10, 27.41, 65.49]
reference  = [4.76, 28.57, 66.67]

# Grouped bar chart
x = np.arange(len(categories))
width = 0.25

fig, ax = plt.subplots(figsize=(6, 4))
bars1 = ax.bar(x - width, initial,   width, label="Initial", color="lightgray")
bars2 = ax.bar(x,         optimized, width, label="Optimized", color="lightblue")
bars3 = ax.bar(x + width, reference, width, label="Reference", color="violet")

# Label the value above each bar
for bars in (bars1, bars2, bars3):
    ax.bar_label(bars, fmt="%.2f%%", fontsize=8, padding=2)

ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.set_ylabel("Abundance (%)")
ax.set_title("")
ax.legend()
ax.set_ylim(0, 100)  # leave headroom for the labels at the top

fig.tight_layout()
fig.savefig("calib_result_test1.png", dpi=200)
