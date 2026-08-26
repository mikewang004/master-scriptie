from analyse9 import polymer, atom_coords, get_time_temp_from_slurm
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from simulation import Simulation, fit_functions
import scienceplots
from matplotlib.lines import Line2D
import pandas as pd


plt.style.use('science')
plt_cm_to_in = 1/2.54
plt_max_x = 15.5
plt_max_y = 22
plt_caption_font = 9 #pt

plt_colours_chain_length = "viridis"
plt_colours_time = "cividis"




# Time grid
short_time = np.linspace(0, 24_000, 500)
long_time = np.linspace(24_000, 1e8, 1500)
time = np.unique(np.concatenate([short_time, long_time]))

# Temperature schedule:
# linear cooling from 1.00 to 0.88 by t = 24,000,
# then remains constant at 0.88
temperature = np.where(
    time <= 24_000,
    1.0 + (0.88 - 1.0) * time / 24_000,
    0.88,
)

# Plot
fig, ax = plt.subplots(figsize=(8.5, 4.8))

ax.plot(time, temperature, color="tab:blue", linewidth=2)
ax.scatter([0, 24_000], [1.0, 0.88], color="tab:blue", zorder=3)

ax.axvline(24_000, color="0.45", linestyle="--", linewidth=1)
ax.axhline(0.88, color="0.45", linestyle="--", linewidth=1)

ax.text(
    12_000, 0.99, "Quick quench",
    ha="center", va="center",
    fontsize=12
)

ax.text(
    62_000, 0.99, "Isothermal crystallisation",
    ha="center", va="center",
    fontsize=12
)

ax.set_xlim(-10000, 1e5)
ax.set_ylim(0.86, 1.02)
ax.set_xlabel(r"$t/\tau$")
ax.set_ylabel(r"T")
ax.set_title("Overview of quench procedure")

#ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
#ax.grid(True, alpha=0.25)

#plt.tight_layout()
plt.savefig("quench_overview.pdf")
plt.show()