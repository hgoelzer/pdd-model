#!/usr/bin/env python3
"""
plot_resolution_comparison.py

Side-by-side bar chart of GrIS SMB for 5 resolutions (CESM2 historical, 1950-1954).
Each year is a group of 5 bars (one per resolution); overlaid line shows the 5-year mean.

Usage (on Betzy, from pdd-model root):
    conda run -n plotting python3 diag/plot_resolution_comparison.py

Output: diag/resolution_comparison.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4 as nc
import glob
import os

# ── Configuration ────────────────────────────────────────────────────────────

BASE   = "/cluster/work/users/heig/pdd"
SETUPS = [
    ("gris_01", 1000,  "1 km"),
    ("gris_02", 2000,  "2 km"),
    ("gris_04", 4000,  "4 km"),
    ("gris_08", 8000,  "8 km"),
    ("gris_16", 16000, "16 km"),
]
YEARS   = list(range(1950, 1955))
PATTERN = "acabf_GIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_{year}.nc"
OUTFILE = os.path.join(os.path.dirname(__file__), "resolution_comparison.pdf")

COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]  # one per resolution

# ── Data loading ─────────────────────────────────────────────────────────────

data = {}   # data[setup_label][year] = SMB Gt/yr
for setup, res, label in SETUPS:
    data[label] = {}
    for year in YEARS:
        path = os.path.join(BASE, setup, "output", PATTERN.format(year=year))
        ds = nc.Dataset(path)
        smb = float(np.sum(ds.variables["acabf"][:]) * res**2 * 30 * 86400 / 1e12)
        ds.close()
        data[label][year] = smb

# ── Layout ───────────────────────────────────────────────────────────────────

labels   = [s[2] for s in SETUPS]
n_res    = len(SETUPS)
n_years  = len(YEARS)
width    = 0.15          # bar width
x        = np.arange(n_years)   # one position per year
offsets  = np.linspace(-(n_res - 1) / 2, (n_res - 1) / 2, n_res) * width

fig, ax = plt.subplots(figsize=(9, 5))

for i, (setup, res, label) in enumerate(SETUPS):
    vals = [data[label][y] for y in YEARS]
    bars = ax.bar(x + offsets[i], vals, width, label=label,
                  color=COLORS[i], alpha=0.85, edgecolor="white", linewidth=0.4)

    # 5-year mean as horizontal line spanning the group
    mean_val = np.mean(vals)
    span = width * n_res / 2
    ax.hlines(mean_val,
              x[0] + offsets[i] - width / 2,
              x[-1] + offsets[i] + width / 2,
              colors=COLORS[i], linewidths=2, linestyles="--", alpha=0.9)

# ── Annotations ──────────────────────────────────────────────────────────────

ax.set_xticks(x)
ax.set_xticklabels([str(y) for y in YEARS])
ax.set_xlabel("Year")
ax.set_ylabel("SMB [Gt yr$^{-1}$]")
ax.set_title("GrIS surface mass balance — resolution comparison\n"
             "CESM2 historical, 1950–1954  |  dashed lines = 5-year mean")
ax.legend(title="Resolution", loc="upper left", framealpha=0.8)
ax.axhline(0, color="black", linewidth=0.6, linestyle="-")
ax.yaxis.grid(True, linestyle=":", alpha=0.5)
ax.set_axisbelow(True)

# Mean values text box
means_str = "5-yr means [Gt/yr]:\n" + "\n".join(
    "%s: %.1f" % (s[2], np.mean([data[s[2]][y] for y in YEARS])) for s in SETUPS
)
ax.text(0.98, 0.04, means_str, transform=ax.transAxes,
        ha="right", va="bottom", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.8))

fig.tight_layout()
fig.savefig(OUTFILE, dpi=150)
print("Saved:", OUTFILE)
