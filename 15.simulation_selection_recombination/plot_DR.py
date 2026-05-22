#!/usr/bin/env python3
"""
plot_DR.py -- Two panels: D' and R2 vs recombination suppression fold
Legend: s=xx (weak/moderate/strong), no 2Ns values.
Ordered weak → strong.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import OrderedDict

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7,
    "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.4, "ytick.major.width": 0.4,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "legend.fontsize": 5.5, "legend.framealpha": 0.9,
})

MM = 1 / 25.4
TARGET_D = 0.963
TARGET_R = 0.595

# ==================================================================
# DATA: ordered weak → strong (legend simplified)
# ==================================================================
all_data = OrderedDict([
    ("s=6.5e-4 (weak)", {
        "color": "#2196F3", "marker": "^",
        "folds":  np.array([30, 100, 300, 1000, 10000]),
        "D_mean": np.array([0.744, 0.828, 0.885, 0.941, 0.969]),
        "D_std":  np.array([0.055, 0.051, 0.062, 0.042, 0.021]),
        "R_mean": np.array([0.270, 0.344, 0.386, 0.445, 0.451]),
        "R_std":  np.array([0.051, 0.055, 0.079, 0.110, 0.132]),
    }),
    ("s=3.2e-3 (moderate)", {
        "color": "#FF9800", "marker": "v",
        "folds":  np.array([30, 100, 300, 1000, 10000]),
        "D_mean": np.array([0.732, 0.872, 0.918, 0.971, 0.991]),
        "D_std":  np.array([0.067, 0.062, 0.053, 0.023, 0.006]),
        "R_mean": np.array([0.256, 0.423, 0.446, 0.532, 0.579]),
        "R_std":  np.array([0.045, 0.094, 0.109, 0.145, 0.077]),
    }),
    ("s=5.8e-3 (strong)", {
        "color": "#4CAF50", "marker": "D",
        "folds":  np.array([10, 30, 50, 70, 100, 300, 1000, 3000, 10000]),
        "D_mean": np.array([0.718, 0.811, 0.787, 0.870, 0.861,
                             0.936, 0.977, 0.983, 0.994]),
        "D_std":  np.array([0.094, 0.059, 0.089, 0.054, 0.078,
                             0.035, 0.022, 0.024, 0.012]),
        "R_mean": np.array([0.271, 0.357, 0.337, 0.456, 0.402,
                             0.537, 0.565, 0.598, 0.696]),
        "R_std":  np.array([0.076, 0.090, 0.127, 0.055, 0.100,
                             0.061, 0.139, 0.138, 0.156]),
    }),
    ("s=2.0e-2 (v.strong)", {
        "color": "#E91E63", "marker": "p",
        "folds":  np.array([30, 100, 300, 1000, 10000]),
        "D_mean": np.array([0.907, 0.936, 0.982, 0.985, 0.998]),
        "D_std":  np.array([0.045, 0.045, 0.023, 0.024, 0.004]),
        "R_mean": np.array([0.529, 0.614, 0.773, 0.702, 0.781]),
        "R_std":  np.array([0.110, 0.124, 0.134, 0.129, 0.180]),
    }),
])

# ==================================================================
# FIGURE: 1 row × 2 cols
# ==================================================================
fig, (ax_d, ax_r) = plt.subplots(1, 2, figsize=(150 * MM, 72 * MM),
                                  gridspec_kw={"wspace": 0.35})
fig.subplots_adjust(left=0.09, right=0.98, top=0.82, bottom=0.16)

# ============================================================
# Panel A: D' vs fold
# ============================================================
for label, d in all_data.items():
    ax_d.errorbar(d["folds"], d["D_mean"], yerr=d["D_std"],
                  fmt=f"{d['marker']}-", color=d["color"],
                  ms=4.5, lw=1.1,
                  capsize=2.5, capthick=0.5, elinewidth=0.5,
                  label=label, zorder=4)

ax_d.axhline(y=TARGET_D, color="gold", ls="--", lw=1.5,
              label=f"Empirical D'={TARGET_D}", zorder=2)
ax_d.axhspan(TARGET_D - 0.03, TARGET_D + 0.03,
              color="gold", alpha=0.12, zorder=1)

ax_d.set_xscale("log")
ax_d.set_xlabel("Recombination suppression fold")
ax_d.set_ylabel("D' (mean +/- std, n=10 reps)")
ax_d.set_title("A  D' vs recombination suppression",
                fontweight="bold", loc="left", fontsize=8)
ax_d.set_ylim(0.55, 1.05)
ax_d.set_xlim(7, 15000)
ax_d.set_xticks([10, 30, 100, 300, 1000, 3000, 10000])
ax_d.set_xticklabels(["10x", "30", "100", "300", "1K", "3K", "10K"])
ax_d.legend(loc="lower right", frameon=True, ncol=1,
             borderpad=0.3, handlelength=1.5)

# ============================================================
# Panel B: R2 vs fold
# ============================================================
for label, d in all_data.items():
    ax_r.errorbar(d["folds"], d["R_mean"], yerr=d["R_std"],
                  fmt=f"{d['marker']}-", color=d["color"],
                  ms=4.5, lw=1.1,
                  capsize=2.5, capthick=0.5, elinewidth=0.5,
                  label=label, zorder=4)

ax_r.axhline(y=TARGET_R, color="gold", ls="--", lw=1.5,
              label=f"Empirical $R^2$={TARGET_R}", zorder=2)
ax_r.axhspan(TARGET_R - 0.05, TARGET_R + 0.05,
              color="gold", alpha=0.12, zorder=1)

ax_r.set_xscale("log")
ax_r.set_xlabel("Recombination suppression fold")
ax_r.set_ylabel(r"$R^2$ (mean +/- std, n=10 reps)")
ax_r.set_title(r"B  $R^2$ vs recombination suppression",
                fontweight="bold", loc="left", fontsize=8)
ax_r.set_ylim(0.05, 1.05)
ax_r.set_xlim(7, 15000)
ax_r.set_xticks([10, 30, 100, 300, 1000, 3000, 10000])
ax_r.set_xticklabels(["10x", "30", "100", "300", "1K", "3K", "10K"])
ax_r.legend(loc="upper left", frameon=True, ncol=1,
             borderpad=0.3, handlelength=1.5)

# ---- Suptitle ----
fig.suptitle(
    "Coalescent simulations (10 reps each): "
    "selection strength vs recombination suppression\n"
    f"NE_LW=500,000 | Grid: 10 Mb SV region | "
    f"Target: D'={TARGET_D}, $R^2$={TARGET_R} | "
    "Error bars = 1 std",
    fontsize=6, fontweight="bold", y=0.96)

plt.savefig("grid_reps_DR.pdf", dpi=300, bbox_inches="tight")
plt.savefig("grid_reps_DR.png", dpi=300, bbox_inches="tight")
print("Saved: grid_reps_DR.pdf / grid_reps_DR.png")
