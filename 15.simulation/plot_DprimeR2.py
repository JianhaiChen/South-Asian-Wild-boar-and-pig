#!/usr/bin/env python3
"""
plot_DR_final.py -- Two panels: D' and R2 vs recombination suppression fold
Legend: s=xx (weak/moderate/strong). Ordered weak → strong.
All curves aligned to 5 data points [30, 100, 300, 1000, 10000] using the final fixed-RNG data.
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

# 设置全局字体和排版参数
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
# 最新修复随机抽样 (fixed RNG) 后的 10 reps 数据
# ==================================================================
all_data = OrderedDict([
    ("s=6.5e-4 (weak)", {
        "color": "#2196F3", "marker": "^",
        "folds":  np.array([30, 100, 300, 1000, 10000]),
        "D_mean": np.array([0.736, 0.825, 0.890, 0.942, 0.970]),
        "D_std":  np.array([0.051, 0.055, 0.059, 0.039, 0.019]),
        "R_mean": np.array([0.265, 0.347, 0.392, 0.444, 0.451]),
        "R_std":  np.array([0.050, 0.058, 0.071, 0.120, 0.134]),
    }),
    ("s=3.2e-3 (moderate)", {
        "color": "#FF9800", "marker": "v",
        "folds":  np.array([30, 100, 300, 1000, 10000]),
        "D_mean": np.array([0.737, 0.870, 0.916, 0.972, 0.992]),
        "D_std":  np.array([0.064, 0.065, 0.059, 0.023, 0.006]),
        "R_mean": np.array([0.256, 0.416, 0.443, 0.529, 0.587]),
        "R_std":  np.array([0.045, 0.092, 0.105, 0.142, 0.073]),
    }),
    ("s=5.8e-3 (strong)", {
        "color": "#4CAF50", "marker": "D",
        "folds":  np.array([30, 100, 300, 1000, 10000]),
        "D_mean": np.array([0.809, 0.863, 0.934, 0.974, 0.994]),
        "D_std":  np.array([0.060, 0.080, 0.033, 0.022, 0.011]),
        "R_mean": np.array([0.354, 0.402, 0.533, 0.563, 0.698]),
        "R_std":  np.array([0.090, 0.103, 0.064, 0.142, 0.150]),
    }),
    ("s=2.0e-2 (v.strong)", {
        "color": "#E91E63", "marker": "p",
        "folds":  np.array([30, 100, 300, 1000, 10000]),
        "D_mean": np.array([0.910, 0.939, 0.982, 0.984, 0.997]),
        "D_std":  np.array([0.045, 0.044, 0.022, 0.026, 0.007]),
        "R_mean": np.array([0.531, 0.614, 0.772, 0.706, 0.776]),
        "R_std":  np.array([0.108, 0.125, 0.136, 0.127, 0.177]),
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
ax_d.set_xlim(20, 15000)
ax_d.set_xticks([30, 100, 300, 1000, 10000])
ax_d.set_xticklabels(["30", "100", "300", "1K", "10K"])
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
ax_r.set_xlim(20, 15000)
ax_r.set_xticks([30, 100, 300, 1000, 10000])
ax_r.set_xticklabels(["30", "100", "300", "1K", "10K"])
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

plt.savefig("grid_reps_DR_final.pdf", dpi=300, bbox_inches="tight")
plt.savefig("grid_reps_DR_final.png", dpi=300, bbox_inches="tight")
print("Saved: grid_reps_DR_final.pdf / grid_reps_DR_final.png")
