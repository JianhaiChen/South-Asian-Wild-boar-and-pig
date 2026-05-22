#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("ils_X43_obs40_rec0p01_NeAnc500k_deep_sweep.summary.csv")
df = df.sort_values("deep_TMRCA_cutoff_Mya")

# Convert Mb to kb for zoomed visualization
for col in [
    "median_longest_deep_tract_Mb",
    "q95_longest_deep_tract_Mb",
    "q99_longest_deep_tract_Mb",
    "max_longest_deep_tract_Mb",
]:
    df[col.replace("_Mb", "_kb")] = df[col] * 1000

# Figure 1: full-scale comparison with observed 40 Mb
plt.figure(figsize=(6.2, 4.2))
plt.plot(df["deep_TMRCA_cutoff_Mya"], df["q95_longest_deep_tract_Mb"], marker="o", label="95th percentile")
plt.plot(df["deep_TMRCA_cutoff_Mya"], df["q99_longest_deep_tract_Mb"], marker="o", label="99th percentile")
plt.plot(df["deep_TMRCA_cutoff_Mya"], df["max_longest_deep_tract_Mb"], marker="o", label="Maximum")
plt.axhline(40, linestyle="--", linewidth=1.5, label="Observed tract >40 Mb")
plt.xlabel("Deep-TMRCA cutoff (Mya)")
plt.ylabel("Longest deep-TMRCA tract under ILS-only model (Mb)")
plt.title("Full-scale comparison with observed tract")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("ils_deep_cutoff_sweep_fullscale.pdf")
plt.savefig("ils_deep_cutoff_sweep_fullscale.png", dpi=300)
plt.close()

# Figure 2: zoomed simulated tract lengths in kb
plt.figure(figsize=(6.2, 4.2))
plt.plot(df["deep_TMRCA_cutoff_Mya"], df["median_longest_deep_tract_kb"], marker="o", label="Median")
plt.plot(df["deep_TMRCA_cutoff_Mya"], df["q95_longest_deep_tract_kb"], marker="o", label="95th percentile")
plt.plot(df["deep_TMRCA_cutoff_Mya"], df["q99_longest_deep_tract_kb"], marker="o", label="99th percentile")
plt.plot(df["deep_TMRCA_cutoff_Mya"], df["max_longest_deep_tract_kb"], marker="o", label="Maximum")
plt.xlabel("Deep-TMRCA cutoff (Mya)")
plt.ylabel("Longest deep-TMRCA tract under ILS-only model (kb)")
plt.title("Zoomed ILS-only tract lengths")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("ils_deep_cutoff_sweep_zoom_kb.pdf")
plt.savefig("ils_deep_cutoff_sweep_zoom_kb.png", dpi=300)
plt.close()

print("[DONE] ils_deep_cutoff_sweep_fullscale.pdf")
print("[DONE] ils_deep_cutoff_sweep_zoom_kb.pdf")