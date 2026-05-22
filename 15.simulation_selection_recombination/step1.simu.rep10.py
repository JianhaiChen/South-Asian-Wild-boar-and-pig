#!/usr/bin/env python3
"""Add 50x/70x with N replicates to reduce stochastic noise."""

import msprime
import numpy as np
import time

MU, REC_BG = 2.7e-9, 1e-8
NE_LW, NE_SUS, NE_GHOST = 500_000, 20_000, 15_000
NE_ROOT, NE_MINI = 30_000, 10_000
T_DOM, T_INTRO, T_DEEP = 3_333, 1_000_000, 2_500_000

L_FLANK = 250_000
L_SV    = 10_000_000
L_TOTAL = 2 * L_FLANK + L_SV
sv_s, sv_e = L_FLANK, L_FLANK + L_SV

N_MINI, N_LARGE = 5, 5
N_PAIRS = 6000
N_REPS  = 10         # <-- KEY: 10 replicates per cell
SCALE   = 4.4

def build_dem(ne_large):
    d = msprime.Demography()
    d.add_population(name="mini",       initial_size=NE_MINI)
    d.add_population(name="large",      initial_size=ne_large)
    d.add_population(name="sus",        initial_size=NE_SUS)
    d.add_population(name="large_wild", initial_size=NE_LW)
    d.add_population(name="ghost",      initial_size=NE_GHOST)
    d.add_population(name="root",       initial_size=NE_ROOT)
    d.set_migration_rate("mini", "large", 5e-5)
    d.set_migration_rate("large", "mini", 5e-5)
    d.add_mass_migration(time=T_DOM, source="mini",
                         dest="sus", proportion=1.0)
    d.add_mass_migration(time=T_DOM+1, source="large",
                         dest="large_wild", proportion=1.0)
    d.add_migration_rate_change(time=T_DOM+2, rate=0)
    d.add_migration_rate_change(time=T_DOM+3, source="sus",
                                dest="large_wild", rate=1e-7)
    d.add_migration_rate_change(time=T_DOM+3, source="large_wild",
                                dest="sus", rate=1e-7)
    d.add_migration_rate_change(time=T_INTRO-1, source="sus",
                                dest="large_wild", rate=0)
    d.add_migration_rate_change(time=T_INTRO-1, source="large_wild",
                                dest="sus", rate=0)
    d.add_mass_migration(time=T_INTRO, source="large_wild",
                         dest="ghost", proportion=1.0)
    d.add_migration_rate_change(time=T_INTRO+1, source="sus",
                                dest="ghost", rate=1e-7)
    d.add_migration_rate_change(time=T_INTRO+1, source="ghost",
                                dest="sus", rate=1e-7)
    d.add_migration_rate_change(time=T_DEEP-1, rate=0)
    d.add_population_split(time=T_DEEP, derived=["sus", "ghost"],
                           ancestral="root")
    d.sort_events()
    return d

# 修改 1：在参数列表中增加了 `seed` 变量
def compute_ld(mts, haps, start, end, seed, max_sites=1200, max_pairs=N_PAIRS):
    gs = []
    for v in mts.variants():
        if v.site.position < start: continue
        if v.site.position >= end: break
        g = v.genotypes[haps]
        f = np.mean(g)
        if 0.05 < f < 0.95:
            gs.append(g.copy())
            
    if len(gs) < 20: return np.nan, np.nan
    gs = np.array(gs)
    
    # 修改 2：使用传入的动态 seed 初始化 NumPy RNG，替代硬编码的 42
    rng = np.random.default_rng(seed)
    
    if len(gs) > max_sites:
        gs = gs[rng.choice(len(gs), max_sites, replace=False)]
        
    dp, rs = [], []
    for _ in range(max_pairs):
        i, j = rng.choice(len(gs), 2, replace=False)
        a, b = gs[i].astype(float), gs[j].astype(float)
        pA, pB = np.mean(a), np.mean(b)
        D = np.mean(a * b) - pA * pB
        pa, pb = 1 - pA, 1 - pB
        dn = pA * pa * pB * pb
        if dn == 0: continue
        Dm = min(pA*pb, pa*pB) if D > 0 else min(pA*pB, pa*pb)
        if Dm == 0: continue
        dp.append(min(abs(D) / Dm, 1.0))
        rs.append(min(D**2 / dn, 1.0))
        
    if not dp: return np.nan, np.nan
    return np.mean(dp), np.mean(rs)

# ---- Grid: your s value (Ne/167x) at all rec folds ----
ne_large = int(NE_LW / 167)  # = 2994 ≈ 3000

# Test all fold values including new 50x, 70x
rec_folds = [10, 30, 50, 70, 100, 300, 1000, 3000, 10000]

print("=" * 72)
print(f" Replicated grid: s=5.8e-03 (Ne={ne_large:,}) x {N_REPS} reps")
print(f" Rec folds: {rec_folds}")
print("=" * 72)

t_start = time.time()

for fold in rec_folds:
    rec_sv = REC_BG / fold
    rho = 4 * NE_LW * rec_sv * L_SV

    rm = msprime.RateMap(
        position=[0, sv_s, sv_e, L_TOTAL],
        rate=[REC_BG, rec_sv, REC_BG])

    d_list, r_list = [], []

    for rep in range(N_REPS):
        seed = 1000 * fold + rep + 1
        ts = msprime.sim_ancestry(
            samples={"mini": N_MINI, "large": N_LARGE},
            demography=build_dem(ne_large),
            recombination_rate=rm,
            sequence_length=L_TOTAL,
            random_seed=seed)
        mts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)
        haps = list(ts.samples())
        
        # 修改 3：调用时传入 seed
        d, r = compute_ld(mts, haps, sv_s, sv_e, seed)
        
        if not np.isnan(d):
            d_list.append(d)
            r_list.append(r)

    d_mean = np.mean(d_list)
    d_std  = np.std(d_list)
    r_mean = np.mean(r_list)
    r_std  = np.std(r_list)

    match = ""
    if d_mean > 0.9 and 0.4 < r_mean < 0.8:
        match = " *MATCH"

    print(f"  {fold:>5}x  rho={rho:>7,.0f}  "
          f"D'={d_mean:.3f}+/-{d_std:.3f}  "
          f"R2={r_mean:.3f}+/-{r_std:.3f}  "
          f"(n={len(d_list)}){match}")

elapsed = time.time() - t_start
print(f"\nDone: {elapsed:.0f}s = {elapsed/60:.1f} min")

# ---- Also run a few other s values at key folds ----
print(f"\n{'='*72}")
print(f" Cross-validation: 3 s values x 5 folds x {N_REPS} reps")
print(f"{'='*72}")

test_configs = [
    (25,  "s=6.5e-04"),   # weak selection
    (100, "s=3.2e-03"),   # moderate
    (500, "s=2.0e-02"),   # very strong
]
test_folds = [30, 100, 300, 1000, 10000]

for ne_fold, s_label in test_configs:
    ne_l = int(NE_LW / ne_fold)
    print(f"\n  --- {s_label} (Ne={ne_l:,}) ---")

    for fold in test_folds:
        rec_sv = REC_BG / fold
        rho = 4 * NE_LW * rec_sv * L_SV

        rm = msprime.RateMap(
            position=[0, sv_s, sv_e, L_TOTAL],
            rate=[REC_BG, rec_sv, REC_BG])

        d_list, r_list = [], []
        for rep in range(N_REPS):
            seed = ne_fold * 10000 + fold * 100 + rep + 1
            ts = msprime.sim_ancestry(
                samples={"mini": N_MINI, "large": N_LARGE},
                demography=build_dem(ne_l),
                recombination_rate=rm,
                sequence_length=L_TOTAL,
                random_seed=seed)
            mts = msprime.sim_mutations(
                ts, rate=MU, random_seed=seed)
            haps = list(ts.samples())
            
            # 修改 4：调用时传入 seed
            d, r = compute_ld(mts, haps, sv_s, sv_e, seed)
            
            if not np.isnan(d):
                d_list.append(d)
                r_list.append(r)

        d_mean = np.mean(d_list)
        d_std  = np.std(d_list)
        r_mean = np.mean(r_list)
        r_std  = np.std(r_list)

        match = ""
        if d_mean > 0.9 and 0.4 < r_mean < 0.8:
            match = " *MATCH"

        print(f"    {fold:>5}x  rho={rho:>7,.0f}  "
              f"D'={d_mean:.3f}+/-{d_std:.3f}  "
              f"R2={r_mean:.3f}+/-{r_std:.3f}{match}")

print(f"""
{'='*72}
 INTERPRETATION GUIDE
{'='*72}

 If std < 0.05:  result is stable (trustworthy)
 If std > 0.10:  result is noisy (single runs unreliable!)

 Previous single-run results that might be noise:
   70x:  D'=0.694 R2=0.222  (suspicious: worse than 50x!)
   100x: D'=0.954 R2=0.513  (the R2 jump might be smaller)
   1000x: R2=0.453 < 300x R2=0.494 (U-shape might be noise)

 After averaging {N_REPS} reps, the TRUE pattern should be:
   monotonically increasing D' with fold
   monotonically increasing R2 with fold (or a REAL plateau)

 sim.py update:
   Use the fold where mean D' > 0.95 AND mean R2 ~ 0.55-0.65
   Scale by x{SCALE:.1f} for full 44 Mb simulation
{'='*72}
""")
