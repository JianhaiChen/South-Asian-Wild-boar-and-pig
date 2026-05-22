#!/usr/bin/env python3
"""
sim_two_scenarios.py -- Two extreme scenarios → VCF → LDBlockShow

Scenario A: Strong selection + weak recombination suppression
  Ne_large=1,000 (s≈0.02), REC_SV=7.57e-11 (full 132x, matching grid 30x)

Scenario B: Weak selection + ultra-strong recombination suppression
  Ne_large=20,000 (s≈6.5e-4), REC_SV=2.27e-13 (full 44000x, matching grid 10000x)

Both should produce heatmaps similar to empirical data in LDBlockShow.
"""

import msprime
import numpy as np
import time
import os

# ==================================================================
# 1. SHARED PARAMETERS
# ==================================================================
MU     = 2.7e-9
REC_BG = 1e-8

L        = 64_000_000
SV_START = 10_000_000
SV_END   = 54_000_000
L_SV     = SV_END - SV_START  # 44 Mb

NE_LW    = 500_000
NE_SUS   = 20_000
NE_GHOST = 15_000
NE_ROOT  = 30_000
NE_MINI  = 10_000

T_DOM   = 3_333
T_INTRO = 1_000_000
T_DEEP  = 2_500_000

N_MINI  = 10
N_LARGE = 10
N_VCF   = 5000
N_PAIRS = 10000

SCALE = L_SV / 10_000_000  # 4.4

# ==================================================================
# 2. TWO SCENARIOS (UPDATED BASED ON GRID SEARCH)
# ==================================================================
scenarios = {
    "A_strong_sel": {
        "label": "Strong selection + weak suppression",
        "ne_large": 1_000,
        "s_approx": 0.02,
        "ns2": 20000,
        "grid_fold": 30,
        "full_fold": int(30 * SCALE), # 132
        "rec_sv": REC_BG / (30 * SCALE), # ~ 7.57e-11
        "grid_D": "0.910+/-0.045",
        "grid_R": "0.531+/-0.108",
        "seed": 42,
    },
    "B_weak_sel": {
        "label": "Weak selection + ultra-strong suppression",
        "ne_large": 20_000,
        "s_approx": 6.5e-4,
        "ns2": 650,
        "grid_fold": 10000,
        "full_fold": int(10000 * SCALE), # 44000
        "rec_sv": REC_BG / (10000 * SCALE), # ~ 2.27e-13
        "grid_D": "0.970+/-0.019",
        "grid_R": "0.451+/-0.134",
        "seed": 123,
    },
}

# ==================================================================
# 3. DEMOGRAPHIC MODEL
# ==================================================================
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

# ==================================================================
# 4. LD COMPUTATION (RNG Fixed)
# ==================================================================
def quick_ld(mts, haps, start, end, seed, max_sites=1500, max_pairs=N_PAIRS):
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
    
    # 使用传入的独立 seed
    rng = np.random.default_rng(seed)
    
    if len(gs) > max_sites:
        gs = gs[rng.choice(len(gs), max_sites, replace=False)]
        
    dp, rs = [], []
    for _ in range(max_pairs):
        i, j = rng.choice(len(gs), 2, replace=False)
        a, b = gs[i].astype(float), gs[j].astype(float)
        pA, pB = np.mean(a), np.mean(b)
        D = np.mean(a*b) - pA*pB
        pa, pb = 1-pA, 1-pB
        dn = pA*pa*pB*pb
        if dn == 0: continue
        Dm = min(pA*pb, pa*pB) if D > 0 else min(pA*pB, pa*pb)
        if Dm == 0: continue
        dp.append(min(abs(D)/Dm, 1.0))
        rs.append(min(D**2/dn, 1.0))
        
    if not dp: return np.nan, np.nan
    return np.mean(dp), np.mean(rs)

# ==================================================================
# 5. RUN BOTH SCENARIOS
# ==================================================================
results = {}

for name, cfg in scenarios.items():
    rho = 4 * NE_LW * cfg["rec_sv"] * L_SV

    print(f"\n{'='*70}")
    print(f" SCENARIO {name}: {cfg['label']}")
    print(f"{'='*70}")
    print(f"   Ne_large = {cfg['ne_large']:,}  "
          f"(s ~ {cfg['s_approx']}, 2Ns ~ {cfg['ns2']})")
    print(f"   REC_SV   = {cfg['rec_sv']:.2e}  "
          f"(full {cfg['full_fold']}x, grid {cfg['grid_fold']}x)")
    print(f"   rho_SV   = {rho:.0f}")
    print(f"   Grid:      D' = {cfg['grid_D']}, R2 = {cfg['grid_R']}")
    print()

    rate_map = msprime.RateMap(
        position=[0, SV_START, SV_END, L],
        rate=[REC_BG, cfg["rec_sv"], REC_BG])

    # ---- Ancestry ----
    print(f"   [1/3] Ancestry ...", flush=True)
    t0 = time.time()
    ts = msprime.sim_ancestry(
        samples={"mini": N_MINI, "large": N_LARGE},
        demography=build_dem(cfg["ne_large"]),
        recombination_rate=rate_map,
        sequence_length=L,
        random_seed=cfg["seed"])
    t1 = time.time()
    print(f"         {ts.num_trees:,} trees ({t1-t0:.0f}s = "
          f"{(t1-t0)/60:.1f} min)")

    # ---- Mutations ----
    print(f"   [2/3] Mutations ...", flush=True)
    mts = msprime.sim_mutations(ts, rate=MU, random_seed=cfg["seed"])
    t2 = time.time()
    print(f"         {mts.num_sites:,} sites ({t2-t1:.0f}s)")

    # ---- LD check ----
    print(f"   [3/3] LD verification ...", flush=True)
    all_haps = list(ts.samples())
    
    # 调用 LD 函数时传入当前的 seed 保证随机独立性
    d_sv, r_sv = quick_ld(mts, all_haps, SV_START, SV_END, cfg["seed"])
    d_fl, r_fl = quick_ld(mts, all_haps, 0, SV_START, cfg["seed"])
    t3 = time.time()

    results[name] = {"D": d_sv, "R": r_sv, "D_fl": d_fl, "R_fl": r_fl}

    print(f"""
         Combined SV:    D' = {d_sv:.3f}  R2 = {r_sv:.3f}
         Flanking:       D' = {d_fl:.3f}  R2 = {r_fl:.3f}
         Target:         D' = 0.963  R2 = 0.595
    """)

    # ---- Write VCF ----
    vcf_file = f"sim_{name}.vcf"
    print(f"   Writing {vcf_file} ...", flush=True)

    all_ids = []
    for v in mts.variants():
        f = np.mean(v.genotypes)
        if 0.05 < f < 0.95:
            all_ids.append(v.site.id)
    print(f"         {len(all_ids):,} informative variants")

    if len(all_ids) > N_VCF:
        # 使用当前场景的 seed 随机抽取指定数量的 VCF 位点，比均匀切片更真实
        rng_vcf = np.random.default_rng(cfg["seed"])
        keep = set(rng_vcf.choice(all_ids, N_VCF, replace=False))
    else:
        keep = set(all_ids)

    remove = [i for i in range(mts.num_sites) if i not in keep]
    mts_sub = mts.delete_sites(remove)

    ind_names = ([f"miniP{i}" for i in range(1, N_MINI+1)] +
                 [f"largeP{i}" for i in range(1, N_LARGE+1)])

    with open(vcf_file, "w") as f:
        mts_sub.write_vcf(f, contig_id="1", individual_names=ind_names)

    t4 = time.time()
    print(f"         Saved: {vcf_file} ({mts_sub.num_sites:,} variants)")
    print(f"         Total: {t4-t0:.0f}s = {(t4-t0)/60:.1f} min")

# ==================================================================
# 6. COMPARISON SUMMARY
# ==================================================================
print(f"""
{'='*70}
 COMPARISON
{'='*70}

 Scenario A (strong s + weak suppression):
   Ne_large = 1,000 | REC_SV = {scenarios['A_strong_sel']['rec_sv']:.2e} ({scenarios['A_strong_sel']['full_fold']}x)
   D' = {results['A_strong_sel']['D']:.3f}  R2 = {results['A_strong_sel']['R']:.3f}
   VCF: sim_A_strong_sel.vcf

 Scenario B (weak s + ultra-strong suppression):
   Ne_large = 20,000 | REC_SV = {scenarios['B_weak_sel']['rec_sv']:.2e} ({scenarios['B_weak_sel']['full_fold']}x)
   D' = {results['B_weak_sel']['D']:.3f}  R2 = {results['B_weak_sel']['R']:.3f}
   VCF: sim_B_weak_sel.vcf

 Target (empirical):
   D' = 0.963  R2 = 0.595

 KEY QUESTION:
   Do the LDBlockShow heatmaps look similar?
   -> If YES: LD alone cannot distinguish the two mechanisms
   -> Independent evidence (FST, SV analysis) is needed

{'='*70}
 Next: bash run_ldblockshow_compare.sh
{'='*70}
""")
