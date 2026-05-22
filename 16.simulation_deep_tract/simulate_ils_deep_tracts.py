#!/usr/bin/env python3
import argparse
import msprime
import numpy as np
import pandas as pd

def parse_list(x, dtype=float):
    return [dtype(i) for i in str(x).split(",") if i.strip()]

def longest_deep_tract(ts, deep_gen):
    """
    Pairwise TMRCA between one haploid sample from L and one from M.
    Return longest continuous interval with TMRCA >= deep_gen.
    """
    pop_names = [p.metadata.get("name", "") for p in ts.populations()]
    pop_id = {name: i for i, name in enumerate(pop_names)}

    L_nodes = list(ts.samples(population=pop_id["L"]))
    M_nodes = list(ts.samples(population=pop_id["M"]))
    if len(L_nodes) == 0 or len(M_nodes) == 0:
        raise RuntimeError("Cannot find samples from populations L and M.")

    u = L_nodes[0]
    v = M_nodes[0]

    longest = 0.0
    current = 0.0

    for tree in ts.trees():
        left, right = tree.interval
        span = right - left
        tmrca = tree.tmrca(u, v)

        if tmrca >= deep_gen:
            current += span
        else:
            if current > longest:
                longest = current
            current = 0.0

    if current > longest:
        longest = current

    return longest

def run_one(length_bp, recomb_rate, split_gen, ne, ne_anc, deep_gen, seed):
    dem = msprime.Demography()
    dem.add_population(name="L", initial_size=ne)
    dem.add_population(name="M", initial_size=ne)
    dem.add_population(name="ANC", initial_size=ne_anc)
    dem.add_population_split(time=split_gen, derived=["L", "M"], ancestral="ANC")
    dem.sort_events()

    ts = msprime.sim_ancestry(
        samples={"L": 1, "M": 1},
        ploidy=1,
        sequence_length=length_bp,
        recombination_rate=recomb_rate,
        demography=dem,
        random_seed=seed,
    )
    return longest_deep_tract(ts, deep_gen)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="ils_tract_sim")
    ap.add_argument("--n-reps", type=int, default=1000)
    ap.add_argument("--length-mb", type=float, default=43.0)
    ap.add_argument("--observed-mb", type=float, default=40.0)
    ap.add_argument("--generation-time", type=float, default=3.0)
    ap.add_argument("--split-kya", type=float, default=695.0)
    ap.add_argument("--deep-mya", type=float, default=4.5)
    ap.add_argument("--ne", type=float, default=20000)
    ap.add_argument("--ne-anc-list", default="50000,100000,200000")
    ap.add_argument(
        "--recomb-cm-mb-list",
        default="0.02,0.05,0.10",
        help="1 cM/Mb = 1e-8 recombination per bp per generation",
    )
    args = ap.parse_args()

    length_bp = int(args.length_mb * 1_000_000)
    observed_bp = args.observed_mb * 1_000_000
    split_gen = args.split_kya * 1000 / args.generation_time
    deep_gen = args.deep_mya * 1_000_000 / args.generation_time

    ne_anc_list = parse_list(args.ne_anc_list, float)
    recomb_cm_mb_list = parse_list(args.recomb_cm_mb_list, float)

    rows = []
    seed = 1000

    for ne_anc in ne_anc_list:
        for cm_mb in recomb_cm_mb_list:
            recomb_rate = cm_mb * 1e-8
            vals = []
            for rep in range(args.n_reps):
                seed += 1
                longest_bp = run_one(
                    length_bp=length_bp,
                    recomb_rate=recomb_rate,
                    split_gen=split_gen,
                    ne=args.ne,
                    ne_anc=ne_anc,
                    deep_gen=deep_gen,
                    seed=seed,
                )
                vals.append(longest_bp / 1_000_000)

            vals = np.array(vals)
            p_emp = (np.sum(vals >= args.observed_mb) + 1) / (len(vals) + 1)

            rows.append({
                "Ne_present": args.ne,
                "Ne_ancestral": ne_anc,
                "recomb_cM_per_Mb": cm_mb,
                "n_reps": args.n_reps,
                "deep_TMRCA_cutoff_Mya": args.deep_mya,
                "observed_tract_Mb": args.observed_mb,
                "median_longest_deep_tract_Mb": np.median(vals),
                "q95_longest_deep_tract_Mb": np.quantile(vals, 0.95),
                "q99_longest_deep_tract_Mb": np.quantile(vals, 0.99),
                "max_longest_deep_tract_Mb": np.max(vals),
                "empirical_P_longest_ge_observed": p_emp,
            })

            pd.DataFrame({
                "longest_deep_tract_Mb": vals,
                "Ne_ancestral": ne_anc,
                "recomb_cM_per_Mb": cm_mb,
            }).to_csv(
                f"{args.out}.NeAnc{int(ne_anc)}.rec{cm_mb}.raw.csv",
                index=False
            )

    df = pd.DataFrame(rows)
    df.to_csv(f"{args.out}.summary.csv", index=False)
    print(df.to_string(index=False))
    print(f"\n[DONE] Summary: {args.out}.summary.csv")

if __name__ == "__main__":
    main()