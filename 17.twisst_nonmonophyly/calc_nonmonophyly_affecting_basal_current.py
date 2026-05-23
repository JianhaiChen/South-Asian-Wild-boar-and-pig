#!/usr/bin/env python3

import argparse
import gzip
import os
from collections import Counter, defaultdict


class Node:
    __slots__ = ("name", "neighbors")

    def __init__(self, name=None):
        self.name = name
        self.neighbors = []


def parse_newick(newick):
    s = newick.strip()
    if s.endswith(";"):
        s = s[:-1]

    nodes = []

    def new_node(name=None):
        node = Node(name)
        nodes.append(node)
        return node

    def add_edge(a, b):
        a.neighbors.append(b)
        b.neighbors.append(a)

    def skip_branch_length(i):
        if i < len(s) and s[i] == ":":
            i += 1
            while i < len(s) and s[i] not in ",()":
                i += 1
        return i

    def parse_label(i):
        start = i
        while i < len(s) and s[i] not in ":,()":
            i += 1
        return s[start:i], i

    def parse_subtree(i):
        if s[i] == "(":
            node = new_node()
            i += 1
            while True:
                child, i = parse_subtree(i)
                add_edge(node, child)
                if i >= len(s):
                    break
                if s[i] == ",":
                    i += 1
                    continue
                if s[i] == ")":
                    i += 1
                    break
                raise ValueError(f"unexpected character in Newick: {s[i]}")

            # Ignore optional internal labels.
            if i < len(s) and s[i] not in ":,()":
                _, i = parse_label(i)
            i = skip_branch_length(i)
            return node, i

        label, i = parse_label(i)
        if not label:
            raise ValueError("empty tip label")
        node = new_node(label)
        i = skip_branch_length(i)
        return node, i

    root, end = parse_subtree(0)
    if end != len(s):
        raise ValueError(f"trailing Newick text at position {end}")

    leaves = {node.name: node for node in nodes if node.name is not None}
    return root, nodes, leaves


def read_group_file(path):
    group_to_tips = defaultdict(set)
    tip_to_group = {}

    with open(path) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            tip, group = parts[0], parts[1]
            group_to_tips[group].add(tip)
            tip_to_group[tip] = group

    return dict(group_to_tips), tip_to_group


def read_trees(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line:
                yield line


def build_tree_context(nodes, leaves):
    root = nodes[0]
    parent = {root: None}
    order = [root]
    stack = [root]

    while stack:
        node = stack.pop()
        for nxt in node.neighbors:
            if nxt is parent.get(node):
                continue
            parent[nxt] = node
            order.append(nxt)
            stack.append(nxt)

    subtree_leaves = {}
    for node in reversed(order):
        node_leaves = set()
        if node.name is not None:
            node_leaves.add(node.name)
        for nxt in node.neighbors:
            if parent.get(nxt) is node:
                node_leaves |= subtree_leaves[nxt]
        subtree_leaves[node] = node_leaves

    all_leaves = set(leaves)

    def side_leaves(start, blocked):
        if parent.get(start) is blocked:
            return subtree_leaves[start]
        if parent.get(blocked) is start:
            return all_leaves - subtree_leaves[blocked]
        raise ValueError("nodes are not connected by a tree edge")

    return {
        "nodes": nodes,
        "leaves": leaves,
        "all_leaves": all_leaves,
        "side_leaves": side_leaves,
    }


def find_split_edge(ctx, target):
    all_leaves = ctx["all_leaves"]
    target = set(target) & all_leaves
    if len(target) < 2 or len(target) == len(all_leaves):
        return None

    for node in ctx["nodes"]:
        for nxt in node.neighbors:
            side = ctx["side_leaves"](node, nxt)
            if side == target:
                return node, nxt

    return None


def strict_monophyly(ctx, target):
    all_leaves = ctx["all_leaves"]
    target = set(target) & all_leaves
    if len(target) < 2 or len(target) == len(all_leaves):
        return "NA"
    return 1 if find_split_edge(ctx, target) is not None else 0


def group_composition(leaves, tip_to_group):
    counts = Counter()
    for tip in leaves:
        counts[tip_to_group.get(tip, "UNKNOWN")] += 1
    return counts


def comp_to_string(counts):
    if not counts:
        return "."
    return ",".join(f"{group}:{counts[group]}" for group in sorted(counts))


def classify_basal_branch(ctx, group_to_tips, tip_to_group, ingroup_groups, outgroup):
    all_leaf_names = ctx["all_leaves"]

    out_tips = set(group_to_tips.get(outgroup, set())) & all_leaf_names
    if len(out_tips) < 2:
        return {
            "call": "outgroup_not_testable",
            "basal_group": "NA",
            "left_comp": "NA",
            "right_comp": "NA",
            "n_children": "NA",
        }

    out_edge = find_split_edge(ctx, out_tips)
    if out_edge is None:
        return {
            "call": "outgroup_nonmonophyletic",
            "basal_group": "NA",
            "left_comp": "NA",
            "right_comp": "NA",
            "n_children": "NA",
        }

    out_side_node, ingroup_root = out_edge
    ingroup_tips = set()
    for group in ingroup_groups:
        ingroup_tips |= group_to_tips[group]
    ingroup_tips &= all_leaf_names

    ingroup_side = all_leaf_names - out_tips
    non_ingroup_inside = [tip for tip in ingroup_side if tip_to_group.get(tip) not in ingroup_groups]
    missing_ingroup = ingroup_tips - ingroup_side
    if non_ingroup_inside or missing_ingroup:
        return {
            "call": "ingroup_nonmonophyletic_relative_to_outgroup",
            "basal_group": "NA",
            "left_comp": "NA",
            "right_comp": "NA",
            "n_children": len([n for n in ingroup_root.neighbors if n is not out_side_node]),
        }

    child_info = []
    single_group_children = []
    ingroup_children = [n for n in ingroup_root.neighbors if n is not out_side_node]

    for child in ingroup_children:
        side = ctx["side_leaves"](child, ingroup_root)
        comp = group_composition(side, tip_to_group)
        ing_comp = {g: comp[g] for g in ingroup_groups if comp.get(g, 0) > 0}

        child_info.append((side, comp, ing_comp))

        if len(ing_comp) == 1:
            single_group_children.append(next(iter(ing_comp)))

    unique_single = sorted(set(single_group_children))

    if len(unique_single) == 1:
        basal_group = unique_single[0]
        call = f"{basal_group}_basal_branch"
    elif len(unique_single) == 0:
        basal_group = "NA"
        call = "mixed_ambiguous_basal_split"
    else:
        basal_group = "NA"
        call = "multiple_single_group_children_ambiguous"

    comps = [comp_to_string(info[1]) for info in child_info]
    while len(comps) < 2:
        comps.append(".")

    return {
        "call": call,
        "basal_group": basal_group,
        "left_comp": comps[0],
        "right_comp": comps[1],
        "n_children": len(ingroup_children),
    }


def safe_int(value):
    if value == "NA":
        return None
    return int(value)


def pct(numerator, denominator):
    if denominator == 0:
        return "NA"
    return f"{numerator / denominator:.6f}"


def summarize(label, records, focal_group, alternative_groups):
    n = len(records)
    if n == 0:
        return []

    def count(field, value=1):
        return sum(1 for rec in records if rec[field] == value)

    focal_compat = count(f"{focal_group}_compat")
    alt_counts = {group: count(f"{group}_compat") for group in alternative_groups}
    alternative_total = sum(alt_counts.values())
    ambiguous = count("ambiguous_focal")
    affects = count("affects_focal")

    nonmono_records = [rec for rec in records if rec["any_nonmono"] == 1]
    n_nonmono = len(nonmono_records)

    rows = []

    def add(metric, numerator, denominator):
        rows.append([label, metric, str(denominator), str(numerator), pct(numerator, denominator)])

    add(f"{focal_group}_basal_branch_compatible", focal_compat, n)
    for group in alternative_groups:
        add(f"{group}_basal_branch_compatible", alt_counts[group], n)
    add(f"alternative_to_{focal_group}_basal", alternative_total, n)
    add(f"ambiguous_for_{focal_group}_basal", ambiguous, n)
    add(f"affects_or_weakens_{focal_group}_basal", affects, n)

    add("any_strict_ingroup_nonmonophyly", n_nonmono, n)
    add(
        f"nonmonophyletic_but_still_{focal_group}_basal",
        sum(1 for rec in nonmono_records if rec[f"{focal_group}_compat"] == 1),
        n_nonmono,
    )
    add(
        f"nonmonophyletic_and_alternative_to_{focal_group}_basal",
        sum(1 for rec in nonmono_records if rec["alternative_focal"] == 1),
        n_nonmono,
    )
    add(
        f"nonmonophyletic_and_ambiguous_for_{focal_group}_basal",
        sum(1 for rec in nonmono_records if rec["ambiguous_focal"] == 1),
        n_nonmono,
    )
    add(
        f"nonmonophyletic_and_affects_or_weakens_{focal_group}_basal",
        sum(1 for rec in nonmono_records if rec["affects_focal"] == 1),
        n_nonmono,
    )

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hap-group", default="id_hap.gr")
    parser.add_argument("--tree-dir", default=".")
    parser.add_argument("--win", default="5000")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--ingroup-groups", default="saw,msea,isea")
    parser.add_argument("--outgroup", default="sus")
    parser.add_argument("--focal-group", default="saw")
    parser.add_argument("--chroms", default=",".join([str(i) for i in range(1, 19)] + ["X"]))
    parser.add_argument(
        "--tree-template",
        default="{chrom}.phyml_bionj.w{win}.trees.gz",
        help="Template relative to --tree-dir; fields: chrom, win.",
    )
    args = parser.parse_args()

    group_to_tips, tip_to_group = read_group_file(args.hap_group)
    ingroup_groups = [x.strip() for x in args.ingroup_groups.split(",") if x.strip()]
    chroms = [x.strip().replace("Chr", "") for x in args.chroms.split(",") if x.strip()]
    focal_group = args.focal_group
    alternative_groups = [g for g in ingroup_groups if g != focal_group]

    for group in ingroup_groups + [args.outgroup]:
        if group not in group_to_tips:
            raise SystemExit(f"[ERROR] group not found in {args.hap_group}: {group}")

    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)

    by_window_path = args.out_prefix + ".by_window.tsv"
    summary_path = args.out_prefix + ".summary.tsv"
    by_chr_path = args.out_prefix + ".by_chr_summary.tsv"
    strict_path = args.out_prefix + ".strict_monophyly.tsv"

    records = []
    strict_counts = defaultdict(lambda: defaultdict(int))

    with open(by_window_path, "w") as out:
        strict_header = "\t".join(f"{group}_strict_mono" for group in ingroup_groups)
        compat_header = "\t".join(f"{group}_basal_branch_compatible" for group in ingroup_groups)
        out.write(
            "chr\twindow\tregion_type\t"
            f"{strict_header}\tany_ingroup_strict_nonmono\t"
            "basal_branch_call\tbasal_group\t"
            f"{compat_header}\t"
            f"alternative_to_{focal_group}_basal\tambiguous_for_{focal_group}_basal\t"
            f"affects_or_weakens_{focal_group}_basal\t"
            "left_child_composition\tright_child_composition\tn_ingroup_root_children\n"
        )

        for chrom in chroms:
            tree_file = os.path.join(
                args.tree_dir,
                args.tree_template.format(chrom=chrom, win=args.win),
            )
            if not os.path.exists(tree_file):
                print(f"[WARN] missing tree file: {tree_file}")
                continue

            region_type = "X" if chrom == "X" else "Autosome"
            processed = 0

            for window_idx, nwk in enumerate(read_trees(tree_file), start=1):
                processed += 1
                _, nodes, leaves = parse_newick(nwk)
                ctx = build_tree_context(nodes, leaves)

                strict = {
                    group: strict_monophyly(ctx, group_to_tips[group])
                    for group in ingroup_groups
                }
                strict_vals = [safe_int(strict[group]) for group in ingroup_groups if strict[group] != "NA"]
                any_nonmono = "NA" if not strict_vals else (1 if any(v == 0 for v in strict_vals) else 0)

                for group in ingroup_groups:
                    strict_counts[(region_type, group)]["den"] += 1
                    if strict[group] == 1:
                        strict_counts[(region_type, group)]["mono"] += 1
                    elif strict[group] == 0:
                        strict_counts[(region_type, group)]["nonmono"] += 1

                basal = classify_basal_branch(
                    ctx,
                    group_to_tips,
                    tip_to_group,
                    ingroup_groups,
                    args.outgroup,
                )
                call = basal["call"]

                compat = {group: 1 if call == f"{group}_basal_branch" else 0 for group in ingroup_groups}
                alternative = 1 if any(compat[group] == 1 for group in alternative_groups) else 0
                ambiguous = 1 if not any(compat[group] == 1 for group in ingroup_groups) else 0
                affects = 0 if compat[focal_group] == 1 else 1

                rec = {
                    "chr": chrom,
                    "region_type": region_type,
                    "any_nonmono": any_nonmono,
                    "alternative_focal": alternative,
                    "ambiguous_focal": ambiguous,
                    "affects_focal": affects,
                }
                for group in ingroup_groups:
                    rec[f"{group}_compat"] = compat[group]
                records.append(rec)

                out.write(
                    f"{chrom}\t{window_idx}\t{region_type}\t"
                    + "\t".join(str(strict[group]) for group in ingroup_groups)
                    + f"\t{any_nonmono}\t"
                    f"{call}\t{basal['basal_group']}\t"
                    + "\t".join(str(compat[group]) for group in ingroup_groups)
                    + f"\t{alternative}\t{ambiguous}\t{affects}\t"
                    f"{basal['left_comp']}\t{basal['right_comp']}\t{basal['n_children']}\n"
                )

            print(f"[INFO] {chrom}: {processed} trees processed")

    with open(summary_path, "w") as out:
        out.write("set\tmetric\tn_denominator\tn_numerator\tfrequency\n")
        for label, subset in [
            ("All", records),
            ("Autosome", [r for r in records if r["region_type"] == "Autosome"]),
            ("X", [r for r in records if r["region_type"] == "X"]),
        ]:
            for row in summarize(label, subset, focal_group, alternative_groups):
                out.write("\t".join(row) + "\n")

    with open(by_chr_path, "w") as out:
        out.write("chr\tmetric\tn_denominator\tn_numerator\tfrequency\n")
        for chrom in chroms:
            subset = [r for r in records if r["chr"] == chrom]
            for row in summarize(chrom, subset, focal_group, alternative_groups):
                out.write("\t".join([chrom] + row[1:]) + "\n")

    with open(strict_path, "w") as out:
        out.write("set\tgroup\tn_denominator\tn_strict_mono\tn_strict_nonmono\tstrict_mono_frequency\tstrict_nonmono_frequency\n")
        for label in ["Autosome", "X"]:
            for group in ingroup_groups:
                vals = strict_counts[(label, group)]
                out.write(
                    f"{label}\t{group}\t{vals['den']}\t{vals['mono']}\t{vals['nonmono']}\t"
                    f"{pct(vals['mono'], vals['den'])}\t{pct(vals['nonmono'], vals['den'])}\n"
                )
        for group in ingroup_groups:
            den = strict_counts[("Autosome", group)]["den"] + strict_counts[("X", group)]["den"]
            mono = strict_counts[("Autosome", group)]["mono"] + strict_counts[("X", group)]["mono"]
            nonmono = strict_counts[("Autosome", group)]["nonmono"] + strict_counts[("X", group)]["nonmono"]
            out.write(
                f"All\t{group}\t{den}\t{mono}\t{nonmono}\t"
                f"{pct(mono, den)}\t{pct(nonmono, den)}\n"
            )

    print(f"[DONE] {by_window_path}")
    print(f"[DONE] {summary_path}")
    print(f"[DONE] {by_chr_path}")
    print(f"[DONE] {strict_path}")


if __name__ == "__main__":
    main()
