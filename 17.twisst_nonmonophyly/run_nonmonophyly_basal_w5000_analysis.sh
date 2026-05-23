#!/bin/bash
set -euo pipefail

# Reproduce the non-monophyly and basal-branch compatibility analysis used for Table S11.
# Input files expected in the current directory:
#   id_hap.gr
#   {1..18,X}.phyml_bionj.w5000.trees.gz
#   {1..18,X}.phyml_bionj.w5000.data.tsv

PREFIX="nonmonophyly_affecting_basal_w5000"

python3 calc_nonmonophyly_affecting_basal_current.py \
  --hap-group id_hap.gr \
  --tree-dir . \
  --win 5000 \
  --out-prefix "$PREFIX" \
  --ingroup-groups saw,msea,isea \
  --outgroup sus \
  --focal-group saw

# Build a coordinate lookup from phyml_sliding_windows data files.
awk 'BEGIN {
    FS=OFS="\t"
    print "chr","window","scaffold","start","end","mid","sites","lnL"
}' > "$PREFIX.window_coordinates.tsv"

for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X
do
  awk -v chr="$c" 'BEGIN {FS=OFS="\t"} NR>1 {
      print chr, NR-1, $1, $2, $3, $4, $5, $6
  }' "$c.phyml_bionj.w5000.data.tsv" >> "$PREFIX.window_coordinates.tsv"
done

# Add genomic coordinates to the full window-level classification table.
awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        if (NR>1) coord[$1 FS $2]=$3 FS $4 FS $5 FS $6 FS $7 FS $8
        next
    }
    FNR==1 {
        print $1,$2,"scaffold","start","end","mid","sites","lnL", \
              $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18
        next
    }
    {
        key=$1 FS $2
        print $1,$2,coord[key],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18
    }' "$PREFIX.window_coordinates.tsv" "$PREFIX.by_window.tsv" \
    > "$PREFIX.by_window.with_coordinates.tsv"

# Table S11: only windows in which at least one ingroup clade is not strictly monophyletic.
# In the coordinate-added table, any_ingroup_strict_nonmono is column 13.
awk 'BEGIN {FS=OFS="\t"} NR==1 || $13==1 {
    print
}' "$PREFIX.by_window.with_coordinates.tsv" \
  > "$PREFIX.any_ingroup_nonmono_windows.with_coordinates.tsv"

# Summary statistics used in the response sentence:
# 63.91% (13,446/21,040) autosomal and 73.49% (352/479) X-linked still South-Asia basal;
# 28.24% and 11.90% support a non-South-Asian basal branch;
# 7.86% and 14.61% are ambiguous.
awk 'BEGIN {
    FS=OFS="\t"
    print "Region", \
          "At_least_one_ingroup_strict_nonmono_windows", \
          "South_Asia_only_basal_windows", "South_Asia_only_basal_pct", \
          "Non_South_Asian_basal_windows", "Non_South_Asian_basal_pct", \
          "Ambiguous_windows", "Ambiguous_pct", \
          "Affects_or_weakens_South_Asia_basal_windows", "Affects_or_weakens_pct"
}
NR>1 && $7==1 {
    region=$3
    den[region]++
    saw_basal[region]+=$10
    non_sa_basal[region]+=$13
    ambiguous[region]+=$14
    affects[region]+=$15
}
END {
    for (i=1; i<=2; i++) {
        region=(i==1 ? "Autosome" : "X")
        printf "%s\t%d\t%d\t%.2f%%\t%d\t%.2f%%\t%d\t%.2f%%\t%d\t%.2f%%\n", \
               region, den[region], \
               saw_basal[region], 100*saw_basal[region]/den[region], \
               non_sa_basal[region], 100*non_sa_basal[region]/den[region], \
               ambiguous[region], 100*ambiguous[region]/den[region], \
               affects[region], 100*affects[region]/den[region]
    }
}' "$PREFIX.by_window.tsv" > "$PREFIX.any_ingroup_nonmono_summary.tsv"

cat "$PREFIX.any_ingroup_nonmono_summary.tsv"
