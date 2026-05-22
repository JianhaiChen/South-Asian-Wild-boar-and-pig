#!/usr/bin/env bash
set -euo pipefail

DEEPS="3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 8.0 10.0"

for D in ${DEEPS}; do
  TAG=$(echo "${D}" | sed 's/\./p/g')
  echo "===== deep-mya ${D} ====="

  python simulate_ils_deep_tracts.py \
    --n-reps 1000 \
    --length-mb 43 \
    --observed-mb 40 \
    --generation-time 3 \
    --split-kya 695 \
    --deep-mya "${D}" \
    --ne 20000 \
    --ne-anc-list 500000 \
    --recomb-cm-mb-list 0.01 \
    --out "ils_X43_obs40_rec0p01_NeAnc500k_deep${TAG}"
done

echo "[DONE] all deep cutoff simulations finished"