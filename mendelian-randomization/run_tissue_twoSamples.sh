#!/bin/bash -l

set -e
set -u
set -o pipefail

tissue=$1

nextflow run mr_twoSamples.nf \
  --gt_gex_stats OUT/${tissue}/gt_gex_stats.tsv \
  --out OUT/${tissue} \
  -with-report "report_${tissue}_2s.html" \
  -with-dag "dag_${tissue}_2s.png" \
  -work-dir "work/${tissue}" \
  -resume
