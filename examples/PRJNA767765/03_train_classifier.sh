#!/usr/bin/env bash
# Step 3: build the taxonomic classifier the example was run with -- UNITE
# (release 2025-02-19, fungi, 99% clusters, no singletons), fetched with
# QIIME2's own `rescript` plugin and trained as a naive-Bayes classifier.
#
# Slow and memory-hungry: on the order of an hour of CPU time and several GB
# of RAM. If you already have a UNITE classifier for your QIIME2 release,
# skip this step and point 04_run_pipeline.sh at it instead.
#
# UNITE is distributed under CC BY-SA 4.0 -- cite it: https://unite.ut.ee/cite.php
#
# Usage (inside your activated QIIME2 conda environment):
#     03_train_classifier.sh [classifier_dir]    (default: ./classifier)
set -euo pipefail

OUT="${1:-classifier}"
mkdir -p "$OUT"

qiime rescript get-unite-data \
    --p-version 2025-02-19 \
    --p-taxon-group fungi \
    --p-cluster-id 99 \
    --p-no-singletons \
    --o-sequences "$OUT/unite-sequences.qza" \
    --o-taxonomy "$OUT/unite-taxonomy.qza"

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$OUT/unite-sequences.qza" \
    --i-reference-taxonomy "$OUT/unite-taxonomy.qza" \
    --o-classifier "$OUT/unite-classifier.qza"

echo "Classifier: $OUT/unite-classifier.qza"
