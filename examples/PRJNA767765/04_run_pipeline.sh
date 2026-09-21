#!/usr/bin/env bash
# Step 4: the actual analysis -- one qiime2-its command.
#
#   -pe                    paired-end reads
#   --extract-its2         the amplicon covers ITS2: ITSxpress trims the
#                          conserved 5.8S/LSU flanks before DADA2
#   --taxa Fungi           ITSxpress HMM profile
#   defaults otherwise     (DADA2 parameters, sampling depth 1000, advanced
#                          statistics and the PDF report all on)
#
# Run it from the folder holding raw_reads/ so the paths recorded in the
# report's provenance pages are the short relative ones below.
#
# Usage (inside your activated QIIME2 conda environment):
#     THREADS=40 PARALLEL=10 04_run_pipeline.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# ./metadata.tsv if step 1 was run here, else the copy committed next to this script
if [ -z "${METADATA:-}" ]; then
    if [ -f metadata.tsv ]; then METADATA=metadata.tsv; else METADATA="$HERE/metadata.tsv"; fi
fi
CLASSIFIER="${CLASSIFIER:-classifier/unite-classifier.qza}"
THREADS="${THREADS:-8}"      # total CPUs
PARALLEL="${PARALLEL:-4}"    # samples processed at once (per-sample steps)

qiime2-its \
    -q "$CONDA_DEFAULT_ENV" \
    -i raw_reads \
    -o output \
    -m "$METADATA" \
    -c "$CLASSIFIER" \
    -pe --extract-its2 --taxa Fungi \
    -t "$THREADS" -p "$PARALLEL"
