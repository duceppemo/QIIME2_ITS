#!/usr/bin/env bash
# Runs qiime2-its end-to-end against small, real, bundled fungal ITS datasets
# (paired-end, single-end, read-length filtering, reverse-complement) to catch
# QIIME2/ITSxpress interface breaks that the mocked-subprocess unit tests
# under tests/ cannot -- those check that *our* code builds the right command,
# not that the command still means the same thing in the QIIME2 version
# actually installed.
#
# Requires: an activated QIIME2 conda environment with this package installed
# (`pip install -e .` from the repo root) and bbduk.sh (BBTools/BBMap) on
# PATH, e.g. `conda install -c bioconda -c conda-forge bbmap` -- the
# paired-end size-filter case below exercises --min-len/--max-len, which
# depends on it.
#
# Usage: validation/run_validation.sh [output_dir]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/data"
OUT="${1:-$HERE/output}"

command -v qiime2-its >/dev/null || {
    echo "qiime2-its not found on PATH. Activate your QIIME2 conda env and run" \
         "'pip install -e .' from the repo root first." >&2
    exit 1
}
command -v qiime >/dev/null || {
    echo "qiime not found on PATH. Activate your QIIME2 conda environment first." >&2
    exit 1
}
command -v bbduk.sh >/dev/null || {
    echo "bbduk.sh not found on PATH. The paired_end_size_filter case below needs it:" \
         "'conda install -c bioconda -c conda-forge bbmap'." >&2
    exit 1
}

QIIME2_ENV="${CONDA_DEFAULT_ENV:-qiime2}"
echo "QIIME2 env : $QIIME2_ENV"
echo "qiime      : $(qiime --version | head -1)"
echo "itsxpress  : $(qiime itsxpress --version 2>&1 | head -1)"
echo "output dir : $OUT"
echo

rm -rf "$OUT"
mkdir -p "$OUT"

echo "== Training a small real classifier (NCBI RefSeq fungal ITS accessions) =="
qiime2-its-train-ncbi \
    -q "$DATA/classifier_training/accessions.list" \
    -o "$OUT/classifier" \
    -t 4 \
    --acc2taxid "$DATA/classifier_training/acc2taxid.tsv" \
    --dead-acc2taxid "$DATA/classifier_training/dead_acc2taxid.tsv" \
    2>&1 | tee "$OUT/train_ncbi.log"

CLASSIFIER="$OUT/classifier/seq_ncbi.qza"

run_case () {
    local name="$1"; shift
    echo
    echo "== $name =="
    qiime2-its \
        -q "$QIIME2_ENV" \
        -m "$DATA/metadata.tsv" \
        -c "$CLASSIFIER" \
        -t 4 -p 2 \
        "$@" \
        2>&1 | tee "$OUT/$name.log"
}

run_case paired_end \
    -i "$DATA/paired_end" -o "$OUT/paired_end" -pe --extract-its2 --taxa Fungi \
    --sampling-depth 10 --max-rarefaction-depth 60

run_case single_end \
    -i "$DATA/single_end" -o "$OUT/single_end" -se --extract-its2 --taxa Fungi \
    --max-ee 4 --allow-one-off --sampling-depth 10 --max-rarefaction-depth 60

run_case paired_end_size_filter \
    -i "$DATA/paired_end" -o "$OUT/paired_end_size_filter" -pe --extract-its2 --taxa Fungi \
    --min-len 120 --max-len 155 --sampling-depth 5 --max-rarefaction-depth 30

# No --extract-its2 here: these ITS3-primed reads are already correctly
# oriented, so reverse-complementing them makes ITSxpress's HMM search fail to
# find ITS boundaries. -rc is validated in isolation (import -> DADA2) rather
# than combined with ITS extraction.
run_case single_end_reverse_complement \
    -i "$DATA/single_end" -o "$OUT/single_end_reverse_complement" -se -rc --taxa Fungi \
    --max-ee 4 --sampling-depth 5 --max-rarefaction-depth 30

echo
echo "All validation scenarios completed successfully."
echo "Logs: $OUT/*.log"
