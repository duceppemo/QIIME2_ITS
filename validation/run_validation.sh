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
# Usage: validation/run_validation.sh [--with-unite] [output_dir]
#
# --with-unite additionally downloads a real, current UNITE release (via the
# `rescript` plugin, no bundled data -- see the README) and trains and
# classifies with it. It's opt-in and not part of the default run: unlike the
# four scenarios above (seconds, bundled ~100KB of data), it downloads tens of
# MB from UNITE and the naive-Bayes fit on the full release can take an hour
# or more of CPU time and several GB of RAM. UNITE data is CC BY-SA 4.0
# (https://unite.ut.ee/cite.php), separate from the CC0/public-domain data
# bundled for the other scenarios (see data/SOURCES.md).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/data"
WITH_UNITE=false
OUT=""
for arg in "$@"; do
    case "$arg" in
        --with-unite) WITH_UNITE=true ;;
        *) OUT="$arg" ;;
    esac
done
OUT="${OUT:-$HERE/output}"

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
if $WITH_UNITE; then
    qiime rescript get-unite-data --help >/dev/null 2>&1 || {
        echo "qiime rescript get-unite-data not found. The rescript plugin ships with QIIME2" \
             "by default -- check 'qiime info' if it's missing." >&2
        exit 1
    }
fi

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

run_unite_scenario () {
    # UNITE's own site now gates release files behind DOI/PlutoF landing
    # pages with no stable direct-download URL (see README.md), so this
    # fetches via the rescript plugin instead, then repackages the result
    # into the "sh_qiime_release_<date>/developer/..." layout a real UNITE
    # QIIME-release archive has, so it exercises qiime2-its-train-unite's own
    # extraction/file-matching/training code (not just rescript).
    local version="${UNITE_VERSION:-2025-02-19}"
    local taxon_group="${UNITE_TAXON_GROUP:-fungi}"
    local cluster_id="${UNITE_CLUSTER_ID:-99}"
    local yyyy mm dd release_date
    IFS='-' read -r yyyy mm dd <<< "$version"
    release_date="${dd}.${mm}.${yyyy}"

    echo
    echo "== unite (heavy: real UNITE $version download, full training, real classification) =="
    local unite_dir="$OUT/unite"
    mkdir -p "$unite_dir"

    echo "Downloading UNITE $version ($taxon_group, ${cluster_id}%, no singletons) via rescript..."
    qiime rescript get-unite-data \
        --p-version "$version" \
        --p-taxon-group "$taxon_group" \
        --p-cluster-id "$cluster_id" \
        --p-no-singletons \
        --o-sequences "$unite_dir/sequences.qza" \
        --o-taxonomy "$unite_dir/taxonomy.qza" \
        2>&1 | tee "$OUT/unite_download.log"

    qiime tools export --input-path "$unite_dir/sequences.qza" --output-path "$unite_dir/seq_export" >/dev/null
    qiime tools export --input-path "$unite_dir/taxonomy.qza" --output-path "$unite_dir/taxo_export" >/dev/null

    local developer="$unite_dir/sh_qiime_release_${release_date}/developer"
    mkdir -p "$developer"
    cp "$unite_dir/seq_export/dna-sequences.fasta" \
        "$developer/sh_refs_qiime_ver10_${cluster_id}_${release_date}.fasta"
    tail -n +2 "$unite_dir/taxo_export/taxonomy.tsv" \
        > "$developer/sh_taxonomy_qiime_ver10_${cluster_id}_${release_date}.txt"

    local archive="$unite_dir/unite_ver10_${cluster_id}_${release_date}.tgz"
    tar -C "$unite_dir" -czf "$archive" "sh_qiime_release_${release_date}"

    echo "Training the UNITE classifier (slow: a naive-Bayes fit on the full release," \
         "expect on the order of an hour or more of CPU time and several GB of RAM)..."
    qiime2-its-train-unite \
        -u "$archive" \
        -o "$OUT/unite_classifier" \
        -q "$QIIME2_ENV" \
        2>&1 | tee "$OUT/unite_train.log"

    # Swap in the just-trained real classifier for this one case, in place of
    # the small NCBI toy classifier run_case otherwise uses.
    local saved_classifier="$CLASSIFIER"
    CLASSIFIER="$(find "$OUT/unite_classifier" -maxdepth 1 -name 'unite-*-classifier-*.qza' | head -1)"

    run_case unite_classification \
        -i "$DATA/single_end" -o "$OUT/unite_pipeline" -se --extract-its2 --taxa Fungi \
        --max-ee 4 --allow-one-off --sampling-depth 10 --max-rarefaction-depth 60

    CLASSIFIER="$saved_classifier"
}

if $WITH_UNITE; then
    run_unite_scenario
fi

echo
echo "All validation scenarios completed successfully."
echo "Logs: $OUT/*.log"
