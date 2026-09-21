#!/usr/bin/env bash
# Runs qiime2-its and its classifier trainers (train_ncbi, train_fasta) end-to-
# end against small, real, bundled fungal ITS datasets (paired-end, single-end,
# read-length filtering, reverse-complement, a multi-sample/near-empty-sample
# run exercising the advanced diversity/classifier/report steps) to catch
# QIIME2/ITSxpress interface breaks that the mocked-subprocess unit tests
# under tests/ cannot -- those check that *our* code builds the right command,
# not that the command still means the same thing in the QIIME2 version
# actually installed, or that a real taxonomy database produces correct output.
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
# six scenarios above (seconds, bundled ~100KB of data), it downloads tens of
# MB from UNITE and the naive-Bayes fit on the full release can take an hour
# or more of CPU time and several GB of RAM. UNITE data is CC BY-SA 4.0
# (https://unite.ut.ee/cite.php), separate from the CC0/public-domain data
# bundled for the other scenarios (see data/SOURCES.md).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/data"
WITH_UNITE=false
OUT=""
usage () {
    echo "Usage: validation/run_validation.sh [--with-unite] [output_dir]" >&2
    exit 2
}
for arg in "$@"; do
    case "$arg" in
        --with-unite) WITH_UNITE=true ;;
        -*) echo "Unknown option: $arg" >&2; usage ;;
        *) [ -z "$OUT" ] || { echo "More than one output_dir given: $OUT, $arg" >&2; usage; }
           OUT="$arg" ;;
    esac
done
DEFAULT_OUT="$HERE/output"
OUT="${OUT:-$DEFAULT_OUT}"
# Written into every output dir this script creates: the only kind of
# existing, non-empty directory (besides the default one) it agrees to wipe.
SENTINEL=".qiime2_its_validation_output"

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

# output_dir is wiped on every run -- but never a directory this script didn't
# create itself: a mistyped path (~/results, validation/data, a typo'd
# "with-unite") used to be deleted without a question.
if [ -e "$OUT" ] && [ "$OUT" != "$DEFAULT_OUT" ] && [ ! -e "$OUT/$SENTINEL" ] \
        && { [ ! -d "$OUT" ] || [ -n "$(ls -A "$OUT")" ]; }; then
    echo "Refusing to delete $OUT: it exists, is not empty, and was not created by this script" \
         "(no $SENTINEL file in it). Choose another output_dir or remove it yourself." >&2
    exit 1
fi
rm -rf -- "$OUT"
mkdir -p "$OUT"
touch "$OUT/$SENTINEL"

# Exit codes alone have already let wrong output through twice (an
# all-"unidentified" taxonomy; classifiers trained on 4 of 20 reference
# sequences) -- every scenario's *output* is checked too.
check () {
    local description="$1"; shift
    if "$@"; then
        echo "  ok: $description"
    else
        echo "VALIDATION FAILED: $description" >&2
        exit 1
    fi
}
one_taxonomy_line_per_sequence () {
    local n_seqs n_lines
    n_seqs="$(grep -c '^>' "$1")"
    n_lines="$(wc -l < "$2")"
    [ "$n_seqs" -gt 0 ] && [ "$n_seqs" -eq "$n_lines" ]
}
every_asv_classified_to_a_fungal_genus () {
    # The bundled reads all come from the genera the toy classifiers are
    # trained on, so anything short of k__Fungi...g__<Genus> on every ASV
    # means taxonomy building or training regressed.
    local tsv="$1/biom_table/taxonomy.tsv" total resolved
    total="$(($(wc -l < "$tsv") - 1))"
    resolved="$(grep -c 'k__Fungi.*g__[A-Z]' "$tsv" || true)"
    [ "$total" -gt 0 ] && [ "$resolved" -eq "$total" ]
}
any_fungal_asv () {
    grep -q 'k__Fungi' "$1/biom_table/taxonomy.tsv"
}
any_file_matches () {
    compgen -G "$1" >/dev/null
}
check_pipeline_output () {
    local dir="$1"
    check "$dir: report.pdf written" test -s "$dir/report.pdf"
    check "$dir: run_metadata.json written" test -s "$dir/run_metadata.json"
    check "$dir: taxonomy (${TAXONOMY_CHECK:-every_asv_classified_to_a_fungal_genus})" \
        "${TAXONOMY_CHECK:-every_asv_classified_to_a_fungal_genus}" "$dir"
}

echo "== Training a small real classifier (NCBI RefSeq fungal ITS accessions) =="
qiime2-its-train-ncbi \
    -q "$DATA/classifier_training/accessions.list" \
    -o "$OUT/classifier" \
    -t 4 \
    --acc2taxid "$DATA/classifier_training/acc2taxid.tsv" \
    --dead-acc2taxid "$DATA/classifier_training/dead_acc2taxid.tsv" \
    2>&1 | tee "$OUT/train_ncbi.log"

CLASSIFIER="$OUT/classifier/seq_ncbi.qza"
check "train_ncbi: one taxonomy line per downloaded sequence" \
    one_taxonomy_line_per_sequence "$OUT/classifier/seq.fasta" "$OUT/classifier/taxonomy.txt"

echo
echo "== Training a classifier from a fasta + id-table (train_fasta) =="
qiime2-its-train-fasta \
    -q "$DATA/classifier_training/seqs.fasta" \
    -i "$DATA/classifier_training/id_table.tsv" \
    -o "$OUT/classifier_fasta" \
    --taxdump "$OUT/classifier/taxdump.tar.gz" \
    2>&1 | tee "$OUT/train_fasta.log"
check "train_fasta: one taxonomy line per input sequence" \
    one_taxonomy_line_per_sequence "$DATA/classifier_training/seqs.fasta" "$OUT/classifier_fasta/taxonomy.txt"

METADATA="$DATA/metadata.tsv"

# run_case <name> <output_dir> <qiime2-its arguments...>
run_case () {
    local name="$1" case_out="$2"; shift 2
    echo
    echo "== $name =="
    qiime2-its \
        -q "$QIIME2_ENV" \
        -m "$METADATA" \
        -c "$CLASSIFIER" \
        -t 4 -p 2 \
        -o "$case_out" \
        "$@" \
        2>&1 | tee "$OUT/$name.log"
    check_pipeline_output "$case_out"
}

run_case paired_end "$OUT/paired_end" \
    -i "$DATA/paired_end" -pe --extract-its2 --taxa Fungi \
    --sampling-depth 10 --max-rarefaction-depth 60

run_case single_end "$OUT/single_end" \
    -i "$DATA/single_end" -se --extract-its2 --taxa Fungi \
    --max-ee 4 --allow-one-off --sampling-depth 10 --max-rarefaction-depth 60

run_case paired_end_size_filter "$OUT/paired_end_size_filter" \
    -i "$DATA/paired_end" -pe --extract-its2 --taxa Fungi \
    --min-len 120 --max-len 155 --sampling-depth 5 --max-rarefaction-depth 30

# No --extract-its2 here: these ITS3-primed reads are already correctly
# oriented, so reverse-complementing them makes ITSxpress's HMM search fail to
# find ITS boundaries. -rc is validated in isolation (import -> DADA2) rather
# than combined with ITS extraction.
run_case single_end_reverse_complement "$OUT/single_end_reverse_complement" \
    -i "$DATA/single_end" -se -rc --taxa Fungi \
    --max-ee 4 --sampling-depth 5 --max-rarefaction-depth 30

# Swap in the fasta-trained classifier for this one case, to confirm it
# actually classifies (not just that training didn't crash), then restore
# the NCBI one for anything run afterward.
saved_classifier="$CLASSIFIER"
CLASSIFIER="$OUT/classifier_fasta/naive-bayes_classifier.qza"
run_case fasta_classifier_classification "$OUT/fasta_classifier_pipeline" \
    -i "$DATA/single_end" -se --extract-its2 --taxa Fungi \
    --max-ee 4 --allow-one-off --sampling-depth 10 --max-rarefaction-depth 60
CLASSIFIER="$saved_classifier"

# 6 real samples (5 healthy + 1 near-empty, siteC-rep2, which DADA2 reduces
# to a zero-read row) with a realistic multi-column metadata file, so the
# advanced-stats/report steps (on by default) get exercised against >2
# samples and the near-empty-sample edge case, not just the 2-sample/
# 1-column metadata.tsv used above.
saved_metadata="$METADATA"
METADATA="$DATA/metadata_multi.tsv"
run_case multi_sample_advanced_stats "$OUT/multi_sample_advanced_stats" \
    -i "$DATA/multi_sample" -pe --extract-its2 --taxa Fungi \
    --sampling-depth 10 --max-rarefaction-depth 15
METADATA="$saved_metadata"
check "multi_sample_advanced_stats: beta group-significance results written" \
    any_file_matches "$OUT/multi_sample_advanced_stats/beta-group-significance-*.qzv"
check "multi_sample_advanced_stats: alpha group-significance results written" \
    any_file_matches "$OUT/multi_sample_advanced_stats/alpha-group-significance-*.qzv"

run_unite_scenario () {
    # UNITE's own site now gates release files behind DOI/PlutoF landing
    # pages with no stable direct-download URL (see README.md), so this
    # fetches via the rescript plugin instead, then repackages the result
    # into the "sh_qiime_release_<date>/developer/..." layout a real UNITE
    # QIIME-release archive has, so it exercises qiime2-its-train-unite's own
    # extraction/file-matching/training code (not just rescript).
    local version="${UNITE_VERSION:-2025-02-19}"
    local taxon_group="${UNITE_TAXON_GROUP:-fungi}"
    # Fixed, not overridable: qiime2-its-train-unite only picks up the 99%
    # files from a UNITE archive.
    local cluster_id="99"
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

    # A real, full-size reference database legitimately leaves some ASVs
    # short of genus -- only require that classification produced fungi.
    TAXONOMY_CHECK=any_fungal_asv run_case unite_classification "$OUT/unite_pipeline" \
        -i "$DATA/single_end" -se --extract-its2 --taxa Fungi \
        --max-ee 4 --allow-one-off --sampling-depth 10 --max-rarefaction-depth 60

    CLASSIFIER="$saved_classifier"
}

if $WITH_UNITE; then
    run_unite_scenario
fi

echo
echo "All validation scenarios completed successfully."
echo "Logs: $OUT/*.log"
