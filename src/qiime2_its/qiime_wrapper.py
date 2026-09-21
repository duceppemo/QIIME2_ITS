"""Thin, individually-testable wrappers around the `qiime` CLI.

Every function builds one `qiime ...` command and runs it with
`subprocess.run(cmd, check=True)` so a failed step stops the pipeline instead
of being silently ignored.
"""
import subprocess


def _run(cmd):
    subprocess.run(cmd, check=True)


# --- Import / sequence & taxonomy prep, shared by all classifier trainers ---

def import_sequences(sequence_file, output_qza):
    """`qiime tools import` a FeatureData[Sequence] fasta file."""
    cmd = ['qiime', 'tools', 'import',
           '--type', 'FeatureData[Sequence]',
           '--input-path', str(sequence_file),
           '--output-path', str(output_qza)]
    _run(cmd)


def import_taxonomy(taxonomy_file, output_qza):
    """`qiime tools import` a headerless TSV FeatureData[Taxonomy] file."""
    cmd = ['qiime', 'tools', 'import',
           '--type', 'FeatureData[Taxonomy]',
           '--input-format', 'HeaderlessTSVTaxonomyFormat',
           '--input-path', str(taxonomy_file),
           '--output-path', str(output_qza)]
    _run(cmd)


def train_naive_bayes_classifier(sequence_qza, taxonomy_qza, classifier_qza, verbose=False):
    cmd = ['qiime', 'feature-classifier', 'fit-classifier-naive-bayes',
           '--i-reference-reads', str(sequence_qza),
           '--i-reference-taxonomy', str(taxonomy_qza),
           '--o-classifier', str(classifier_qza)]
    if verbose:
        cmd.append('--verbose')
    _run(cmd)


# --- Main ITS pipeline ---

def import_fastq_se(fastq_folder, reads_qza):
    """Import demultiplexed single-end fastq (Casava layout)."""
    cmd = ['qiime', 'tools', 'import',
           '--type', 'SampleData[SequencesWithQuality]',
           '--input-format', 'CasavaOneEightSingleLanePerSampleDirFmt',
           '--input-path', str(fastq_folder),
           '--output-path', str(reads_qza)]
    _run(cmd)


def import_fastq_pe(fastq_folder, reads_qza):
    """Import demultiplexed paired-end fastq (Casava layout)."""
    cmd = ['qiime', 'tools', 'import',
           '--type', 'SampleData[PairedEndSequencesWithQuality]',
           '--input-format', 'CasavaOneEightSingleLanePerSampleDirFmt',
           '--input-path', str(fastq_folder),
           '--output-path', str(reads_qza)]
    _run(cmd)


def demux_summary(reads_qza, output_qzv, n=1000):
    cmd = ['qiime', 'demux', 'summarize',
           '--p-n', str(n),
           '--i-data', str(reads_qza),
           '--o-visualization', str(output_qzv)]
    _run(cmd)


# Defaults match `qiime dada2 denoise-single/-paired`'s own defaults, tuned for
# Illumina data. Override max_ee/allow_one_off in particular for noisier
# single-end platforms (e.g. IonTorrent): its higher per-base error rate and
# homopolymer-driven indels mean more legitimate reads get discarded by the
# default max-expected-errors threshold, and more true ASVs get miscalled as
# one-off bimeras by the default (non-one-off) chimera search.
def dada2_denoise_single(reads_qza, repseq_qza, table_qza, stats_qza, base_transition_stats_qza, n_threads=0,
                          max_ee=2.0, trunc_q=2, pooling_method='independent', chimera_method='consensus',
                          min_fold_parent_over_abundance=1.0, allow_one_off=False, n_reads_learn=1000000):
    cmd = ['qiime', 'dada2', 'denoise-single',
           '--p-n-threads', str(n_threads),
           '--p-trim-left', '0',
           '--p-trunc-len', '0',
           '--p-max-ee', str(max_ee),
           '--p-trunc-q', str(trunc_q),
           '--p-pooling-method', pooling_method,
           '--p-chimera-method', chimera_method,
           '--p-min-fold-parent-over-abundance', str(min_fold_parent_over_abundance),
           '--p-allow-one-off' if allow_one_off else '--p-no-allow-one-off',
           '--p-n-reads-learn', str(n_reads_learn),
           '--i-demultiplexed-seqs', str(reads_qza),
           '--o-representative-sequences', str(repseq_qza),
           '--o-table', str(table_qza),
           '--o-denoising-stats', str(stats_qza),
           '--o-base-transition-stats', str(base_transition_stats_qza)]
    _run(cmd)


def dada2_denoise_paired(reads_qza, repseq_qza, table_qza, stats_qza, base_transition_stats_qza, n_threads=0,
                          max_ee_f=2.0, max_ee_r=2.0, trunc_q=2, pooling_method='independent',
                          chimera_method='consensus', min_fold_parent_over_abundance=1.0, allow_one_off=False,
                          n_reads_learn=1000000):
    cmd = ['qiime', 'dada2', 'denoise-paired',
           '--p-n-threads', str(n_threads),
           '--p-trim-left-f', '0',
           '--p-trim-left-r', '0',
           '--p-trunc-len-f', '0',
           '--p-trunc-len-r', '0',
           '--p-max-ee-f', str(max_ee_f),
           '--p-max-ee-r', str(max_ee_r),
           '--p-trunc-q', str(trunc_q),
           '--p-pooling-method', pooling_method,
           '--p-chimera-method', chimera_method,
           '--p-min-fold-parent-over-abundance', str(min_fold_parent_over_abundance),
           '--p-allow-one-off' if allow_one_off else '--p-no-allow-one-off',
           '--p-n-reads-learn', str(n_reads_learn),
           '--i-demultiplexed-seqs', str(reads_qza),
           '--o-representative-sequences', str(repseq_qza),
           '--o-table', str(table_qza),
           '--o-denoising-stats', str(stats_qza),
           '--o-base-transition-stats', str(base_transition_stats_qza)]
    _run(cmd)


def metadata_tabulate(input_qza, output_qzv):
    cmd = ['qiime', 'metadata', 'tabulate',
           '--m-input-file', str(input_qza),
           '--o-visualization', str(output_qzv)]
    _run(cmd)


def export(qza, output_folder):
    cmd = ['qiime', 'tools', 'export',
           '--input-path', str(qza),
           '--output-path', str(output_folder)]
    _run(cmd)


def sample_summarize(metadata_file, table_qza, table_qzv, feature_frequencies_qza, sample_frequencies_qza):
    cmd = ['qiime', 'feature-table', 'summarize',
           '--m-metadata-file', str(metadata_file),
           '--i-table', str(table_qza),
           '--o-summary', str(table_qzv),
           '--o-feature-frequencies', str(feature_frequencies_qza),
           '--o-sample-frequencies', str(sample_frequencies_qza)]
    _run(cmd)


def seq_summary(repseqs_qza, repseqs_qzv):
    cmd = ['qiime', 'feature-table', 'tabulate-seqs',
           '--i-data', str(repseqs_qza),
           '--o-visualization', str(repseqs_qzv)]
    _run(cmd)


def phylogeny(repseqs_qza, align_repseqs_qza, masked_align_repseqs_qza, unrooted_tree_qza, rooted_tree_qza,
              n_threads='auto'):
    """`n_threads`: a thread count, or 'auto' for every available core."""
    cmd = ['qiime', 'phylogeny', 'align-to-tree-mafft-fasttree',
           '--p-n-threads', str(n_threads),
           '--i-sequences', str(repseqs_qza),
           '--o-alignment', str(align_repseqs_qza),
           '--o-masked-alignment', str(masked_align_repseqs_qza),
           '--o-tree', str(unrooted_tree_qza),
           '--o-rooted-tree', str(rooted_tree_qza)]
    _run(cmd)


def core_diversity(cpu, metadata_file, rooted_tree_qza, table_qza, output_folder, sampling_depth=1000):
    cmd = ['qiime', 'diversity', 'core-metrics-phylogenetic',
           '--p-n-jobs-or-threads', str(cpu),
           '--p-sampling-depth', str(sampling_depth),
           '--i-phylogeny', str(rooted_tree_qza),
           '--i-table', str(table_qza),
           '--m-metadata-file', str(metadata_file),
           '--output-dir', f'{output_folder}/core-metrics-results']
    _run(cmd)


def rarefaction(metadata_file, rooted_tree_qza, table_qza, rare_qzv, max_depth=4000):
    cmd = ['qiime', 'diversity', 'alpha-rarefaction',
           '--p-max-depth', str(max_depth),
           '--i-phylogeny', str(rooted_tree_qza),
           '--i-table', str(table_qza),
           '--m-metadata-file', str(metadata_file),
           '--o-visualization', str(rare_qzv)]
    _run(cmd)


def classify(classifier_qza, repseqs_qza, taxonomy_qza, n_jobs=0):
    """`n_jobs`: 0 uses all CPUs, 1 disables parallelism (q2-feature-classifier's
    own convention -- unlike scikit-learn, it does not accept -1 for "all")."""
    cmd = ['qiime', 'feature-classifier', 'classify-sklearn',
           '--p-n-jobs', str(n_jobs),
           '--i-classifier', str(classifier_qza),
           '--i-reads', str(repseqs_qza),
           '--o-classification', str(taxonomy_qza)]
    _run(cmd)


def taxa_barplot(table_qza, taxonomy_qza, metadata_file, output_qzv):
    cmd = ['qiime', 'taxa', 'barplot',
           '--i-table', str(table_qza),
           '--i-taxonomy', str(taxonomy_qza),
           '--m-metadata-file', str(metadata_file),
           '--o-visualization', str(output_qzv)]
    _run(cmd)


# --- Group-significance, taxonomy collapse, and sample classification ---

def alpha_group_significance(alpha_qza, metadata_file, output_qzv):
    """Kruskal-Wallis test of one alpha-diversity vector against every
    categorical column in `metadata_file` (all in one call)."""
    cmd = ['qiime', 'diversity', 'alpha-group-significance',
           '--i-alpha-diversity', str(alpha_qza),
           '--m-metadata-file', str(metadata_file),
           '--o-visualization', str(output_qzv)]
    _run(cmd)


def beta_group_significance(distance_matrix_qza, metadata_file, column, output_qzv, method='permanova'):
    """PERMANOVA (by default) test of one distance matrix against one
    categorical metadata column."""
    cmd = ['qiime', 'diversity', 'beta-group-significance',
           '--i-distance-matrix', str(distance_matrix_qza),
           '--m-metadata-file', str(metadata_file),
           '--m-metadata-column', column,
           '--p-method', method,
           '--o-visualization', str(output_qzv)]
    _run(cmd)


def taxa_collapse(table_qza, taxonomy_qza, level, output_qza):
    """Collapse a feature table to the given taxonomic level (6 = genus, 7 = species
    in this pipeline's k;p;c;o;f;g;s lineage strings)."""
    cmd = ['qiime', 'taxa', 'collapse',
           '--i-table', str(table_qza),
           '--i-taxonomy', str(taxonomy_qza),
           '--p-level', str(level),
           '--o-collapsed-table', str(output_qza)]
    _run(cmd)


def relative_frequency(table_qza, output_qza):
    cmd = ['qiime', 'feature-table', 'relative-frequency',
           '--i-table', str(table_qza),
           '--o-relative-frequency-table', str(output_qza)]
    _run(cmd)


def classify_samples(table_qza, metadata_file, column, output_dir, cv, n_estimators=100):
    """Train a random-forest classifier predicting `column` from `table_qza`.

    `cv` should be capped by the caller to what the smallest class in
    `column` can support (scikit-learn's stratified k-fold cross-validation
    requires cv <= the smallest class size).
    """
    cmd = ['qiime', 'sample-classifier', 'classify-samples',
           '--i-table', str(table_qza),
           '--m-metadata-file', str(metadata_file),
           '--m-metadata-column', column,
           '--p-cv', str(cv),
           '--p-n-estimators', str(n_estimators),
           '--output-dir', str(output_dir)]
    _run(cmd)


def qiime_info():
    """Returns `qiime info`'s raw stdout: framework/Python versions and every
    installed plugin's version. Used for the report's QA/provenance pages
    (provenance.parse_qiime_info() does the actual parsing)."""
    return subprocess.run(['qiime', 'info'], check=True, capture_output=True, text=True).stdout
