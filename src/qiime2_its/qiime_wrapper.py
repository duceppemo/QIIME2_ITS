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


def dada2_denoise_single(reads_qza, repseq_qza, table_qza, stats_qza, n_threads=0):
    cmd = ['qiime', 'dada2', 'denoise-single',
           '--p-n-threads', str(n_threads),
           '--p-trim-left', '0',
           '--p-trunc-len', '0',
           '--i-demultiplexed-seqs', str(reads_qza),
           '--o-representative-sequences', str(repseq_qza),
           '--o-table', str(table_qza),
           '--o-denoising-stats', str(stats_qza)]
    _run(cmd)


def dada2_denoise_paired(reads_qza, repseq_qza, table_qza, stats_qza, n_threads=0):
    cmd = ['qiime', 'dada2', 'denoise-paired',
           '--p-n-threads', str(n_threads),
           '--p-trim-left-f', '0',
           '--p-trim-left-r', '0',
           '--p-trunc-len-f', '0',
           '--p-trunc-len-r', '0',
           '--i-demultiplexed-seqs', str(reads_qza),
           '--o-representative-sequences', str(repseq_qza),
           '--o-table', str(table_qza),
           '--o-denoising-stats', str(stats_qza)]
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


def sample_summarize(metadata_file, table_qza, table_qzv):
    cmd = ['qiime', 'feature-table', 'summarize',
           '--m-sample-metadata-file', str(metadata_file),
           '--i-table', str(table_qza),
           '--o-visualization', str(table_qzv)]
    _run(cmd)


def seq_summary(repseqs_qza, repseqs_qzv):
    cmd = ['qiime', 'feature-table', 'tabulate-seqs',
           '--i-data', str(repseqs_qza),
           '--o-visualization', str(repseqs_qzv)]
    _run(cmd)


def phylogeny(repseqs_qza, align_repseqs_qza, masked_align_repseqs_qza, unrooted_tree_qza, rooted_tree_qza):
    cmd = ['qiime', 'phylogeny', 'align-to-tree-mafft-fasttree',
           '--p-n-threads', 'auto',
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


def classify(classifier_qza, repseqs_qza, taxonomy_qza, n_jobs=-1):
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
