"""Wrapper around the `qiime itsxpress` plugin (ITSxpress >= 2.0).

ITSxpress v2 merged the old standalone `itsxpress` CLI into a QIIME2 plugin.
It now operates on imported `.qza` artifacts rather than raw fastq files, so
callers must `qiime tools import` first (see qiime_wrapper.import_fastq_se/pe)
and can feed the trimmed output straight into DADA2 or re-export it to fastq
for further processing (e.g. size filtering, empty-record cleanup).
"""
import subprocess

# CLI-friendly taxon names (as used by this pipeline's own --taxa flag) mapped
# to the single-letter codes the `qiime itsxpress` plugin expects.
TAXA_CODES = {
    'Alveolata': 'A',
    'Bryophyta': 'B',
    'Bacillariophyta': 'C',
    'Amoebozoa': 'D',
    'Euglenozoa': 'E',
    'Fungi': 'F',
    'Chlorophyta': 'G',
    'Rhodophyta': 'H',
    'Phaeophyceae': 'I',
    'Marchantiophyta': 'L',
    'Metazoa': 'M',
    'Oomycota': 'O',
    'Haptophyceae': 'P',
    'Raphidophyceae': 'Q',
    'Rhizaria': 'R',
    'Synurophyceae': 'S',
    'Tracheophyta': 'T',
    'Eustigmatophyceae': 'U',
    'Parabasalia': 'Y',
    'All': 'ALL',
}


def _run(cmd):
    subprocess.run(cmd, check=True)


def trim_single(reads_qza, trimmed_qza, region, taxa, threads=1, cluster_id=1.0):
    """`qiime itsxpress trim-single`. `taxa` is a full name (see TAXA_CODES)."""
    cmd = ['qiime', 'itsxpress', 'trim-single',
           '--i-per-sample-sequences', str(reads_qza),
           '--p-region', region,
           '--p-taxa', TAXA_CODES[taxa],
           '--p-threads', str(threads),
           '--p-cluster-id', str(cluster_id),
           '--o-trimmed', str(trimmed_qza)]
    _run(cmd)


def trim_pair_unmerged(reads_qza, trimmed_qza, region, taxa, threads=1, cluster_id=1.0):
    """`qiime itsxpress trim-pair-output-unmerged` (unmerged output, for DADA2).

    `taxa` is a full name (see TAXA_CODES).
    """
    cmd = ['qiime', 'itsxpress', 'trim-pair-output-unmerged',
           '--i-per-sample-sequences', str(reads_qza),
           '--p-region', region,
           '--p-taxa', TAXA_CODES[taxa],
           '--p-threads', str(threads),
           '--p-cluster-id', str(cluster_id),
           '--o-trimmed', str(trimmed_qza)]
    _run(cmd)
