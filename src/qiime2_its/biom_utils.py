"""Helpers for merging QIIME2 taxonomy assignments into a biom feature table."""
import subprocess
from pathlib import Path

TAXONOMY_HEADER = '#OTUID\ttaxonomy\tconfidence\n'


def rewrite_taxonomy_header(taxonomy_tsv):
    """Replace the first line of a QIIME2-exported taxonomy.tsv with the header
    biom's `add-metadata` expects (#OTUID / taxonomy / confidence).
    """
    taxonomy_tsv = Path(taxonomy_tsv)
    tmp_path = taxonomy_tsv.with_name(taxonomy_tsv.name + '.tmp')
    with open(taxonomy_tsv) as in_f, open(tmp_path, 'w') as out_f:
        out_f.write(TAXONOMY_HEADER)
        next(in_f, None)  # skip original header
        out_f.writelines(in_f)
    tmp_path.replace(taxonomy_tsv)


def add_metadata(input_biom, taxonomy_tsv, output_biom):
    cmd = ['biom', 'add-metadata',
           '--sc-separated', 'taxonomy',
           '-i', str(input_biom),
           '--observation-metadata-fp', str(taxonomy_tsv),
           '-o', str(output_biom)]
    subprocess.run(cmd, check=True)


def convert_to_tsv(input_biom, taxonomy_tsv, output_tsv):
    cmd = ['biom', 'convert',
           '--to-tsv',
           '--header-key', 'taxonomy',
           '-i', str(input_biom),
           '--observation-metadata-fp', str(taxonomy_tsv),
           '-o', str(output_tsv)]
    subprocess.run(cmd, check=True)
