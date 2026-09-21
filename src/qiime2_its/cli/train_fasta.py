"""Build a QIIME2 classifier from a user-supplied fasta file + accession-to-taxid table."""
import argparse
from pathlib import Path
from time import time

from qiime2_its import downloader, env_checks, qiime_wrapper, taxonomy, timing
from qiime2_its._version import __version__

TAXDUMP_URL = 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'


def run(fasta_query, id_table, output_folder, taxdump):
    env_checks.check_qiime2_env_active()

    fasta_query = Path(fasta_query)
    id_table = Path(id_table)
    output_folder = Path(output_folder)

    if not fasta_query.is_file():
        raise ValueError('Your "query" is not a file or does not exist.')
    if not id_table.is_file():
        raise ValueError('Your "id-table" is not a file or does not exist.')

    # Checked up front (before a taxdump download and a QIIME2 import): every
    # sequence needs a taxonomy line, and a mismatch otherwise only surfaces
    # as an error from deep inside `qiime feature-classifier`.
    id_dict = taxonomy.parse_id_table(id_table)  # {accession: taxid}
    fasta_ids = taxonomy.read_fasta_ids(fasta_query)
    if not fasta_ids:
        raise ValueError(f'No sequences found in {fasta_query}.')
    missing = [seq_id for seq_id in fasta_ids if seq_id not in id_dict]
    if missing:
        shown = ', '.join(missing[:10]) + (f', ... ({len(missing)} in total)' if len(missing) > 10 else '')
        raise ValueError(f'These sequence IDs from {fasta_query.name} are missing from the id-table '
                         f'{id_table.name}: {shown}. IDs (everything before the first whitespace in each '
                         f'fasta header) must match the table\'s first column exactly.')
    output_folder.mkdir(parents=True, exist_ok=True)

    t_zero = time()

    start = time()
    if taxdump:
        print('Extracting taxdump.tar.gz...', end='', flush=True)
        downloader.extract_targz(taxdump, output_folder)
    else:
        print('Downloading taxdump.tar.gz...', end='', flush=True)
        downloader.download(TAXDUMP_URL, output_folder / 'taxdump.tar.gz')
        downloader.extract_targz(output_folder / 'taxdump.tar.gz', output_folder)
    print(f' took {timing.format_elapsed(time() - start)}')

    start = time()
    print('Writing taxonomy...', end='', flush=True)
    taxonomy_file = output_folder / 'taxonomy.txt'
    taxonomy.write_taxonomy_file(id_dict, taxonomy_file, output_folder / 'nodes.dmp',
                                  output_folder / 'names.dmp', output_folder / 'merged.dmp')
    print(f' took {timing.format_elapsed(time() - start)}')

    qiime2_seq = output_folder / (fasta_query.stem + '.qza')
    qiime2_taxo = taxonomy_file.with_suffix(taxonomy_file.suffix + '.qza')
    classifier_file = output_folder / 'naive-bayes_classifier.qza'

    print('Importing sequences into QIIME2...')
    qiime_wrapper.import_sequences(fasta_query, qiime2_seq)
    print('Importing taxonomy into QIIME2...')
    qiime_wrapper.import_taxonomy(taxonomy_file, qiime2_taxo)
    print('Training classifier...')
    qiime_wrapper.train_naive_bayes_classifier(qiime2_seq, qiime2_taxo, classifier_file, verbose=True)

    print(f'Done (total time: {timing.format_elapsed(time() - t_zero)}).')


def build_parser():
    parser = argparse.ArgumentParser(
        description='Prep a QIIME2 classifier from a fasta file and corresponding "acc to taxid" table.')
    parser.add_argument('-q', '--query', metavar='my_sequences.fasta', required=True, type=str,
                         help='A fasta file. Mandatory.')
    parser.add_argument('-i', '--id-table', metavar='/path/to/acc2taxid_table.tsv', required=True, type=str,
                         help='Tab-separated text file with 2 columns (accession, taxid) matching the input '
                              'fasta. Accession numbers (everything before the first whitespace in each fasta '
                              'header) must match exactly between the fasta and the table. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/path/to/output_folder/', required=True, type=str,
                         help='Output folder. Mandatory.')
    parser.add_argument('--taxdump', metavar='/path/to/taxdump.tar.gz', default=None, type=str,
                         help='Path to a downloaded taxdump.tar.gz. Downloaded otherwise. Optional.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser


def main():
    args = build_parser().parse_args()
    run(args.query, args.id_table, args.output, args.taxdump)


if __name__ == '__main__':
    main()
