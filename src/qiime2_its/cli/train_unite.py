"""Train a QIIME2 naive-Bayes classifier from a UNITE QIIME-release archive."""
import argparse
import os
from multiprocessing import cpu_count
from pathlib import Path

from qiime2_its import downloader, env_checks, qiime_wrapper
from qiime2_its._version import __version__


def _find_unite_files(search_dir):
    """Locate the UNITE "developer" QIIME-release sequence and taxonomy files."""
    seq_file = taxo_file = None
    for root, _dirs, files in os.walk(search_dir):
        for name in files:
            path = Path(root) / name
            path_str = str(path)
            if 'developer' in path_str and 'sh_refs_qiime' in path_str and '_99_' in path_str \
                    and path_str.endswith('.fasta'):
                seq_file = path
            elif 'developer' in path_str and 'sh_taxonomy_qiime' in path_str and '_99_' in path_str \
                    and path_str.endswith('.txt'):
                taxo_file = path
    if seq_file is None or taxo_file is None:
        raise FileNotFoundError('Could not find UNITE "developer" QIIME-release sequence/taxonomy files '
                                 'under {}.'.format(search_dir))
    return seq_file, taxo_file


def fix_fasta(input_fasta, output_fasta):
    """Strip whitespace and upper-case sequence lines so QIIME2 accepts the UNITE fasta."""
    with open(input_fasta) as in_f, open(output_fasta, 'w') as out_f:
        for line in in_f:
            line = line.rstrip().replace(' ', '')
            if not line.startswith('>'):
                line = line.upper()
            out_f.write(f'{line}\n')


def run(url, output_folder, threads):
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    env_checks.check_qiime2_env_active()

    if os.path.isfile(url):
        archive = Path(url)
    else:
        print('Downloading UNITE database...')
        archive = output_folder / 'unite.tar.gz'
        downloader.download(url, archive)

    print('Extracting database...')
    downloader.extract_targz(archive, output_folder)
    unite_seq, unite_taxo = _find_unite_files(output_folder)

    print('Proofreading database...')
    unite_seq_fixed = unite_seq.with_name(unite_seq.stem + '_upper.fasta')
    fix_fasta(unite_seq, unite_seq_fixed)

    qiime2_seq = unite_seq_fixed.with_suffix(unite_seq_fixed.suffix + '.qza')
    qiime2_taxo = unite_taxo.with_suffix(unite_taxo.suffix + '.qza')

    print('Preparing sequences for QIIME2...')
    qiime_wrapper.import_sequences(unite_seq_fixed, qiime2_seq)
    print('Preparing taxonomy for QIIME2...')
    qiime_wrapper.import_taxonomy(unite_taxo, qiime2_taxo)

    name_fields = unite_seq.name.split('_')
    unite_version, unite_clustering, unite_release_date = name_fields[3], name_fields[4], name_fields[5]
    classifier_file = output_folder / f'unite-{unite_version}-{unite_clustering}-classifier-{unite_release_date}.qza'

    print('Training UNITE classifier for QIIME2...')
    qiime_wrapper.train_naive_bayes_classifier(qiime2_seq, qiime2_taxo, classifier_file)

    print(f'DONE. Classifier written to {classifier_file}')


def build_parser():
    cpu = cpu_count()
    parser = argparse.ArgumentParser(description='Train a UNITE classifier for QIIME2')
    parser.add_argument('-u', '--url', metavar='http://www.fileserver.com/file.txt', required=True, type=str,
                         help='URL of the UNITE QIIME2 database (compressed tar file with sequences and '
                              'taxonomy), or the path to an already-downloaded copy. Mandatory.')
    parser.add_argument('-o', '--output_folder', metavar='/unite_folder/', required=True, type=str,
                         help='Output folder for the classifier. Mandatory.')
    parser.add_argument('-q', '--qiime2', metavar='rachis-qiime2-2026.7', required=True, type=str,
                         help='Name of your QIIME2 conda environment. Mandatory.')
    parser.add_argument('-t', '--threads', metavar=str(cpu), default=cpu, type=int,
                         help=f'Number of CPUs. Default is {cpu}.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser


def main():
    args = build_parser().parse_args()
    run(args.url, args.output_folder, args.threads)


if __name__ == '__main__':
    main()
