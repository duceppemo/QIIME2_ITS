"""Train a QIIME2 naive-Bayes classifier from a UNITE QIIME-release archive."""
import argparse
import os
from pathlib import Path

from qiime2_its import downloader, env_checks, qiime_wrapper
from qiime2_its._version import __version__


def _find_unite_files(search_dir, clustering='99'):
    """Locate the UNITE "developer" QIIME-release sequence and taxonomy files
    for the given percent-similarity clustering threshold ('97', '99', or
    'dynamic').

    Real UNITE archives also ship "_all_" (singleton-inclusive) variants and
    other clustering thresholds side by side, all sharing the substrings a
    looser match would key on -- e.g. both "sh_refs_qiime_ver10_99_<date>.fasta"
    and "sh_refs_qiime_ver10_99_all_<date>.fasta" contain "sh_refs_qiime" and
    "_99_". Files are matched on their exact field layout instead, since the
    "_all_" variant has one extra "_"-separated field.
    """
    seq_file = taxo_file = None
    for root, _dirs, files in os.walk(search_dir):
        if 'developer' not in Path(root).parts:
            continue
        for name in files:
            path = Path(root) / name
            fields = path.stem.split('_')
            if len(fields) != 6 or fields[4] != clustering:
                continue
            if fields[:3] == ['sh', 'refs', 'qiime'] and name.endswith('.fasta'):
                seq_file = path
            elif fields[:3] == ['sh', 'taxonomy', 'qiime'] and name.endswith('.txt'):
                taxo_file = path
    if seq_file is None or taxo_file is None:
        raise FileNotFoundError('Could not find UNITE "developer" QIIME-release sequence/taxonomy files '
                                 'under {} for clustering "{}".'.format(search_dir, clustering))
    return seq_file, taxo_file


def parse_unite_filename(seq_filename):
    """Parse "sh_refs_qiime_<version>_<clustering>_<release_date>.fasta" into its
    (version, clustering, release_date) fields, extension-safe.
    """
    fields = Path(seq_filename).stem.split('_')
    return fields[3], fields[4], fields[5]


def fix_fasta(input_fasta, output_fasta):
    """Strip whitespace and upper-case sequence lines so QIIME2 accepts the UNITE fasta."""
    with open(input_fasta) as in_f, open(output_fasta, 'w') as out_f:
        for line in in_f:
            line = line.rstrip().replace(' ', '')
            if not line.startswith('>'):
                line = line.upper()
            out_f.write(f'{line}\n')


def run(url, output_folder):
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

    unite_version, unite_clustering, unite_release_date = parse_unite_filename(unite_seq.name)
    classifier_file = output_folder / f'unite-{unite_version}-{unite_clustering}-classifier-{unite_release_date}.qza'

    print('Training UNITE classifier for QIIME2...')
    qiime_wrapper.train_naive_bayes_classifier(qiime2_seq, qiime2_taxo, classifier_file)

    print(f'DONE. Classifier written to {classifier_file}')


def build_parser():
    parser = argparse.ArgumentParser(description='Train a UNITE classifier for QIIME2')
    parser.add_argument('-u', '--url', metavar='http://www.fileserver.com/file.txt', required=True, type=str,
                         help='URL of the UNITE QIIME2 database (compressed tar file with sequences and '
                              'taxonomy), or the path to an already-downloaded copy. Mandatory.')
    parser.add_argument('-o', '--output_folder', metavar='/unite_folder/', required=True, type=str,
                         help='Output folder for the classifier. Mandatory.')
    parser.add_argument('-q', '--qiime2', metavar='rachis-qiime2-2026.7', required=True, type=str,
                         help='Name of your QIIME2 conda environment. Mandatory.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser


def main():
    args = build_parser().parse_args()
    run(args.url, args.output_folder)


if __name__ == '__main__':
    main()
