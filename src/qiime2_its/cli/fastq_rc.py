"""Standalone utility: reverse-complement every fastq file in a folder."""
import argparse
from multiprocessing import cpu_count
from pathlib import Path

from qiime2_its import env_checks, fastq_utils
from qiime2_its._version import __version__


def run(input_folder, output_folder, threads):
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)

    if not input_folder.exists() or not input_folder.is_dir():
        raise ValueError('Input folder does not exist or is not a directory.')
    # Not inside it either: the fastq search below is recursive, so a second
    # run would find (and try to re-process) the first run's own output.
    if output_folder.resolve() == input_folder.resolve() or input_folder.resolve() in output_folder.resolve().parents:
        raise ValueError('Please choose an output folder different from (and not inside) the input folder.')
    output_folder.mkdir(parents=True, exist_ok=True)

    fastq_list = fastq_utils.list_fastq(input_folder)
    if not fastq_list:
        raise ValueError('No fastq files found in the input folder.')

    print('Processing:')
    fastq_utils.rc_fastq_parallel(fastq_list, output_folder, env_checks.clamp_cpu(threads))


def build_parser():
    cpu = cpu_count()
    parser = argparse.ArgumentParser(description='Reverse complement all entries in fastq file(s).')
    parser.add_argument('-i', '--input', metavar='/input_folder/', required=True, type=str,
                         help='Input folder with fastq file(s), gzipped or not. Accepted extensions are '
                              '".fastq", ".fastq.gz", ".fq" and ".fq.gz". Searched recursively. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/modified_fastq/', required=True, type=str,
                         help='Output folder. Must be different from (and not inside) the input folder. Mandatory.')
    parser.add_argument('-t', '--threads', metavar=str(cpu), default=cpu, type=int,
                         help=f'Number of CPUs. Default is {cpu}.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser


def main():
    args = build_parser().parse_args()
    run(args.input, args.output, args.threads)


if __name__ == '__main__':
    main()
