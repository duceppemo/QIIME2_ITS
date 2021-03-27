#!/usr/local/env python3

import os
import subprocess
from argparse import ArgumentParser
from multiprocessing import cpu_count
from qiime2_methods import Qiime2Methods


class Qiime2(object):
    def __init__(self, args):
        self.args = args
        self.input_folder = args.input
        self.output_folder = args.output
        self.qiime2_env = args.qiime2
        self.cpu = args.threads
        self.unite_qiime2_file = args.unite
        self.metadata_file = args.metadata
        self.revers_complement = args.reverse_complement
        self.min_len = args.min_len
        self.max_len = args.max_len

        # Data
        # Create a list of all the fastq files in the input folder
        self.fastq_list = list()

        # Run
        self.run()

    def run(self):
        # List fastq files from input folder
        self.fastq_list = Qiime2Methods.list_fastq(self.input_folder)

        # Check a few things before starting processing data
        self.checks()

        # Check for reverse complement flag
        if self.revers_complement:
            print('Reverse complementing reads...')
            rc_folder = self.output_folder + '/rc_reads'
            Qiime2Methods.make_folder(rc_folder)
            Qiime2.run_fastq_rc(self.input_folder, rc_folder)
            self.input_folder = rc_folder

        # Create ITSxpress output folder and log folder
        its_out = self.output_folder + '/ITSxpress'
        Qiime2Methods.make_folder(its_out + '/log')
        # Extract Fungi ITS1 in parallel
        Qiime2Methods.extract_its_parallel(self.fastq_list, its_out, self.cpu)

        # Remove empty sequences


    def checks(self):
        # Check input folder
        if not self.fastq_list:
            raise Exception('No fastq files found in provided input folder.')

        # Check output folder
        if not os.path.exists(self.output_folder):
            Qiime2Methods.make_folder(self.output_folder)  # Create output folder is does not exists already

        # Check if qiime2 conda environment is activated
        if 'qiime2' not in os.environ['CONDA_DEFAULT_ENV']:
            raise Exception('You must activate your QIIME2 conda environment to run this script. '
                            '"conda activate qiime2-2020.8"')

        # Check file names
        # Valid QIIME2 file names:
        # L2S357_15_L001_R1_001.fastq.gz. The underscore-separated fields in this file name are:
        #     1. the sample identifier,
        #     2. the barcode sequence or a barcode identifier,
        #     3. the lane number,
        #     4. the direction of the read (i.e. only R1, because these are single-end reads), and
        #     5. the set number.

        # Check metadata
        # Very basic check. Only checking that sample names are the same as fastq files

    @staticmethod
    def run_fastq_rc(input_folder, output_folder):
        cmd = ['python3', 'fastq_rc.py',
               '-i', input_folder,
               '-o', output_folder]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_validate_metadata(metadata_file):
        """

        :param metadata_file:
        :return:
        """
        import csv
        with open(metadata_file) as f:
            tsv_reader = csv.reader(f, delimiter="\t")
            for line in tsv_reader:
                pass


if __name__ == '__main__':
    cpu = cpu_count()

    parser = ArgumentParser(description='Run QIIME2 on IonTorrent sequencing data using the UNITE database')
    parser.add_argument('-q', '--qiime2', metavar='qiime2-2020.8',
                        required=True,
                        type=str,
                        help='Name of your QIIME2 conda environment.'
                             'Mandatory.')
    parser.add_argument('-i', '--input', metavar='/input_folder/',
                        required=True,
                        type=str,
                        help='Input folder where the fastq reads are located.'
                             'Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        type=str,
                        help='Output folder for QIIME2 files.'
                             'Mandatory.')
    parser.add_argument('-m', '--metadata', metavar='qiime2_metadata.tsv',
                        required=True,
                        type=str,
                        help='Validated QIIME2 metadata file (samples description).'
                             'Mandatory.')
    parser.add_argument('-u', '--unite', metavar='unite_classifier_qiime2.qza',
                        required=True,
                        type=str,
                        help='UNITE classifier for QIIME2. See script "train_unite_classifier_qiime2.py" to compile it.'
                             'Mandatory.')
    parser.add_argument('-t', '--threads', metavar='{}'.format(cpu),
                        required=False, default=cpu,
                        type=int,
                        help='Number of CPU. Default is {}'.format(cpu))
    parser.add_argument('-rc', '--reverse_complement',
                        action='store_true',  # 'store_true' means False by default
                        required=False,
                        help='Use this flag is your reads are in reverse complement. '
                             'For example if you sequenced from 5.8S to 18S.')
    parser.add_argument('--min_len',
                        type=int, default=0,
                        required=False,
                        help='Minimum read length to keep. Default is 0 (no min length).')
    parser.add_argument('--max_len',
                        type=int, default=0,
                        required=False,
                        help='Maximum read length to keep. Default is 0 (no max length).')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Qiime2(arguments)
