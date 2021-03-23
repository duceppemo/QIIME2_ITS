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

        # Run
        self.run()

    def run(self):
        self.checks()

        # Create a list of all the fastq files in the input folder
        fastq_list = Qiime2Methods.list_fastq(self.input_folder)

        # Check for reverse complement flag
        if self.revers_complement:
            rc_folder = self.output_folder + '/rc_reads'
            Qiime2Methods.make_folder(rc_folder)
            Qiime2.run_fastq_rc(self.input_folder, rc_folder)
            self.input_folder = rc_folder

        #

    def checks(self):
        # Check output folder
        if not os.path.exists(self.output_folder):
            Qiime2Methods.make_folder(self.output_folder)  # Create output folder is does not exists already

        # Check if qiime2 conda environment is activated
        if 'qiime2' not in os.environ['CONDA_DEFAULT_ENV']:
            raise Exception('You must activate your QIIME2 conda environment to run this script. '
                            '"conda activate qiime2-2020.8"')

    @staticmethod
    def run_fastq_rc(input_folder, output_folder):
        cmd = ['python3', 'fastq_rc.py',
               '-i', input_folder,
               '-o', output_folder]
        subprocess.run(cmd)


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

    # Get the arguments into an object
    arguments = parser.parse_args()

    Qiime2(arguments)
