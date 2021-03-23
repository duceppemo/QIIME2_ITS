#!/usr/local/env python3

import os
import gzip
import pathlib
from argparse import ArgumentParser
from concurrent import futures
from multiprocessing import cpu_count
from qiime2_methodss import Qiime2Methods


class FastqRC(object):
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        "N": 'N'
    }

    def __init__(self, args):
        # Get options
        self.args = args
        self.input_folder = args.input
        self.output_folder = args.output
        self.cpu = args.threads

        # Dictionary to store all the fastq file(s)
        self.master_dict = dict()

        # Run
        self.run()

    def run(self):
        """
        This is where all the main program parts are being called
        :return:
        """
        # Run all checks to make sure the script runs smoothly
        self.checks()
        # Create a list of all the fastq files in the input folder
        fastq_list = Qiime2Methods.list_fastq(self.input_folder)

        # Parallel process fastq files
        print('Processing:')
        FastqRC.rc_fastq_parallel(fastq_list, self.output_folder, self.cpu)

        # # Process each fastq file
        # print('Processing:')
        # for fq in fastq_list:
        #     FastqRC.rc_fastq(fq, self.output_folder)

    def checks(self):
        """
        Check if all options have been filled properly
        :return:
        """
        # Check input folder
        if not os.path.exists(self.input_folder):
            raise Exception('Input folder does not exists')
        if not os.path.isdir(self.input_folder):
            raise Exception('Input is not a folder')

        # Check output folder
        if self.output_folder == self.input_folder:
            raise Exception('Please choose a different output folder than the input folder')
        if not os.path.exists(self.output_folder):
            # Create output folder is does not exists already
            Qiime2Methods.make_folder() .make_folder(self.output_folder)

    @staticmethod
    def rc_fastq_parallel(fastq_list, out_folder, ncpu):
        """
        Parallel process the fastq files
        :param fastq_list: string. List of fastq file paths
        :param out_folder: string. Output folder path
        :param ncpu: int. Number of files to process in parallel
        :return:
        """
        with futures.ThreadPoolExecutor(max_workers=ncpu) as executor:
            args = ((fastq, out_folder) for fastq in fastq_list)
            for results in executor.map(lambda p: FastqRC.rc_fastq(*p), args):  # (*p) does the unpacking part
                pass

    @staticmethod
    def rc_fastq(fastq, out_folder):
        """
        Parse fastq file into a dictionary
        :param fastq: string. Input fastq file path
        :param out_folder: string. Output folder path
        :return:
        """
        print('\t{}'.format(os.path.basename(fastq)))

        # Create file path of output file
        outfile = out_folder + '/' + os.path.basename(fastq)
        # Since all output files are gzipped, make sure to add ".gz" extensions to the input files
        # that were not compressed
        if not outfile.endswith('.gz'):
            outfile += '.gz'

        with gzip.open(outfile, 'wb') as out_f:
            with gzip.open(fastq, 'rb') if fastq.endswith('.gz') else open(fastq, 'r') as in_f:
                counter = 0
                for line in in_f:
                    counter += 1
                    if fastq.endswith('.gz'):
                        line = line.decode('ascii')  # Convert to ASCII if files were gzipped (binary)
                    line = line.rstrip()

                    if counter == 1:  # Header line line
                        if not line.startswith('@'):  # first line of 4 should always start with "@"
                            raise Exception('Fastq file invalid.')
                        out_f.write('{}\n'.format(line).encode())
                    elif counter == 2:  # Sequence line
                        line = line[::-1]  # reverse
                        seq_rc = ''
                        for bp in line:
                            seq_rc += FastqRC.complement[bp]
                        out_f.write('{}\n'.format(seq_rc).encode())
                    elif counter == 3:  # that's the "+" line, we can skip it
                        out_f.write('{}\n'.format(line).encode())
                    elif counter == 4:  # Quality line
                        out_f.write('{}\n'.format(line[::-1]).encode())  # reverse
                        counter = 0  # Reset counter because each fastq entry in a file always have 4 lines


if __name__ == '__main__':
    cpu = cpu_count()

    parser = ArgumentParser(description='Reverse complement all entries in a fastq file.')
    parser.add_argument('-i', '--input', metavar='/input_folder/',
                        required=True,
                        type=str,
                        help='Input folder with fastq file(s), gzipped or not. '
                             'Accepted file extensions are ".fastq", ".fastq.gz", ".fq" and ".fq.gz". '
                             'Will look for files recursively. '
                             'Mandatory.')
    parser.add_argument('-o', '--output', metavar='/modified_fastq/',
                        required=True,
                        type=str,
                        help='Output folder. Must be different than input folder. '
                             'Mandatory.')
    parser.add_argument('-t', '--threads', metavar='{}'.format(cpu),
                        required=False, default=cpu,
                        type=int,
                        help='Number of CPU. Default is {}'.format(cpu))

    # Get the arguments into an object
    arguments = parser.parse_args()

    FastqRC(arguments)
