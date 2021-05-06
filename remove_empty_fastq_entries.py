#!/usr/local/env python3

import os
import gzip
from argparse import ArgumentParser
from collections import OrderedDict


class RemoveEmpties(object):
    extensions = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']

    def __init__(self, args):
        # Arguments
        self.args = args
        self.fastq = args.fastq
        self.fastq2 = args.fastq2

        # Run
        self.run()

    def run(self):
        self.checks()
        if not self.fastq2:
            RemoveEmpties.remove_empties_se(self.fastq)
        else:
            RemoveEmpties.remove_empties_pe(self.fastq, self.fastq2)

    def checks(self):
        # Check if input fastq(s) exist(s)
        if not os.path.exists(self.fastq):
            raise Exception('Input fastq doest not exists.')

        # Check if input fastq is a file
        if not os.path.isfile(self.fastq):
            raise Exception('You must provide a file as input.')

        # Check if input fastq has the right file extension
        if not self.fastq.endswith(tuple(RemoveEmpties.extensions)):
            raise Exception('Accepted file extensions are ".fastq", ".fastq.gz", ".fq" and ".fq.gz".')

    @staticmethod
    def remove_empties_se(input_fastq):
        out_fastq = input_fastq + '.clean'

        with gzip.open(out_fastq, 'wb') as out_f:
            with gzip.open(input_fastq, 'rb') if input_fastq.endswith('gz') else open(input_fastq, 'w') as in_f:
                seq_list = list()
                counter = 0
                seq_counter = 0
                good_counter = 0
                empties_counter = 0
                for line in in_f:
                    counter += 1
                    line = line.rstrip()
                    seq_list.append(line)

                    if counter == 4:
                        seq_counter += 1
                        counter = 0
                        if seq_list[1] == b'':
                            # print('Sequence {} is empty'.format(seq_list[0].decode()))
                            seq_list = list()
                            empties_counter += 1
                        else:
                            out_f.write(b'\n'.join(seq_list) + b'\n')
                            seq_list = list()
                            good_counter += 1
                print('{}: total sequences: {}, good sequences: {}, empty sequences: {}'.format(
                    os.path.basename(input_fastq), seq_counter, good_counter, empties_counter))

        # Replace old fastq file with new one
        os.replace(out_fastq, input_fastq)

    @staticmethod
    def remove_empties_pe(r1, r2):
        out_r1 = r1 + '.clean'
        out_r2 = r2 + '.clean'

        seq_dict = OrderedDict()

        # Parse r1 into dictionary
        with gzip.open(r1, 'rb') if r1.endswith('gz') else open(r1, 'w') as in_f:
            seq_list = list()
            counter = 0
            seq_counter = 0

            for line in in_f:
                counter += 1
                line = line.rstrip()
                seq_list.append(line)

                if counter == 4:
                    seq_dict[seq_counter] = [seq_list]
                    seq_counter += 1
                    counter = 0
                    seq_list = list()

        # Add r2 to dictionary
        with gzip.open(r2, 'rb') if r2.endswith('gz') else open(r2, 'w') as in_f:
            seq_list = list()
            counter = 0
            seq_counter = 0

            for line in in_f:
                counter += 1
                line = line.rstrip()
                seq_list.append(line)

                if counter == 4:
                    seq_dict[seq_counter].append(seq_list)
                    seq_counter += 1
                    counter = 0
                    seq_list = list()

        with gzip.open(out_r1, 'wb') as fh_r1:
            with gzip.open(out_r2, 'wb') as fh_r2:
                empties_counter = 0
                for read, reads_list in seq_dict.items():
                    if seq_dict[read][0][1] == b'' or seq_dict[read][1][1] == b'':
                        empties_counter += 1
                        continue
                    else:
                        fh_r1.write(b'\n'.join(seq_dict[read][0]) + b'\n')
                        fh_r2.write(b'\n'.join(seq_dict[read][1]) + b'\n')

                print('{}: input sequences: {}, good: {}, empty: {}'.format(
                    os.path.basename(r1).split('_')[0], len(seq_dict.keys()),
                    len(seq_dict.keys()) - empties_counter, empties_counter))

        # Replace old fastq files with new one
        os.replace(out_r1, r1)
        os.replace(out_r2, r2)


if __name__ == '__main__':
    parser = ArgumentParser(description='Remove empty entries in a fastq file.')
    parser.add_argument('-f', '--fastq', metavar='sample_R1.fastq.gz',
                        required=True,
                        type=str,
                        help='Input (R1) fastq file, gzipped or not. '
                             'Accepted file extensions are ".fastq", ".fastq.gz", ".fq" and ".fq.gz". '
                             'Mandatory.')
    parser.add_argument('-f2', '--fastq2', metavar='sample_R1.fastq.gz',
                        required=False,
                        type=str,
                        help='Input R2 fastq file, gzipped or not. '
                             'Accepted file extensions are ".fastq", ".fastq.gz", ".fq" and ".fq.gz". '
                             'Mandatory.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    RemoveEmpties(arguments)
