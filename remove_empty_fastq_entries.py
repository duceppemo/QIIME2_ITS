#!/usr/local/env python3

import os
import gzip
from argparse import ArgumentParser


class RemoveEmpties(object):
    extensions = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']

    def __init__(self, args):
        self.args = args
        self.input_fastq = args.input_fastq

        # Run
        self.run()

    def run(self):

        RemoveEmpties.remove_empties(self.input_fastq)

    def checks(self):
        # Check if input fastq exists
        if not os.path.exists(self.input_fastq):
            raise Exception('Input fastq doest not exists.')

        # Check if input fastq is a file
        if not os.path.isfile(self.input_fastq):
            raise Exception('You must provide a file as input.')

        # Check if input fastq has the right file extension
        if not self.input_fastq.endswith(tuple(RemoveEmpties.extensions)):
            raise Exception('Accepted file extensions are ".fastq", ".fastq.gz", ".fq" and ".fq.gz".')

    @staticmethod
    def remove_empties(input_fastq):
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
                        if seq_list[1] == '':
                            print('Sequence {} is empty'.format(seq_list[0]))
                            seq_list = list()
                            empties_counter += 1
                        else:
                            out_f.write('\n'.join(seq_list) + '\n')
                            seq_list = list()
                            good_counter += 1
                print('{}: total sequences: {}, good sequences: {}, empty sequences: {}'.format(
                    os.path.basename(input_fastq), seq_counter, good_counter, empties_counter))

        # Replace old fastq file with new one
        os.replace(out_fastq, input_fastq)


if __name__ == '__main__':
    parser = ArgumentParser(description='Remove empty entries in a fastq file.')
    parser.add_argument('-i', '--input_fastq', metavar='sample.fastq.gz',
                        required=True,
                        type=str,
                        help='Input fastq file, gzipped or not. '
                             'Accepted file extensions are ".fastq", ".fastq.gz", ".fq" and ".fq.gz". '
                             'Mandatory.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    RemoveEmpties(arguments)
