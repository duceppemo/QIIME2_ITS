#!/usr/local/env python3

import os
import subprocess
import shutil
import urllib.request
import tarfile
from argparse import ArgumentParser
from multiprocessing import cpu_count
from qiime2_methods import Qiime2Methods


__author__ = 'duceppemo'
__version__ = '0.1'


class Qiime2Trainer(object):
    def __init__(self, args):
        self.args = args
        self.url = args.url
        self.output_folder = args.output_folder
        self.qiime2_env = args.qiime2
        self.cpu = args.threads

        # Run
        self.run()

    def run(self):
        self.checks()
        tmp_file = self.output_folder + '/tmp.tar.gz'

        if os.path.exists(self.url) and os.path.isfile(self.url):
            tmp_file = self.url
        else:
            print('Downloading UNITE database')
            Qiime2Trainer.download(self.url, tmp_file)

        # Uncompress DB
        print('Extracting database...')
        (unite_seq, unite_taxo) = Qiime2Trainer.untargz(tmp_file)  # tuple
        unite_seq_fixed = '.'.join(unite_seq.split('.')[:-1]) + '_upper.fasta'
        print('Proofreading database...')
        Qiime2Trainer.fix_fasta(unite_seq, unite_seq_fixed)

        qiime2_seq = unite_seq_fixed + '.qza'
        qiime2_taxo = unite_taxo + '.qza'

        print('Preparing sequences for QIIME2...')
        Qiime2Trainer.qiime2_import_unite_sequences(unite_seq_fixed, qiime2_seq)
        print('Preparing taxonomy for QIIME2...')
        Qiime2Trainer.qiime2_import_unite_taxonomy(unite_taxo, qiime2_taxo)

        unite_version = os.path.basename(unite_seq).split("_")[3]
        unite_clustering = os.path.basename(unite_seq).split("_")[4]
        unite_release_date = os.path.basename(unite_seq).split("_")[5]
        qiime2_classifier = '{}/unite-{}-{}-classifier-{}.qza'.format(self.output_folder, unite_version,
                                                                      unite_clustering, unite_release_date)

        print('Training UNITE classifier for QIIME2...')
        Qiime2Trainer.qiime2_train_unite_classifier(qiime2_seq, qiime2_taxo, qiime2_classifier)

        # Cleanup temporary files
        print('Removing temporary files...')
        pass

        print('DONE')

    def checks(self):
        # Check output folder
        if not os.path.exists(self.output_folder):
            Qiime2Methods.make_folder(self.output_folder)  # Create output folder is does not exists already

        # Check if qiime2 conda environment is activated
        if 'qiime2' not in os.environ['CONDA_DEFAULT_ENV']:
            raise Exception('You must activate your QIIME2 conda environment to run this script. '
                            '"conda activate qiime2-2020.8"')

    @staticmethod
    def download(url, file_path):
        """
        Download the file from `url` and save it locally under `file_name`
        https://stackoverflow.com/questions/7243750/download-file-from-web-in-python-3
        https://unite.ut.ee/repository.php
        :param url: string. URL of file (remotely)
        :param file_path: string. Path of downloaded file (locally)
        :return:
        """
        with open(file_path, 'wb') as out_file:
            with urllib.request.urlopen(url) as response:
                shutil.copyfileobj(response, out_file)

    @staticmethod
    def untargz(targz_file):
        """
        Decompress and return path of UNITE sequence and taxonomy files
        :param targz_file:
        :return: sting tuple. Path of sequence and taxonomy files
        """
        output_path = os.path.dirname(targz_file)

        if targz_file.endswith(tuple(['.tar.gz', '.tgz'])):
            with tarfile.open(targz_file, "r:gz") as f:  # Open for reading with gzip compression
                def is_within_directory(directory, target):
                    
                    abs_directory = os.path.abspath(directory)
                    abs_target = os.path.abspath(target)
                
                    prefix = os.path.commonprefix([abs_directory, abs_target])
                    
                    return prefix == abs_directory
                
                def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                
                    for member in tar.getmembers():
                        member_path = os.path.join(path, member.name)
                        if not is_within_directory(path, member_path):
                            raise Exception("Attempted Path Traversal in Tar File")
                
                    tar.extractall(path, members, numeric_owner=numeric_owner) 
                    
                
                safe_extract(f, path=output_path)

        # Retrieve path of created directory
        seq_file = ''
        taxo_file = ''
        for (root, dirs, files) in os.walk(output_path):
            for name in files:
                file_path = os.path.join(root, name)
                if ('developer' in file_path) and ('sh_refs_qiime' in file_path) \
                        and ('_99_' in file_path) and ('.fasta' in file_path):
                    seq_file = file_path
                elif ('developer' in file_path) and ('sh_taxonomy_qiime' in file_path) \
                        and ('_99_' in file_path) and ('.txt' in file_path):
                    taxo_file = file_path
        return seq_file, taxo_file

    @staticmethod
    def fix_fasta(input_fasta, output_fasta):
        """
        Fix formatting errors that prevent importation of the reference sequences into QIIME2.
        There are white spaces that interfere, and possibly some lower case letters that need
        to be converted to upper case
        :param input_fasta: sting. Path to input fasta file
        :param output_fasta: sting. Path to modified fasta file
        :return:
        """
        with open(output_fasta, 'w') as out_f:
            with open(input_fasta, 'r') as in_f:
                for line in in_f:
                    line = line.rstrip().replace(' ', '')  # Remove trailing carriage return
                    if not line.startswith('>'):
                        line = line.upper()  # To upper case and remove space
                    out_f.write('{}\n'.format(line))

    @staticmethod
    def qiime2_import_unite_sequences(sequence_file, qiime2_sequence_file):
        """

        :param sequence_file:
        :param qiime2_sequence_file:
        :return:
        """
        cmd = ['qiime', 'tools', 'import',
               '--type', 'FeatureData[Sequence]',
               '--input-path', sequence_file,
               '--output-path', qiime2_sequence_file]

        subprocess.run(cmd)

    @staticmethod
    def qiime2_import_unite_taxonomy(taxonomy_file, qiime2_taxonomy_file):
        """

        :param taxonomy_file:
        :param qiime2_taxonomy_file:
        :return:
        """
        cmd = ['qiime', 'tools', 'import',
               '--type', 'FeatureData[Taxonomy]',
               '--input-format', 'HeaderlessTSVTaxonomyFormat',
               '--input-path', taxonomy_file,
               '--output-path', qiime2_taxonomy_file]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_train_unite_classifier(sequence_file, taxonomy_file, classifier_file):
        """
        # Training the QIIME2 Classifier with UNITE ITS Reference Sequences
        # http://john-quensen.com/tutorials/training-the-qiime2-classifier-with-unite-its-reference-sequences/
        :param sequence_file: string. QIIME2 sequence file (.qza)
        :param taxonomy_file: string. QIIME2 taxonomy file (.qza)
        :param classifier_file: string. QIIME2 classifier file (.qza)
        :return:
        """
        cmd = ['qiime', 'feature-classifier', 'fit-classifier-naive-bayes',
               '--i-reference-reads', sequence_file,
               '--i-reference-taxonomy', taxonomy_file,
               '--o-classifier', classifier_file]

        subprocess.run(cmd)


if __name__ == '__main__':
    cpu = cpu_count()

    parser = ArgumentParser(description='Train Unite classifier for QIIME2')
    parser.add_argument('-u', '--url', metavar='http://www.fileserver.com/file.txt',
                        required=True,
                        type=str,
                        help='URL of UNITE QIIME2 database (compressed tar file with sequences and taxonomy) OR the path the '
                             'the file already downloaded. Mandatory.')
    parser.add_argument('-o', '--output_folder', metavar='/unite_folder/',
                        required=True,
                        type=str,
                        help='Output folder for classifier. Mandatory.')
    parser.add_argument('-q', '--qiime2', metavar='qiime2-2020.8',
                        required=True,
                        type=str,
                        help='Name of your QIIME2 conda environment. Mandatory.')
    parser.add_argument('-t', '--threads', metavar='{}'.format(cpu),
                        required=False, default=cpu,
                        type=int,
                        help='Number of CPU. Default is {}'.format(cpu))

    # Get the arguments into an object
    arguments = parser.parse_args()

    Qiime2Trainer(arguments)
