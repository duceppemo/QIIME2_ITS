#!/usr/local/env python3
import http.client
from argparse import ArgumentParser
from multiprocessing import cpu_count
from time import time
import shutil
from urllib.request import urlopen
from pathlib import Path
import os
import tarfile
from Bio import Entrez
from concurrent import futures
from time import sleep
import gzip
import numpy as np
import subprocess


__author__ = 'duceppemo'
__version__ = '0.1'

from bokeh.io import output_file
from debugpy.common.messaging import MessageDict


# TODO: add download progressbar
#  (https://stackoverflow.com/questions/41106599/python-3-5-urllib-request-urlopen-progress-bar-available)
# TODO: check if taxonomy download files are already present and skip download
# TODO: resume failed taxonomy file download and/or try many times
# TODO: Parallel parse taxonomy files


class Methods(object):

    @staticmethod
    def elapsed_time(seconds):
        """
        Transform a time value into a string
        :param seconds: Elapsed time in seconds
        :return: Formatted time string
        """
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(np.round(value)), name) for name, value in periods if value)
        return time_string

    @staticmethod
    def make_folder(folder):
        """
        Create output folder.
        :param folder: string. Output folder path.
        :return:
        """
        # Will create parent directories if don't exist and will not return error if already exists
        Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def parse_id_table(id_table):
        id_dict = dict()
        with gzip.open(id_table, 'rt') if id_table.endswith('.gz') else open(id_table, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:  # Empty line
                    continue
                field_list = line.split('\t')
                if len(field_list) != 2:
                    raise Exception('"id_list" should have two tab-separated columns')
                else:
                    # taxid: acc
                    id_dict[field_list[1]] = field_list[0]
        return id_dict

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
            with urlopen(url) as response:
                shutil.copyfileobj(response, out_file)

    @staticmethod
    def untargz(targz_file, output_path):
        """
        Decompress and return path of UNITE sequence and taxonomy files
        :param targz_file: string. Path to .tar.gz file
        :param output_path: string. Path to uncompressed file
        :return:
        """

        if targz_file.endswith('.tar.gz'):
            with tarfile.open(targz_file, "r:gz") as f:  # Open for reading with gzip compression
                f.extractall(path=output_path)

    @staticmethod
    def expand_taxonomy_from_taxid(id_dict, taxonomy_file, nodes_file, names_file, merged_file):
        # Parse nodes.dmp file
        node_dict = dict()
        with open(nodes_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                fields = line.split('\t')
                taxid = fields[0]
                parent = fields[2]
                rank = fields[4]

                node_dict[taxid] = (parent, rank)

        # Parse names.dmp file
        names_dict = dict()
        with open(names_file, 'r') as f:
            for line in f:
                if 'scientific name' in line:
                    fields = line.split('\t')
                    taxid = fields[0]
                    name = fields[2]
                    if taxid not in names_dict:
                        names_dict[taxid] = name
                    else:
                        continue

        # Parse meged.dmp
        with open(merged_file, 'r') as f:
            for line in f:
                fields = line.split('\t')
                old_taxid = fields[0]
                new_taxid = fields[2]
                if old_taxid in id_dict:
                    id_dict[new_taxid] = id_dict[old_taxid]
                    del id_dict[old_taxid]

        with open(taxonomy_file, 'w') as f:
            # for taxid in taxid_list:
            for taxid, acc in id_dict.items():
                taxo_list = list()
                while taxid != '1':
                    (parent, rank) = node_dict[taxid]
                    taxo_list.insert(0, (rank, names_dict[taxid]))
                    taxid = parent

                # Write taxonomy to file
                # Kingdom;Phylum;Class;Order;Family;Genus;Species
                # k__;p__;c__;o__;f__;g__;s__
                # unidentified
                taxo_dict = {'k': 'unidentified', 'p': 'unidentified', 'c': 'unidentified',
                             'o': 'unidentified', 'f': 'unidentified', 'g': 'unidentified', 's': 'unidentified'}

                for i in taxo_list:
                    if i[0] == 'species':
                        taxo_dict['s'] = i[1].replace(' ', '_')
                    elif i[0] == 'genus':
                        taxo_dict['g'] = i[1]
                    elif i[0] == 'family':
                        taxo_dict['f'] = i[1]
                    elif i[0] == 'order':
                        taxo_dict['o'] = i[1]
                    elif i[0] == 'clade':
                        taxo_dict['c'] = i[1]
                    elif i[0] == 'phylum':
                        taxo_dict['p'] = i[1]
                    elif i[0] == 'superkingdom':
                        taxo_dict['k'] = i[1]
                    else:
                        continue

                f.write('{}\tk__{};p__{};c__{};o__{};f__{};g__{};s__{}\n'.format(
                    acc, taxo_dict['k'], taxo_dict['p'], taxo_dict['c'],
                    taxo_dict['o'], taxo_dict['f'], taxo_dict['g'], taxo_dict['s']))

    @staticmethod
    def qiime2_import_ncbi_sequences(sequence_file, qiime2_sequence_file):
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
    def qiime2_import_ncbi_taxonomy(taxonomy_file, qiime2_taxonomy_file):
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
    def qiime2_train_ncbi_classifier(sequence_file, taxonomy_file, classifier_file):
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
               '--o-classifier', classifier_file, '--verbose']

        subprocess.run(cmd)


class DbBuilder(object):
    taxdump_url = 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
    acc2taxid_url = 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
    dead_acc2taxid_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz'

    def __init__(self, args):
        self.args = args

        # Args
        self.query = args.query
        self.id_table = args.id_table
        self.output_folder = args.output
        self.cpu = args.threads
        self.email = args.email
        self.api = args.api_key
        self.taxdump = args.taxdump
        self.acc2taxid = args.acc2taxid
        self.dead_acc2taxid = args.dead_acc2taxid

        # Run
        self.run()

    def run(self):
        # Time tracking
        t_zero = time()

        # Check a few things before starting processing data
        self.checks()

        seq_file = self.output_folder + '/seq.fasta'

        # Download and extract taxdump
        start_time = time()
        if self.taxdump:
            print('Extracting taxdump.tar.gz... ', end="", flush=True)
            Methods.untargz(self.taxdump, self.output_folder)
            end_time = time()
        else:
            print('Downloading taxdump.tar.gz... ', end="", flush=True)
            Methods.download(DbBuilder.taxdump_url, self.output_folder + '/taxdump.tar.gz')
            Methods.untargz(self.output_folder + '/taxdump.tar.gz', self.output_folder)
            end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        # Parse id_list
        acc_dict = Methods.parse_id_table(self.id_table)

        # Expand taxonomy from taxid
        print('Writting taxonomy...', end="", flush=True)
        start_time = time()
        taxonomy_file = self.output_folder + '/taxonomy.txt'
        nodes_file = self.output_folder + '/nodes.dmp'
        names_file = self.output_folder + '/names.dmp'
        merged_file = self.output_folder + '/merged.dmp'
        Methods.expand_taxonomy_from_taxid(acc_dict, taxonomy_file, nodes_file, names_file, merged_file)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        # QIIME2
        qiime2_seq = self.query + '.qza'
        qiime2_taxo = taxonomy_file + '.qza'
        classifier_file = self.output_folder + '/naive-bayes_classifier.qza'

        print('Importing fasta file in QIIME2...')
        Methods.qiime2_import_ncbi_sequences(self.query, qiime2_seq)

        print('Importing  taxonomy file in QIIME2...')
        Methods.qiime2_import_ncbi_taxonomy(taxonomy_file, qiime2_taxo)

        print('Training classifier...')
        Methods.qiime2_train_ncbi_classifier(qiime2_seq, qiime2_taxo, classifier_file)

        end_time = time()
        interval = end_time - t_zero
        print("Done (total time: %s)." % Methods.elapsed_time(interval))

    def checks(self):
        # Check query string
        if not os.path.isfile(self.query) or not os.path.exists(self.query):
            raise Exception('Your "query" is not a file or does not exist.')

        # Check id_list
        if not os.path.isfile(self.id_table) or not os.path.exists(self.id_table):
            raise Exception('Your "id-table" is not a file or does not exist.')

        # Check output folder
        if not os.path.exists(self.output_folder):
            Methods.make_folder(self.output_folder)  # Create output folder if it does not exist already

        # Check cpu and parallel processes
        cpu = cpu_count()
        if 1 > self.cpu > cpu:  # smaller than 1 or greater than available cpu
            self.cpu = cpu


if __name__ == '__main__':

    parser = ArgumentParser(description='Prep a QIIME2 classifier from a fasta file.')
    parser.add_argument('-q', '--query', metavar='my_sequences.fasta',
                        required=True,
                        type=str,
                        help='A fasta file. Mandatory.')
    parser.add_argument('-i', '--id-table', metavar='/path/to/acc2taxid_table.tsv',
                        required=True,
                        type=str,
                        help='Tab-separated text file with 2 columns (accession + taxid) matching your input fasta. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        type=str,
                        help='Output folder. Mandatory.')
    parser.add_argument('-t', '--threads', metavar='4',
                        required=False, default=4,
                        type=int,
                        help='Number of CPU. Default is 4. Optional.')
    parser.add_argument('-e', '--email', metavar='your.email@example.org',
                        required=False, default='\'your.email@example.org\'',
                        type=str,
                        help='Your email address. Optional.')
    parser.add_argument('-a', '--api-key', metavar='',
                        required=False,
                        type=str,
                        help='Your NCBI API key. Allows up to 10 requests per second instead of 3. Optional.')
    parser.add_argument('--taxdump', metavar='/path/to/taxdump.tar.gz',
                        required=False,
                        type=str,
                        help='Path to downloaded taxdump.tar.gz. Optional.')
    parser.add_argument('--acc2taxid', metavar='/path/to/nucl_gb.accession2taxid.gz',
                        required=False,
                        type=str,
                        help='Path to downloaded nucl_gb.accession2taxid.gz. Optional.')
    parser.add_argument('--dead-acc2taxid', metavar='/path/to/dead_nucl.accession2taxid.gz',
                        required=False,
                        type=str,
                        help='Path to downloaded dead_nucl.accession2taxid.gz. Optional.')


    # Get the arguments into an object
    arguments = parser.parse_args()

    DbBuilder(arguments)
