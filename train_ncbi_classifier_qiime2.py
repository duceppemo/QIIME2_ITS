#!/usr/local/env python3
import http.client
from argparse import ArgumentParser
from multiprocessing import cpu_count
from time import time

import shutil
from itertools import islice
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
import multiprocessing as mp


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
                # Methods.copyfile(response, out_file)

    @staticmethod
    def download_parallel(url_list, path_list, cpu):
        with futures.ThreadPoolExecutor(max_workers=cpu) as executor:
            args = ((u, p) for u in url_list for p in path_list)
            for results in executor.map(lambda p: Methods.download(*p), args):  # (*p) unpacks arguments
                pass

    @staticmethod
    def untargz(targz_file):
        """
        Decompress and return path of UNITE sequence and taxonomy files
        :param targz_file: string. Path to .tar.gz file
        :return:
        """
        # Path of uncompressed file is same as compressed file
        output_path = os.path.dirname(targz_file)

        if targz_file.endswith('.tar.gz'):
            with tarfile.open(targz_file, "r:gz") as f:  # Open for reading with gzip compression
                f.extractall(path=output_path)

    @staticmethod
    def download_seq_from_query(query, seq_file, email, api):
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec167

        Entrez.email = email
        Entrez.api_key = api
        Entrez.sleep_between_tries = 15
        Entrez.max_tries = 3

        # Search
        print('Searching NCBI for \'{}\'.'.format(query))
        search_handle = Entrez.esearch(db='nucleotide', term=query, idtype="acc", usehistory='y')
        search_results = Entrez.read(search_handle)
        search_handle.close()

        count = int(search_results["Count"])
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        print('Found {} records.'.format(count))

        # Fetch
        print('Downloading sequences.')
        batch_size = 300
        with open(seq_file, "w") as out_handle:
            for start in range(0, count, batch_size):
                end = min(count, start + batch_size)
                print("Downloading record %i to %i" % (start + 1, end))
                try:  # Trying to deal with "http.client.IncompleteRead
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        rettype="fasta",
                        retmode="text",
                        retstart=start,
                        retmax=batch_size,
                        webenv=webenv,
                        query_key=query_key,
                        idtype="acc")
                except (http.client.IncompleteRead, ValueError) as e:
                    print('Network error ({}). Last attempt.'.format(e))
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        rettype="fasta",
                        retmode="text",
                        retstart=start,
                        retmax=batch_size,
                        webenv=webenv,
                        query_key=query_key,
                        idtype="acc")

                data = fetch_handle.read()
                fetch_handle.close()
                out_handle.write(data)
                sleep(0.5)

    @staticmethod
    def extract_acc_from_fasta(input_fasta, acc_file):
        acc_dict = dict()
        with open(acc_file, 'w') as out_fh:
            with open(input_fasta, 'r') as in_fh:
                for line in in_fh:
                    line = line.rstrip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        acc = line.split()[0].replace('>', '')
                        out_fh.write('{}\n'.format(acc))
                        acc_dict[acc] = ''

        return acc_dict

    @staticmethod
    def parse_acc2taxid_file(acc2taxid):
        acc2taxid_dict = dict()
        # Parse acc2taxid gzippped file
        with gzip.open(acc2taxid, 'rb') as f:
            f.readline()  # Skip header
            for line in f:
                field_list = line.decode().split('\t')
                acc = field_list[0]
                taxid = field_list[2]
                acc2taxid_dict[acc] = taxid

        return acc2taxid_dict

    # @staticmethod
    # def parse_acc2taxid_chunk(chunk):
    #     acc2taxid_dict = dict()
    #     for line in chunk:
    #         field_list = line.decode().split('\t')
    #         acc = field_list[0]
    #         taxid = field_list[2]
    #         acc2taxid_dict[acc] = taxid
    #
    #     return acc2taxid_dict
    #
    # @staticmethod
    # def make_chunks(file_handle, size):
    #     file_handle.readline()  # skip first line (header)
    #     while True:
    #         chunk = list(islice(file_handle, size))
    #         if not chunk:
    #             break
    #         yield chunk
    #
    # @staticmethod
    # def parse_acc2taxid_file_parallel(acc2taxid, cpu):
    #     # Chunk taxonomy file and run chunks in parallel
    #     acc2taxid_dict = dict()
    #     with gzip.open(acc2taxid, 'rb') as f:
    #         pool = mp.Pool(cpu)
    #         jobs = [pool.apply_async(Methods.parse_acc2taxid_chunk, args=(chunk,))
    #                 for chunk in Methods.make_chunks(f, 1000000)]
    #         results = [j.get() for j in jobs]
    #         pool.close()
    #         pool.join()
    #         pool.terminate()  # Needed to do proper garbage collection?
    #
    #         # Update self.sample_dict with results from every chunk
    #         for d in results:
    #             acc2taxid_dict.update(d)  # Do the merge
    #
    #     return acc2taxid_dict

    @staticmethod
    def acc_to_taxid(acc_dict, acc2taxid_dict, output_file):
        # Parse acc2taxid gzipped file
        missing_acc2taxid = list()
        with open(output_file, 'w') as f:
            for acc, x in acc_dict.items():
                try:
                    f.write('{}\t{}\n'.format(acc, acc2taxid_dict[acc.split('.')[0]]))
                except KeyError as e:
                    # print('The following accession was not found in the taxdump files: {}'.format(e))
                    missing_acc2taxid.append(e.args[0])
                    f.write('{}\t{}\n'.format(acc, '12908'))  # unclassified sequence
        print('The following accession was not found in the taxdump files: {}'.format(', '.join(missing_acc2taxid)))

    @staticmethod
    def expand_taxonomy_from_taxid(taxid_file, taxonomy_file, nodes_file, names_file, merged_file):
        # Parse taxid_file
        taxid_dict = dict()
        with open(taxid_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                acc, taxid = line.split('\t')
                taxid_dict[taxid] = acc

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
                if old_taxid in taxid_dict:
                    taxid_dict[new_taxid] = taxid_dict[old_taxid]
                    del taxid_dict[old_taxid]

        with open(taxonomy_file, 'w') as f:
            # for taxid in taxid_list:
            for taxid, acc in taxid_dict.items():
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


class DbBuilder(object):
    taxdump_url = 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
    acc2taxid_url = 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
    dead_acc2taxid_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz'

    def __init__(self, args):
        self.args = args

        # Args
        self.query = args.query
        self.output_folder = args.output
        self.cpu = args.threads
        self.email = args.email
        self.api = args.api_key

        # Run
        self.run()

    def run(self):
        # time tracking
        t_zero = time()

        # Check a few things before starting processing data
        self.checks()

        # Download and extract taxdump
        print('Downloading taxdump.tar.gz... ', end="", flush=True)
        start_time = time()
        # Methods.download(DbBuilder.taxdump_url, self.output_folder + '/taxdump.tar.gz')
        Methods.untargz(self.output_folder + '/taxdump.tar.gz')
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        # Download nucleotide accession2taxid
        print('Downloading accession2taxid.gz... ', end="", flush=True)
        start_time = time()
        # Methods.download(DbBuilder.acc2taxid_url, self.output_folder + '/nucl_gb.accession2taxid.gz')
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        # Download dead nucleotide accession2taxid
        print('Downloading dead_nucl.accession2taxid.gz... ', end="", flush=True)
        start_time = time()
        # Methods.download(DbBuilder.acc2taxid_url, self.output_folder + '/dead_nucl.accession2taxid.gz')
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        # Download nucleotide sequences of the query result
        seq_file = self.output_folder + '/seq.fasta'
        # Methods.download_seq_from_query(self.query, seq_file, self.email, self.api)

        # Extract accession numbers from downloaded fasta sequences
        print('Extracting accession numbers from fasta file...', end="", flush=True)
        start_time = time()
        acc_file = self.output_folder + '/acc.list'
        acc_dict = Methods.extract_acc_from_fasta(seq_file, acc_file)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        print('Parsing nucl_gb.accession2taxid.gz... ', end="", flush=True)
        start_time = time()
        acc2taxid_dict = Methods.parse_acc2taxid_file(self.output_folder + '/nucl_gb.accession2taxid.gz')
        # acc2taxid_dict = Methods.parse_acc2taxid_file_parallel(self.output_folder + '/nucl_gb.accession2taxid.gz', self.cpu)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        print('Parsing dead_nucl_gb.accession2taxid.gz... ', end="", flush=True)
        start_time = time()
        acc2taxid_dict.update(Methods.parse_acc2taxid_file(self.output_folder + '/dead_nucl.accession2taxid.gz'))
        # acc2taxid_dict.update(Methods.parse_acc2taxid_file_parallel(self.output_folder + '/dead_nucl.accession2taxid.gz', self.cpu))
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        # Find taxid for acc
        print('Finding taxID for accessions... ', end="", flush=True)
        start_time = time()
        taxid_file = self.output_folder + '/taxid.list'
        Methods.acc_to_taxid(acc_dict, acc2taxid_dict, taxid_file)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        # Expand taxonomy from taxid
        print('Writting taxonomy...', end="", flush=True)
        start_time = time()
        taxonomy_file = self.output_folder + '/taxonomy.txt'
        nodes_file = self.output_folder + '/nodes.dmp'
        names_file = self.output_folder + '/names.dmp'
        merged_file = self.output_folder + '/merged.dmp'
        Methods.expand_taxonomy_from_taxid(taxid_file, taxonomy_file, nodes_file, names_file, merged_file)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % Methods.elapsed_time(interval))

        # QIIME2
        qiime2_seq = seq_file + '.qza'
        qiime2_taxo = taxonomy_file + '.qza'

        print('Importing fasta file in QIIME2...')
        Methods.qiime2_import_unite_sequences(seq_file, qiime2_seq)

        print('Importing fasta file in QIIME2...')
        Methods.qiime2_import_unite_taxonomy(taxonomy_file, qiime2_taxo)

        print('Training classifier...')
        classifier_file = self.output_folder + '/seq_ncbi.qza'
        Methods.qiime2_train_unite_classifier(qiime2_seq, qiime2_taxo, classifier_file)

        end_time = time()
        interval = end_time - t_zero
        print("Done (total time: %s)." % Methods.elapsed_time(interval))

    def checks(self):
        # Check query string
        if not self.query:
            raise Exception('Your query is empty.')

        # Check output folder
        if not os.path.exists(self.output_folder):
            Methods.make_folder(self.output_folder)  # Create output folder is does not exists already

        # Check cpu and parallel processes
        cpu = cpu_count()
        if 1 > self.cpu > cpu:  # smaller than 1 or greater than available cpu
            self.cpu = cpu


if __name__ == '__main__':

    parser = ArgumentParser(description='Download DNA sequence from NCBI and add taxonomy for QIIME2.')
    parser.add_argument('-q', '--query', metavar='\"txid4762[Organism:exp] AND (\"internal transcribed spacer\"[Title]) NOT uncultured[Title]\"',
                        required=True,
                        type=str,
                        help='NCBI query string.'
                             ' Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        type=str,
                        help='Output folder.'
                             ' Mandatory.')
    parser.add_argument('-t', '--threads', metavar='4',
                        required=False, default=4,
                        type=int,
                        help='Number of CPU. Default is 4.')
    parser.add_argument('-e', '--email', metavar='your.email@example.org',
                        required=False, default='\'your.email@example.org\'',
                        type=str,
                        help='Your email address.')
    parser.add_argument('-a', '--api-key', metavar='',
                        required=False,
                        type=str,
                        help='Your NCBI API key. Allows up to 10 requests per second instead of 3.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    DbBuilder(arguments)
