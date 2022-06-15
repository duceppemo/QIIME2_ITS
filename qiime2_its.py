#!/usr/local/env python3

import os
import subprocess
from argparse import ArgumentParser
from multiprocessing import cpu_count
from qiime2_methods import Qiime2Methods
from collections import defaultdict


class Qiime2(object):
    def __init__(self, args):
        self.args = args

        # Input
        self.input_folder = args.input
        self.qiime2_classifier = args.classifier
        self.metadata_file = args.metadata

        # Read type
        self.single = args.se
        self.paired = args.pe

        # Output
        self.output_folder = args.output

        # Computing resources
        self.cpu = args.threads
        self.parallel = args.parallel_processes

        # QIIME2 environment
        self.qiime2_env = args.qiime2

        # Options
        self.revers_complement = args.reverse_complement
        self.min_len = args.min_len
        self.max_len = args.max_len

        # ITS
        self.its1 = args.extract_its1
        self.its2 = args.extract_its2
        self.taxa = args.taxa
        if self.its1:
            self.region = 'ITS1'
        elif self.its2:
            self.region = 'ITS2'

        # Data
        self.fastq_list = list()  # List of all the fastq files in the input folder
        self.sample_dict = dict  # { sample: [read location])

        # Install path
        self.install_path = os.path.dirname(__file__)

        # Run
        self.run()

    def run(self):
        # List fastq files from input folder
        self.fastq_list = Qiime2Methods.list_fastq(self.input_folder)

        # Check a few things before starting processing data
        self.checks()

        # Parse samples and their locations in a dictionary
        self.sample_dict = self.parse_fastq_list(self.fastq_list)

        if self.its1 or self.its2:  # Extracting ITS region from reads
            # Check for reverse complement flag
            if self.revers_complement:
                print('Reverse complementing reads...')
                rc_folder = self.output_folder + '/rc_reads'
                Qiime2Methods.make_folder(rc_folder)
                Qiime2.run_fastq_rc(self.input_folder, rc_folder)
                # Update fastq list abd input folder
                self.fastq_list = Qiime2Methods.list_fastq(rc_folder)

            # Create ITSxpress output folder and log folder
            print('Extracting ITS...')
            its_folder = self.output_folder + '/ITSxpress'
            itsxpress_log_folder = self.output_folder + '/itsxpress_log'
            Qiime2Methods.make_folder(its_folder)
            Qiime2Methods.make_folder(itsxpress_log_folder)

            # Extract Fungi ITS1 in parallel
            if self.single:
                Qiime2Methods.extract_its_se_parallel(self.fastq_list, its_folder, itsxpress_log_folder,
                                                      self.cpu, self.parallel, self.taxa, self.region)
            else:  # if self.paired:
                Qiime2Methods.extract_its_pe_parallel(self.sample_dict, its_folder, itsxpress_log_folder,
                                                      self.cpu, self.parallel, self.taxa, self.region)

            # Remove empty sequences. This is an artifact from ITSxpress.
            print('Checking for empty entries...')
            its_fastq_list = Qiime2Methods.list_fastq(its_folder)
            if self.single:
                Qiime2Methods.fix_fastq_se_parallel(its_fastq_list, self.cpu)

            else:
                its_sample_dict = self.parse_fastq_list(its_fastq_list)
                Qiime2Methods.fix_fastq_pe_parallel(its_sample_dict, self.cpu)
        else:  # Not extracting ITS from reads
            its_folder = self.input_folder

        # Filter reads by size
        if self.min_len or self.max_len:  # discard reads outside of range
            print('Filtering reads based on sizes...')
            size_folder = self.output_folder + '/size_filtered'
            Qiime2Methods.make_folder(size_folder)

            if self.single:
                my_fastq_list = Qiime2Methods.list_fastq(its_folder)
                Qiime2Methods.size_select_se_parallel(my_fastq_list, size_folder, self.min_len, self.max_len,
                                                      self.cpu, self.parallel)
            else:
                my_sample_dict = self.parse_fastq_list(its_fastq_list)
                Qiime2Methods.size_select_pe_parallel(my_sample_dict, size_folder, self.min_len, self.max_len,
                                                      self.cpu, self.parallel)
            # Make sure to use size filtered read folder for QIIME2
            its_folder = size_folder

        # Run QIIME2
        print('Running QIIME2...')
        print('\tImporting data...')
        if self.paired:
            Qiime2Methods.qiime2_import_fastq_pe(its_folder, self.output_folder + '/demux-seqs.qza')
        else:  # elif self.single
            Qiime2Methods.qiime2_import_fastq_se(its_folder, self.output_folder + '/demux-seqs.qza')
        Qiime2Methods.qiime2_demux_summary(self.output_folder + '/demux-seqs.qza',
                                           self.output_folder + '/demux-seqs.qzv')

        # Denoise with DADA2
        print('\tDenoising data with DADA2...')
        if self.paired:
            Qiime2Methods.qiime2_dada2_denoise_paired(self.output_folder + '/demux-seqs.qza',
                                                      self.output_folder + '/rep-seqs.qza',
                                                      self.output_folder + '/table.qza',
                                                      self.output_folder + '/stats.qza')
        else:
            Qiime2Methods.qiime2_dada2_denoise_single(self.output_folder + '/demux-seqs.qza',
                                                      self.output_folder + '/rep-seqs.qza',
                                                      self.output_folder + '/table.qza',
                                                      self.output_folder + '/stats.qza')
        # Tabulate
        Qiime2Methods.qiime2_metadata_tabulate(self.output_folder + '/stats.qza',
                                               self.output_folder + '/stats.qzv')

        # Export BIOM table
        print('\tExporting BIOM table...')
        Qiime2Methods.qiime2_export(self.output_folder + '/table.qza', self.output_folder + '/biom_table')

        # Feature table
        print('\tCreating representative sequences...')
        Qiime2Methods.qiime2_sample_summarize(self.metadata_file,
                                              self.output_folder + '/table.qza',
                                              self.output_folder + '/table.qzv')

        Qiime2Methods.qiime2_seq_sumamry(self.output_folder + '/rep-seqs.qza',
                                         self.output_folder + '/rep-seqs.qzv')

        # Phylogeny
        print('\tAligning representative sequences and making phylogenic tree...')
        Qiime2Methods.qiime2_phylogeny(self.output_folder + '/rep-seqs.qza',
                                       self.output_folder + '/aligned-rep-seqs.qza',
                                       self.output_folder + '/masked-aligned-rep-seqs.qza',
                                       self.output_folder + '/unrooted-tree.qza',
                                       self.output_folder + '/rooted-tree.qza')

        # Export phylogenetic tree
        Qiime2Methods.qiime2_export(self.output_folder + '/unrooted-tree.qza',
                                    self.output_folder)

        # Alpha and Beta diversity
        print('\tAnalysing Alpha and Beta diversity...')
        Qiime2Methods.qiime2_core_diversity(self.cpu,
                                            self.metadata_file,
                                            self.output_folder + '/rooted-tree.qza',
                                            self.output_folder + '/table.qza',
                                            self.output_folder + '/core-metrics-results')

        # Alpha rarefaction plotting
        print('\tCreating rarefaction plot...')
        Qiime2Methods.qiime2_rarefaction(self.metadata_file,
                                         self.output_folder + '/rooted-tree.qza',
                                         self.output_folder + '/table.qza',
                                         self.output_folder + '/alpha-rarefaction.qzv')

        # Taxonomic analysis
        print('\tAssigning taxonomy to representative sequences...')
        Qiime2Methods.qiime2_classify(self.qiime2_classifier,
                                      self.output_folder + '/rep-seqs.qza',
                                      self.output_folder + '/taxonomy.qza')

        # Create taxonomy table
        print('\tCreating taxonomy table...')
        Qiime2Methods.qiime2_metadata_tabulate(self.output_folder + '/taxonomy.qza',
                                               self.output_folder + '/taxonomy.qzv')

        # Export taxonomy table
        print('\tExporting taxonomy...')
        Qiime2Methods.qiime2_export(self.output_folder + '/taxonomy.qza',
                                    self.output_folder + '/biom_table')

        # Merge taxonomy yo bio table
        print('\tIncorporating taxonomy to BIOM table...')
        Qiime2Methods.change_taxonomy_file_header(self.output_folder + '/biom_table/taxonomy.tsv')
        Qiime2Methods.biom_add_metadata(self.output_folder + '/biom_table/feature-table.biom',
                                        self.output_folder + '/biom_table/taxonomy.tsv',
                                        self.output_folder + '/biom_table/table-with-taxonomy.biom')

        # Convert to biom with taxonomy to tsv
        print('\tExporting BIOM table with taxonomy...')
        Qiime2Methods.biom_convert_taxo(self.output_folder + '/biom_table/table-with-taxonomy.biom',
                                        self.output_folder + '/biom_table/taxonomy.tsv',
                                        self.output_folder + '/biom_table/table-with-taxonomy.biom.tsv')

        # Bar plot of sample composition
        print('Creating bar plot of sample composition...')
        Qiime2Methods.qiime2_taxa_barplot(self.output_folder + '/table.qza',
                                          self.output_folder + '/taxonomy.qza',
                                          self.metadata_file,
                                          self.output_folder + '/taxa-bar-plots.qzv')

        print('DONE!')

    @staticmethod
    def parse_fastq_list(fastq_list):
        sample_dict = defaultdict(list)
        for fq in fastq_list:
            sample = os.path.basename(fq).split('_')[0]
            # sample = ''
            # if fq.endswith('.gz'):
            #     sample = '.'.join(os.path.basename(fq).split('.')[:-2])
            # else:
            #     sample = '.'.join(os.path.basename(fq).split('.')[:-1])
            if 'R1' in os.path.basename(fq).split('_')[3]:
                sample_dict[sample].insert(0, fq)
            elif 'R2' in os.path.basename(fq).split('_')[3]:
                sample_dict[sample].insert(1, fq)

        return sample_dict

    def checks(self):
        # Check if paired-end of single-end flag used
        if not (self.paired or self.single):
            raise Exception('You must state if reads are single-end or paired end ("-pe" or "-se")')

        # Check output folder
        if not os.path.exists(self.output_folder):
            Qiime2Methods.make_folder(self.output_folder)  # Create output folder is does not exists already

        # Check input folder
        if not self.fastq_list:
            raise Exception('No fastq files found in provided input folder.')

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
        error_message = ['File name must be as is:',
                         '\t1. the sample identifier,',
                         '\t2. the barcode sequence or a barcode identifier,',
                         '\t3. the lane number starting with "L" followed by 3 digits,',
                         '\t4. the direction of the read ("R1" or "R2"), and',
                         '\t5. the set number (always "001").']

        for f in self.fastq_list:
            base_name = os.path.basename(f)
            fields = base_name.split('.')[0].split('_')
            ext = ''
            if base_name.endswith('gz'):
                ext = base_name.split(".")[-2]
            else:
                ext = base_name.split(".")[-1]

            if ext != 'fastq':
                raise Exception('\n'.join(error_message))

            if len(fields) != 5:
                raise Exception('\n'.join(error_message))
            else:
                if not fields[2].startswith('L'):
                    raise Exception('\n'.join(error_message))
                elif fields[3] not in ['R1', 'R2']:
                    raise Exception('\n'.join(error_message))
                elif not fields[4] == '001':
                    raise Exception('\n'.join(error_message))

        # Check metadata
        # Very basic check. Only checking that sample names are the same as fastq files

        # Paired end
        # Check if there are 2 files per sample
        if self.paired:
            pass

        # Check that input folder does not contain both read types (paired-end and single-end)

        # Check cpu and parallel processes
        cpu = cpu_count()
        if 1 > self.cpu > cpu:  # smaller than 1 or greater than available cpu
            self.cpu = cpu
        if self.parallel > self.cpu:
            self.parallel = self.cpu  # Will use only 1 cpu per parallel process

        # Check extract-its flag
        if self.its1 and self.its2:
            raise Exception('You cannot choose both ITS1 and ITS2 during the same analysis.')

        # Check taxa and region
        if self.its1 or self.its2:
            if self.region:
                if self.region not in ['ITS2', 'ITS1', 'ALL']:
                    raise Exception('Taxa must be "ITS2", "ITS1" or "ALL".')
            if self.taxa:
                if self.taxa not in ['Alveolata', 'Bryophyta', 'Bacillariophyta', 'Amoebozoa', 'Euglenozoa',
                                     'Fungi', 'Chlorophyta', 'Rhodophyta', 'Phaeophyceae', 'Marchantiophyta',
                                     'Metazoa', 'Oomycota', 'Haptophyceae', 'Raphidophyceae', ' Rhizaria',
                                     'Synurophyceae', 'Tracheophyta', 'Eustigmatophyceae', 'All']:
                    raise Exception('Region must be "Alveolata", "Bryophyta", "Bacillariophyta", "Amoebozoa", '
                                    '"Euglenozoa", "Fungi", "Chlorophyta", "Rhodophyta", "Phaeophyceae", '
                                    '"Marchantiophyta", "Metazoa", "Oomycota", "Haptophyceae", "Raphidophyceae", '
                                    '" Rhizaria", "Synurophyceae", "Tracheophyta", "Eustigmatophyceae" or "All".')

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

    parser = ArgumentParser(description='Run QIIME2 on IonTorrent sequencing data using the UNITE database')
    parser.add_argument('-q', '--qiime2', metavar='qiime2-2020.8',
                        required=True,
                        type=str,
                        help='Name of your QIIME2 conda environment.'
                             ' Mandatory.')
    parser.add_argument('-i', '--input', metavar='/input_folder/',
                        required=True,
                        type=str,
                        help='Input folder where the fastq reads are located.'
                             ' Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        type=str,
                        help='Output folder for QIIME2 files.'
                             ' Mandatory.')
    parser.add_argument('-m', '--metadata', metavar='qiime2_metadata.tsv',
                        required=True,
                        type=str,
                        help='Validated QIIME2 metadata file (samples description).'
                             ' Mandatory.')
    parser.add_argument('-c', '--classifier', metavar='unite_classifier_qiime2.qza',
                        required=True,
                        type=str,
                        help='Classifier for QIIME2. See script "train_unite_classifier_qiime2.py" to compile it.'
                             ' Mandatory.')
    parser.add_argument('-t', '--threads', metavar='4',
                        required=False, default=4,
                        type=int,
                        help='Number of CPU. Default is 4.')
    parser.add_argument('-p', '--parallel-processes', metavar='1',
                        required=False, default=1,
                        type=int,
                        help='Processes to run in parallel. Adjust according the number of threads. '
                             'For example, if 16 threads, using 4 parallel processes will run 4 samples in parallel '
                             'using 4 threads each (16/4). Default is 1.')
    parser.add_argument('-rc', '--reverse_complement',
                        action='store_true',  # 'store_true' means False by default
                        required=False,
                        help='Use this flag is your reads are in reverse complement. '
                             'For example if you sequenced from 5.8S to 18S.'
                             ' Optional.')
    parser.add_argument('--min-len',
                        type=int, default=0,
                        required=False,
                        help='Minimum read length to keep. Default is 0 (no min length).'
                             ' Optional.')
    parser.add_argument('--max-len',
                        type=int, default=0,
                        required=False,
                        help='Maximum read length to keep. Default is 0 (no max length).'
                             ' Optional.')
    parser.add_argument('-se',
                        required=False,
                        action='store_true',
                        help='Reads are single-end. One fastq file per sample.'
                             ' "-se" or "-pe" mandatory.')
    parser.add_argument('-pe',
                        required=False,
                        action='store_true',
                        help='Reads are paired-end. Two fastq file per sample.'
                             ' "-se" or "-pe" mandatory.')
    parser.add_argument('--extract-its1',
                        required=False,
                        action='store_true',
                        help='Extract ITS1 sequence from reads with ITSxpress. Cannot be used with "--extract-its2".')
    parser.add_argument('--extract-its2',
                        required=False,
                        action='store_true',
                        help='Extract ITS2 sequence from reads with ITSxpress. Cannot be used with "--extract-its1".')
    parser.add_argument('--taxa', metavar='Fungi',
                        required=False, default='Fungi',
                        type=str,
                        help='Select taxa of interest for ITSxpress: {Alveolata,Bryophyta,Bacillariophyta,Amoebozoa,'
                             'Euglenozoa,Fungi,Chlorophyta,Rhodophyta,Phaeophyceae,Marchantiophyta,Metazoa,Oomycota,'
                             'Haptophyceae,Raphidophyceae, Rhizaria,Synurophyceae,Tracheophyta,Eustigmatophyceae,All}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Qiime2(arguments)
