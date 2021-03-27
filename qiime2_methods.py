#!/usr/local/env python3

import os
import pathlib
import subprocess
from concurrent import futures


class Qiime2Methods(object):

    @staticmethod
    def make_folder(folder):
        """
        Create output folder.
        :param folder: string. Output folder path.
        :return:
        """
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def list_fastq(my_path):
        """
        Walk input directory and list all the fastq files. Accepted file extensions are '.fastq', '.fastq.gz',
        '.fq' and '.fq.gz'.
        :param my_path: string. Input folder path
        :return: list of strings. Fastq files in input folder
        """
        # Create empty list to hold the file paths
        fastq_list = list()
        # Walk the input directory recursively and look for fastq files
        for root, directories, filenames in os.walk(my_path):
            for filename in filenames:
                absolute_path = os.path.join(root, filename)
                if os.path.isfile(absolute_path) and filename.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                    fastq_list.append(absolute_path)  # Add fastq file path to the list
        return fastq_list

    @staticmethod
    def extract_its(fastq_file, output_folder, cpu):
        cmd = ['itsxpress',
               '--threads', str(cpu),
               '--single_end',
               '--fastq', fastq_file,
               '--region', 'ITS1',
               '--taxa', 'Fungi',
               '--cluster_id', str(0.99),
               '--outfile', output_folder + '/' + os.path.basename(fastq_file),
               '--log',  output_folder + '/log/' + os.path.basename(fastq_file).split('_')[0] + '.log',
               '--threads', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def extract_its_parallel(fastq_list, output_folder, cpu):
        with futures.ThreadPoolExecutor(max_workers=cpu / 4) as executor:
            args = ((fastq, output_folder, cpu / 4) for fastq in fastq_list)
            for results in executor.map(lambda p: Qiime2Methods.extract_its(*p), args):  # (*p) does the unpacking part
                pass

    @staticmethod
    def fix_fastq(fastq_file):
        cmd = ['python', ]

    @staticmethod
    def qiime2_import_fastq(input_path, output_path):
        """
        Import fastq files
        https://docs.qiime2.org/2020.8/tutorials/importing/
        :param input_path:
        :param output_path:
        :return:
        """
        cmd = ['qiime', 'tools', 'import',
               '--type', 'SampleData[SequencesWithQuality]',
               '--input-format', 'CasavaOneEightSingleLanePerSampleDirFmt',  # For demultiplexed single end fastq
               '--input-path', input_path,
               '--output-path', output_path]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_demux_summary(input_folder, output_folder):
        """
        Make summary of samples
        Subsample 10,000 reads by default, only use 1000 instead (faster)
        :param input_folder:
        :param output_folder:
        :return:
        """
        cmd = ['qiime', 'demux', 'summary',
               '--p-n', str(1000),
               '--i-data', input_folder,
               '--o-visualization', output_folder + '/demux-single-end.qzv']
        subprocess.run(cmd)

    @staticmethod
    def qiime2_dada2_denoise_single(qiime2_reads, output_folder):
        """
        Denoise reads with DADA2
        :param qiime2_reads: string. qiime2 qza file containing reads info
        :param output_folder: string. Output folder
        :return:
        """
        cmd = ['qiime', 'dada2', 'denoise-single',
               '--p-n-threads', str(0),
               '--p-trim-left', str(0),  # No trimming
               '--p-trunc-len', str(0),  # No trimming
               '--i-demultiplexed-seqs', qiime2_reads,
               '--o-representative-sequences', output_folder + '/rep-seqs.qza',
               '--o-table', output_folder + '/table.qza',
               '--o-denoising-stats', output_folder + '/stats.qza']
        subprocess.run(cmd)

    @staticmethod
    def qiime2_metadata_tabulate(qiime2_stats, output_folder):
        """

        :param qiime2_stats:
        :param output_folder:
        :return:
        """
        cmd = ['qiime', 'metadata', 'tabulate',
               '--m-input-file',  qiime2_stats,
               '--o-visualization',  output_folder + '/stats.qzv']
        subprocess.run(cmd)

    @staticmethod
    def eqiime2_export(input, output_folder):
        """
        Export biom table
        :param input: string. QIIME2 table
        :param output_folder: sting. Output folder
        :return:
        """
        cmd = ['qiime', 'tools', 'export',
               '--input-path',  input,
               '--output-path',  output_folder + '/biom_table']  # this is a folder
        subprocess.run(cmd)

    @staticmethod
    def qiime2_sample_summary(metadata_file, qiime2_table, output_folder):
        """

        :param metadata_file:
        :param qiime2_table:
        :param output_folder:
        :return:
        """
        cmd = ['qiime', 'feature-table', 'summarize',
               '--m-sample-metadata-file', metadata_file,
               '--i-table', qiime2_table,
               '--o-visualization',  output_folder + '/table.qzv']
        subprocess.run(cmd)

    @staticmethod
    def qiime2_seq_sumamry(qiime2_repseqs, output_folder):
        """

        :param qiime2_repseqs:
        :param output_folder:
        :return:
        """
        cmd = ['qiime', 'feature-table', 'tabulate-seqs',
               '--i-data', qiime2_repseqs,
               '--o-visualization',  output_folder + '/rep-seqs.qzv']
        subprocess.run(cmd)

    @staticmethod
    def qiime2_phylogeny(qiime2_repseqs, output_folder):
        """

        :param qiime2_repseqs:
        :param output_folder:
        :return:
        """
        cmd = ['qiime', 'phylogeny', 'align-to-tree-mafft-fasttree',
               '--p-n-threads', 'auto',
               '--i-sequences', qiime2_repseqs,
               '--o-alignment', output_folder + '/aligned-rep-seqs.qza',
               '--o-masked-alignment', output_folder + '/masked-aligned-rep-seqs.qza',
               '--o-tree', output_folder + '/unrooted-tree.qza',
               '--o-rooted-tree', output_folder + '/rooted-tree.qza']
        subprocess.run(cmd)

    @staticmethod
    def qiime2_core_diversity(cpu, metadata_file, qiime2_rooted_tree, qiime2_table, output_folder):
        # Alpha and Beta analysis
        cmd = ['qiime', 'diversity', 'core-metrics-phylogenetic',
               '--p-n-jobs-or-threads', str(cpu),
               '--p-sampling-depth', str(1000),
               '--i-phylogeny', qiime2_rooted_tree,
               '--i-table', qiime2_table,
               '--m-metadata-file', metadata_file,
               '--output-dir', output_folder + '/core-metrics-results']

    @staticmethod
    def qiime2_rarefaction(metadata_file, qiime2_rooted_tree, qiime2_table, output_folder):
        cmd = ['qiime', 'diversity', 'alpha-rarefaction',
               '--p-max-depth', str(4000),
               '--i-phylogeny', qiime2_rooted_tree,
               '--i-table', qiime2_table,
               '--m-metadata-file', metadata_file,
               '--o-visualization', output_folder + '/alpha-rarefaction.qzv']
        subprocess.run(cmd)

    @staticmethod
    def qiime2_classify(qiime2_classifier, qiime2_repseqs, output_folder): # Taxonomic analysis
        cmd = ['qiime', 'feature-classifier', 'classify-sklearn',
               '--p-n-jobs', str(-1),
               '--i-classifier', qiime2_classifier,
               '--i-reads', qiime2_repseqs,
               '--o-classification',  output_folder + '/taxonomy.qza']
        subprocess.run(cmd)

    @staticmethod
    def change_taxonomy_file_header(input_taxo, output_taxo):
        with open(output_taxo, 'w') as out_f:
            out_f.write('#OTUID\ttaxonomy\tconfidence\n')  # write new header
            with open(input_taxo, 'r') as in_f:
                next(in_f)  # skip header
                for line in in_f:
                    out_f.write(line)  # dump rest of file

    @staticmethod
    def biom_add_metadata(input_biom, taxonomy_tsv, output_biom):
        cmd = ['biom', 'add-metadata',
               '--sc-separated', 'taxonomy',
               '-i', input_biom,
               '--observation-metadata-fp', taxonomy_tsv,
               '-o', output_biom]
        subprocess.run(cmd)

    @staticmethod
    def biom_convert_taxo(input_biom_taxo, taxonomy_tsv, output_tsv_taxo):
        cmd = ['biom', 'convert',
               '--to-tsv',
               '--header-key taxonomy',
               '-i', input_biom_taxo,
               '--observation-metadata-fp', taxonomy_tsv,
               '-o', output_tsv_taxo]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_taxa_barplot(qiime2_table, qiime2_taxonomy, metadata_file, output_folder):
        cmd = ['qiime', 'taxa', 'barplot',
               '--i-table', qiime2_table,
               '--i-taxonomy', qiime2_taxonomy,
               '--m-metadata-file', metadata_file,
               '--o-visualization', output_folder + '/taxa-bar-plots.qzv']
        subprocess.run(cmd)
