#!/usr/local/env python3

import os
import pathlib
import subprocess
from concurrent import futures
from collections import OrderedDict
import gzip


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
    def rc_fastq(input_fastq, output_fastq):
        """

        :param input_fastq:
        :param output_fastq:
        :return:
        """
        cmd = ['reformat.sh',
               'ow=t',
               'rcomp=t',
               'in={}'.format(input_fastq),
               'out={}'.format(output_fastq)]
        subprocess.run(cmd)

    @staticmethod
    def rc_fastq_parallel(fastq_list, output_folder, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((fastq, output_folder, int(cpu/parallel)) for fastq in fastq_list)
            for results in executor.map(lambda p: Qiime2Methods.extract_its_se(*p), args):  # (*p) unpacks
                pass

    @staticmethod
    def extract_its_se(fastq_file, output_folder, log_folder, cpu, taxa, region):
        """
        Extract Fungi ITS1 sequence from fastq files. Use ITSxpress program.

        :param fastq_file: string. Fastq file path.
        :param output_folder: string. Path of output folder.
        :param log_folder: sting. Path of log folder.
        :param cpu: int. Number of CPU to use.
        :param taxa: str. Taxonomy of interest.
        :param region: string. Region of interest.
        :return:
        """
        cmd = ['itsxpress',
               '--threads', str(cpu),
               '--single_end',
               '--fastq', fastq_file,
               '--region', region,
               '--taxa', taxa,
               '--cluster_id', str(0.99),
               '--outfile', output_folder + '/' + os.path.basename(fastq_file),
               '--log',  log_folder + '/' + os.path.basename(fastq_file).split('_')[0] + '.log',
               '--threads', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def extract_its_se_parallel(fastq_list, output_folder, log_folder, cpu, parallel, taxa, region):
        """
        Run "extract_its" in parallel using 4 cores per instance.
        :param fastq_list: string. A list of fastq file paths.
        :param output_folder: sting. Path of output folder
        :param log_folder: sting. Path of log folder.
        :param cpu: int. Number of cpu to use.
        :param parallel: int. Number of samples to process in parallel
        :param taxa: str. Taxonomy of interest.
        :param region: string. Region of interest.
        :return:
        """
        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((fastq, output_folder, log_folder, int(cpu/parallel), taxa, region) for fastq in fastq_list)
            for results in executor.map(lambda p: Qiime2Methods.extract_its_se(*p), args):  # (*p) unpacks
                pass

    @staticmethod
    def extract_its_pe(fastq_r1, fastq_r2, output_folder, log_folder, cpu, taxa, region):
        """
        Extract Fungi ITS1 sequence from fastq files. Use ITSxpress program.
        :param fastq_r1: string. Fastq file path
        :param fastq_r2: string. Fastq file path
        :param output_folder: string. Path of output folder.
        :param log_folder: sting. Path of log folder.
        :param cpu: int. number of CPU to use.
        :param taxa: str. Taxonomy of interest.
        :param region: string. Region of interest.
        :return:
        """
        cmd = ['itsxpress',
               '--threads', str(cpu),
               '--fastq', fastq_r1,
               '--fastq2', fastq_r2,
               '--region', region,
               '--taxa', taxa,
               '--cluster_id', str(0.99),
               '--outfile', output_folder + '/' + os.path.basename(fastq_r1),
               '--outfile2', output_folder + '/' + os.path.basename(fastq_r2),
               '--log',  log_folder + '/' + os.path.basename(fastq_r1).split('_')[0] + '.log',
               '--threads', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def extract_its_pe_parallel(sample_dict, output_folder, log_folder, cpu, parallel, taxa, region):
        """
        Run "extract_its" in parallel using 4 cores per instance.
        :param sample_dict: string. A dictionary of fastq file paths.
        :param output_folder: sting. Path of output folder
        :param log_folder: sting. Path of log folder.
        :param cpu: int. Number of cpu to use.
        :param parallel: int. Number of samples to process in parallel
        :param taxa: str. Taxonomy of interest.
        :param region: string. Region of interest.
        :return:
        """
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((fastq_list[0], fastq_list[1], output_folder, log_folder, int(cpu/parallel), taxa, region)
                    for sample, fastq_list in sample_dict.items())
            for results in executor.map(lambda p: Qiime2Methods.extract_its_pe(*p), args):  # (*p) unpacks arguments
                pass

    @staticmethod
    def remove_empties_se(input_fastq):
        """
        Remove empty entries (header present, but sequence and q phred scores empty) from single fastq file.
        :param input_fastq: string. Path to single-end fastq file.
        :return:
        """
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
        """
        Remove empty entries (header present, but sequence and q phred scores empty) from paired-end fastq file.
        Keeps both files synchronized.
        :param r1: string. Path to R1 fastq file from pair.
        :param r2: string. Path to R2 fastq file from pair.
        :return:
        """
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

    @staticmethod
    def fix_fastq_se_parallel(fastq_list, parallel):
        """
        Run "fix_fastq" in parallel using all the threads, one file per thread
        :param fastq_list: string. list of fastq file paths
        :param parallel: int. number of samples to process in parallel
        :return:
        """
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            for results in executor.map(Qiime2Methods.remove_empties_se, fastq_list):
                pass

    @staticmethod
    def fix_fastq_pe_parallel(sample_dict, parallel):
        """
        Run "fix_fastq" in parallel using all the threads, one file per thread
        :param sample_dict: string. Dictionary of fastq file paths
        :param parallel: int. number of samples to process in parallel
        :return:
        """
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((fastq_list[0], fastq_list[1]) for sample, fastq_list in sample_dict.items())
            for results in executor.map(lambda p: Qiime2Methods.remove_empties_pe(*p), args):  # (*p) unpacks
                pass

    @staticmethod
    def size_select_se(r1, out_folder, min_len, max_len, cpu):
        cmd = ['bbduk.sh',
               'in={}'.format(r1),
               'out={}'.format(out_folder + '/' + os.path.basename(r1)),
               'minlength={}'.format(min_len),
               'maxlength={}'.format(max_len),
               'threads={}'.format(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def size_select_pe(r1, r2, out_folder, min_len, max_len, cpu):
        cmd = ['bbduk.sh',
               'in={}'.format(r1),
               'in2={}'.format(r2),
               'out={}'.format(out_folder + '/' + os.path.basename(r1)),
               'out2={}'.format(out_folder + '/' + os.path.basename(r2)),
               'minlength={}'.format(min_len),
               'maxlength={}'.format(max_len),
               'threads={}'.format(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def size_select_se_parallel(fastq_list, out_folder, min_len, max_len, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((fastq, out_folder, min_len, max_len, int(cpu/parallel)) for fastq in fastq_list)
            for results in executor.map(lambda p: Qiime2Methods.size_select_se(*p), args):  # (*p) unpacks arguments
                pass

    @staticmethod
    def size_select_pe_parallel(sample_dict, out_folder, min_len, max_len, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((fastq_list[0], fastq_list[1], out_folder, min_len, max_len, int(cpu/parallel))
                    for sample, fastq_list in sample_dict.items())
            for results in executor.map(lambda p: Qiime2Methods.size_select_pe(*p), args):  # (*p) unpacks arguments
                pass

    @staticmethod
    def qiime2_import_fastq_se(fastq_folder, reads_qza):
        """
        Import single-end fastq files
        https://docs.qiime2.org/2020.8/tutorials/importing/
        :param fastq_folder:
        :param reads_qza:
        :return:
        """
        cmd = ['qiime', 'tools', 'import',
               '--type', 'SampleData[SequencesWithQuality]',
               '--input-format', 'CasavaOneEightSingleLanePerSampleDirFmt',  # For demultiplexed single end fastq
               '--input-path', fastq_folder,
               '--output-path', reads_qza]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_import_fastq_pe(fastq_folder, reads_qza):
        """
        Import paired-end fastq files
        https://docs.qiime2.org/2020.8/tutorials/importing/
        :param fastq_folder:
        :param reads_qza:
        :return:
        """
        cmd = ['qiime', 'tools', 'import',
               '--type', 'SampleData[PairedEndSequencesWithQuality]',
               '--input-format', 'CasavaOneEightSingleLanePerSampleDirFmt',  # For demultiplexed paired-end fastq
               '--input-path', fastq_folder,
               '--output-path', reads_qza]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_demux_summary(reads_qza, output_qzv):
        """
        Make summary of samples
        Subsample 10,000 reads by default, only use 1000 instead (faster)
        :param reads_qza:
        :param output_qzv:
        :return:
        """
        cmd = ['qiime', 'demux', 'summarize',
               '--p-n', str(1000),
               '--i-data', reads_qza,
               '--o-visualization', output_qzv]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_dada2_denoise_single(reads_qza, repseq_qza, table_qza, stats_qza):
        """
        Denoise single-end reads with DADA2. Optimized for Illumina.
        :param reads_qza:
        :param repseq_qza:
        :param table_qza:
        :param stats_qza:
        :return:
        """
        cmd = ['qiime', 'dada2', 'denoise-single',
               '--p-n-threads', str(0),
               '--p-trim-left', str(0),  # No trimming
               '--p-trunc-len', str(0),  # No trimming
               '--i-demultiplexed-seqs', reads_qza,
               '--o-representative-sequences', repseq_qza,
               '--o-table', table_qza,
               '--o-denoising-stats', stats_qza]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_dada2_denoise_paired(reads_qza, repseq_qza, table_qza, stats_qza):
        """
        Denoise paired-end reads with DADA2. Optimized for Illumina.
        :param reads_qza:
        :param repseq_qza:
        :param table_qza:
        :param stats_qza:
        :return:
        """
        cmd = ['qiime', 'dada2', 'denoise-paired',
               '--p-n-threads', str(0),
               '--p-trim-left-f', str(0),  # No trimming
               '--p-trim-left-r', str(0),  # No trimming
               '--p-trunc-len-f', str(0),  # No trimming
               '--p-trunc-len-r', str(0),  # No trimming
               '--i-demultiplexed-seqs', reads_qza,
               '--o-representative-sequences', repseq_qza,
               '--o-table', table_qza,
               '--o-denoising-stats', stats_qza]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_metadata_tabulate(stats_qza, stats_qzv):
        """

        :param stats_qza:
        :param stats_qzv:
        :return:
        """
        cmd = ['qiime', 'metadata', 'tabulate',
               '--m-input-file',  stats_qza,
               '--o-visualization',  stats_qzv]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_export(qza, output_folder):
        """
        Export biom table
        :param qza: string. QIIME2 table
        :param output_folder: sting. Output folder
        :return:
        """
        cmd = ['qiime', 'tools', 'export',
               '--input-path',  qza,
               '--output-path', output_folder]  # this is a folder
        subprocess.run(cmd)

    @staticmethod
    def qiime2_sample_summarize(metadata_file, table_qza, table_qzv):
        """

        :param metadata_file:
        :param table_qza:
        :param table_qzv:
        :return:
        """
        cmd = ['qiime', 'feature-table', 'summarize',
               '--m-sample-metadata-file', metadata_file,
               '--i-table', table_qza,
               '--o-visualization',  table_qzv]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_seq_sumamry(repseqs_qza, repseqs_qzv):
        """

        :param repseqs_qza:
        :param repseqs_qzv:
        :return:
        """
        cmd = ['qiime', 'feature-table', 'tabulate-seqs',
               '--i-data', repseqs_qza,
               '--o-visualization',  repseqs_qzv]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_phylogeny(repseqs_qza, align_repseqs_qza, masked_align_repseqs_qza, unrooted_tree_qza, rooted_tree_qza):
        """

        :param repseqs_qza:
        :param align_repseqs_qza:
        :param masked_align_repseqs_qza:
        :param unrooted_tree_qza:
        :param rooted_tree_qza:
        :return:
        """
        cmd = ['qiime', 'phylogeny', 'align-to-tree-mafft-fasttree',
               '--p-n-threads', 'auto',
               '--i-sequences', repseqs_qza,
               '--o-alignment', align_repseqs_qza,
               '--o-masked-alignment', masked_align_repseqs_qza,
               '--o-tree', unrooted_tree_qza,
               '--o-rooted-tree', rooted_tree_qza]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_core_diversity(cpu, metadata_file, rooted_tree_qza, table_qza, output_folder):
        """
        Alpha and Beta analysis
        :param cpu:
        :param metadata_file:
        :param rooted_tree_qza:
        :param table_qza:
        :param output_folder:
        :return:
        """
        cmd = ['qiime', 'diversity', 'core-metrics-phylogenetic',
               '--p-n-jobs-or-threads', str(cpu),
               '--p-sampling-depth', str(1000),
               '--i-phylogeny', rooted_tree_qza,
               '--i-table', table_qza,
               '--m-metadata-file', metadata_file,
               '--output-dir', output_folder + '/core-metrics-results']

    @staticmethod
    def qiime2_rarefaction(metadata_file, rooted_tree_qza, table_qza, rare_qzv):
        """

        :param metadata_file:
        :param rooted_tree_qza:
        :param table_qza:
        :param rare_qzv:
        :return:
        """
        cmd = ['qiime', 'diversity', 'alpha-rarefaction',
               '--p-max-depth', str(4000),
               '--i-phylogeny', rooted_tree_qza,
               '--i-table', table_qza,
               '--m-metadata-file', metadata_file,
               '--o-visualization', rare_qzv]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_classify(qiime2_classifier, repseqs_qza, taxonomy_qza):
        """
        Taxonomic analysis
        :param qiime2_classifier:
        :param repseqs_qza:
        :param taxonomy_qza:
        :return:
        """
        cmd = ['qiime', 'feature-classifier', 'classify-sklearn',
               '--p-n-jobs', str(-1),
               '--i-classifier', qiime2_classifier,
               '--i-reads', repseqs_qza,
               '--o-classification',  taxonomy_qza]
        subprocess.run(cmd)

    @staticmethod
    def change_taxonomy_file_header(input_taxo):
        """

        :param input_taxo:
        :return:
        """
        tmp = input_taxo + '.tmp'
        with open(tmp, 'w') as out_f:
            out_f.write('#OTUID\ttaxonomy\tconfidence\n')  # write new header
            with open(input_taxo, 'r') as in_f:
                next(in_f)  # skip header
                for line in in_f:
                    out_f.write(line)  # dump rest of file
        # overwrite original taxonomy file
        os.replace(tmp, input_taxo)

    @staticmethod
    def biom_add_metadata(input_biom, taxonomy_tsv, output_biom):
        """

        :param input_biom:
        :param taxonomy_tsv:
        :param output_biom:
        :return:
        """
        cmd = ['biom', 'add-metadata',
               '--sc-separated', 'taxonomy',
               '-i', input_biom,
               '--observation-metadata-fp', taxonomy_tsv,
               '-o', output_biom]
        subprocess.run(cmd)

    @staticmethod
    def biom_convert_taxo(input_biom_taxo, taxonomy_tsv, output_tsv_taxo):
        """

        :param input_biom_taxo:
        :param taxonomy_tsv:
        :param output_tsv_taxo:
        :return:
        """
        cmd = ['biom', 'convert',
               '--to-tsv',
               '--header-key', 'taxonomy',
               '-i', input_biom_taxo,
               '--observation-metadata-fp', taxonomy_tsv,
               '-o', output_tsv_taxo]
        subprocess.run(cmd)

    @staticmethod
    def qiime2_taxa_barplot(table_qza, taxonomy_qza, metadata_file, taxo_bar_plot_qzv):
        """

        :param table_qza:
        :param taxonomy_qza:
        :param metadata_file:
        :param taxo_bar_plot_qzv:
        :return:
        """
        cmd = ['qiime', 'taxa', 'barplot',
               '--i-table', table_qza,
               '--i-taxonomy', taxonomy_qza,
               '--m-metadata-file', metadata_file,
               '--o-visualization', taxo_bar_plot_qzv]
        subprocess.run(cmd)
