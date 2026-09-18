"""Main ITS pipeline: import -> (optional) ITSxpress trim -> (optional) size
filter -> DADA2 denoise -> phylogeny -> diversity -> taxonomy -> barplot ->
(optional) diversity/composition/classifier stats -> (optional) PDF report.
"""
import argparse
import getpass
import platform
import socket
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from qiime2_its import (biom_utils, env_checks, fastq_utils, itsxpress_wrapper, metadata_utils, provenance,
                         qiime_wrapper, report, report_data, size_filter, taxonomy)
from qiime2_its._version import __version__
from qiime2_its.itsxpress_wrapper import TAXA_CODES


class Pipeline:
    def __init__(self, args):
        self.start_time = datetime.now().astimezone()
        self.cli_args = vars(args)

        self.input_folder = Path(args.input)
        self.qiime2_classifier = args.classifier
        self.metadata_file = args.metadata

        self.single = args.se
        self.paired = args.pe

        self.output_folder = Path(args.output)

        self.cpu = env_checks.clamp_cpu(args.threads)
        self.parallel = env_checks.clamp_parallel(args.parallel_processes, self.cpu)

        self.qiime2_env = args.qiime2

        self.reverse_complement = args.reverse_complement
        self.min_len = args.min_len
        self.max_len = args.max_len

        self.its1 = args.extract_its1
        self.its2 = args.extract_its2
        self.taxa = args.taxa
        self.region = 'ITS1' if self.its1 else 'ITS2' if self.its2 else None

        # DADA2 denoising. Defaults match qiime2's own (Illumina-tuned)
        # defaults; single-end/IonTorrent users typically want a higher
        # --max-ee (noisier reads) and --allow-one-off (homopolymer indels
        # otherwise get miscalled as one-off bimeras).
        self.max_ee = args.max_ee
        self.max_ee_r = args.max_ee_r if args.max_ee_r is not None else args.max_ee
        self.trunc_q = args.trunc_q
        self.pooling_method = args.pooling_method
        self.chimera_method = args.chimera_method
        self.min_fold_parent_over_abundance = args.min_fold_parent_over_abundance
        self.allow_one_off = args.allow_one_off
        self.n_reads_learn = args.n_reads_learn

        # Diversity analysis
        self.sampling_depth = args.sampling_depth
        self.max_rarefaction_depth = args.max_rarefaction_depth

        # Group-significance/taxonomy-collapse/classifier stats and the PDF report
        self.skip_advanced_stats = args.skip_advanced_stats
        self.skip_report = args.skip_report
        self.report_column = args.report_metadata_column

        self.fastq_list = []
        self.sample_dict = {}

        self.run()

    def run(self):
        self.fastq_list = fastq_utils.list_fastq(self.input_folder)
        self.checks()
        self.output_folder.mkdir(parents=True, exist_ok=True)

        input_folder = self.input_folder
        if self.reverse_complement:
            print('Reverse complementing reads...')
            rc_folder = self.output_folder / 'rc_reads'
            rc_folder.mkdir(parents=True, exist_ok=True)
            fastq_utils.rc_fastq_parallel(self.fastq_list, rc_folder, self.parallel)
            input_folder = rc_folder
            self.fastq_list = fastq_utils.list_fastq(rc_folder)

        # Built from self.fastq_list as it stands now (i.e. after -rc, if
        # used, has already pointed it at the reverse-complemented copies)
        # so the QA report's "Input sample files" page and provenance record
        # the files actually fed into the pipeline, not the pre-RC originals.
        self.sample_dict = fastq_utils.parse_fastq_list(self.fastq_list)

        demux_qza = self.output_folder / 'demux-seqs.qza'
        needs_fastq_roundtrip = bool(self.its1 or self.its2) or self.min_len > 0 or self.max_len > 0

        if not needs_fastq_roundtrip:
            print('Importing data into QIIME2...')
            self._import(input_folder, demux_qza)
        else:
            self._run_with_roundtrip(input_folder, demux_qza)

        qiime_wrapper.demux_summary(demux_qza, self.output_folder / 'demux-seqs.qzv')

        print('Denoising data with DADA2...')
        repseq_qza = self.output_folder / 'rep-seqs.qza'
        table_qza = self.output_folder / 'table.qza'
        stats_qza = self.output_folder / 'stats.qza'
        base_transition_stats_qza = self.output_folder / 'base-transition-stats.qza'
        if self.paired:
            qiime_wrapper.dada2_denoise_paired(
                demux_qza, repseq_qza, table_qza, stats_qza, base_transition_stats_qza,
                n_threads=self.cpu, max_ee_f=self.max_ee, max_ee_r=self.max_ee_r, trunc_q=self.trunc_q,
                pooling_method=self.pooling_method, chimera_method=self.chimera_method,
                min_fold_parent_over_abundance=self.min_fold_parent_over_abundance,
                allow_one_off=self.allow_one_off, n_reads_learn=self.n_reads_learn)
        else:
            qiime_wrapper.dada2_denoise_single(
                demux_qza, repseq_qza, table_qza, stats_qza, base_transition_stats_qza,
                n_threads=self.cpu, max_ee=self.max_ee, trunc_q=self.trunc_q,
                pooling_method=self.pooling_method, chimera_method=self.chimera_method,
                min_fold_parent_over_abundance=self.min_fold_parent_over_abundance,
                allow_one_off=self.allow_one_off, n_reads_learn=self.n_reads_learn)
        qiime_wrapper.metadata_tabulate(stats_qza, self.output_folder / 'stats.qzv')

        print('Exporting BIOM table...')
        biom_folder = self.output_folder / 'biom_table'
        qiime_wrapper.export(table_qza, biom_folder)

        print('Summarizing feature table and representative sequences...')
        qiime_wrapper.sample_summarize(self.metadata_file, table_qza, self.output_folder / 'table.qzv',
                                        self.output_folder / 'feature-frequencies.qza',
                                        self.output_folder / 'sample-frequencies.qza')
        qiime_wrapper.seq_summary(repseq_qza, self.output_folder / 'rep-seqs.qzv')

        print('Aligning representative sequences and building phylogenetic tree...')
        aligned_qza = self.output_folder / 'aligned-rep-seqs.qza'
        masked_qza = self.output_folder / 'masked-aligned-rep-seqs.qza'
        unrooted_qza = self.output_folder / 'unrooted-tree.qza'
        rooted_qza = self.output_folder / 'rooted-tree.qza'
        qiime_wrapper.phylogeny(repseq_qza, aligned_qza, masked_qza, unrooted_qza, rooted_qza)
        qiime_wrapper.export(unrooted_qza, self.output_folder)

        print('Analyzing alpha and beta diversity...')
        qiime_wrapper.core_diversity(self.cpu, self.metadata_file, rooted_qza, table_qza, self.output_folder,
                                      sampling_depth=self.sampling_depth)

        print('Creating rarefaction plot...')
        qiime_wrapper.rarefaction(self.metadata_file, rooted_qza, table_qza,
                                   self.output_folder / 'alpha-rarefaction.qzv',
                                   max_depth=self.max_rarefaction_depth)

        print('Assigning taxonomy to representative sequences...')
        taxonomy_qza = self.output_folder / 'taxonomy.qza'
        qiime_wrapper.classify(self.qiime2_classifier, repseq_qza, taxonomy_qza)
        qiime_wrapper.metadata_tabulate(taxonomy_qza, self.output_folder / 'taxonomy.qzv')

        print('Exporting taxonomy...')
        qiime_wrapper.export(taxonomy_qza, biom_folder)

        print('Incorporating taxonomy into BIOM table...')
        taxonomy_tsv = biom_folder / 'taxonomy.tsv'
        biom_utils.rewrite_taxonomy_header(taxonomy_tsv)
        feature_table_biom = biom_folder / 'feature-table.biom'
        table_with_taxonomy_biom = biom_folder / 'table-with-taxonomy.biom'
        biom_utils.add_metadata(feature_table_biom, taxonomy_tsv, table_with_taxonomy_biom)

        print('Exporting BIOM table with taxonomy...')
        biom_utils.convert_to_tsv(table_with_taxonomy_biom, taxonomy_tsv,
                                   self.output_folder / 'biom_table' / 'table-with-taxonomy.biom.tsv')

        print('Creating bar plot of sample composition...')
        qiime_wrapper.taxa_barplot(table_qza, taxonomy_qza, self.metadata_file,
                                    self.output_folder / 'taxa-bar-plots.qzv')

        print('Exporting sample read-frequency table...')
        qiime_wrapper.export(self.output_folder / 'sample-frequencies.qza',
                              self.output_folder / 'sample_frequencies')
        print('Exporting DADA2 denoising stats...')
        qiime_wrapper.export(stats_qza, self.output_folder / 'dada2_stats')

        if not self.skip_advanced_stats:
            self._run_advanced_stats(table_qza, taxonomy_qza)

        print('Writing run provenance metadata...')
        self._write_run_metadata()

        if not self.skip_report:
            print('Building PDF report...')
            report_path = report.build_report(self.output_folder, self.metadata_file, self.report_column)
            print(f'Report written to {report_path}')

        print('DONE!')

    def _run_advanced_stats(self, table_qza, taxonomy_qza):
        """Diversity group-significance tests, a genus-level composition
        table, and best-effort sample classification. Not part of the core
        pipeline the way DADA2/phylogeny/diversity are: every step here is
        either read-only reporting or explicitly skippable per metadata
        column when the data can't support it.
        """
        core_metrics_dir = self.output_folder / 'core-metrics-results'

        sample_freq_path = self.output_folder / 'sample_frequencies' / 'metadata.tsv'
        sample_frequencies = (report_data.parse_sample_frequencies(sample_freq_path)
                               if sample_freq_path.exists() else {})
        # A near-empty input sample can survive DADA2 as a zero-read row
        # (retain-all-samples defaults to True); it must not count toward
        # group eligibility for alpha/beta-group-significance/classify-samples.
        final_sample_ids = [sid for sid, freq in sample_frequencies.items() if freq > 0] \
            or list(self.sample_dict)

        if metadata_utils.has_alpha_group_significance_column(self.metadata_file, final_sample_ids):
            print('Testing alpha diversity group significance...')
            for metric in ('faith_pd', 'observed_features', 'shannon', 'evenness'):
                alpha_qza = core_metrics_dir / f'{metric}_vector.qza'
                if alpha_qza.exists():
                    qiime_wrapper.alpha_group_significance(
                        alpha_qza, self.metadata_file,
                        self.output_folder / f'alpha-group-significance-{metric}.qzv')
        else:
            print('Skipping alpha diversity group significance: no metadata column has both a '
                  'repeated and a varying value.')

        eligible_columns = metadata_utils.eligible_categorical_columns(self.metadata_file, final_sample_ids)

        print('Testing beta diversity group significance...')
        for column in eligible_columns:
            for metric in ('bray_curtis', 'unweighted_unifrac'):
                distance_qza = core_metrics_dir / f'{metric}_distance_matrix.qza'
                if distance_qza.exists():
                    qiime_wrapper.beta_group_significance(
                        distance_qza, self.metadata_file, column,
                        self.output_folder / f'beta-group-significance-{column}-{metric}.qzv')

        print('Exporting PCoA ordinations...')
        for metric in ('bray_curtis', 'unweighted_unifrac'):
            pcoa_qza = core_metrics_dir / f'{metric}_pcoa_results.qza'
            if pcoa_qza.exists():
                qiime_wrapper.export(pcoa_qza, core_metrics_dir / f'{metric}_pcoa_export')

        # Genus (level 6) if the classifier resolved that deep for enough of
        # this dataset; qiime taxa collapse otherwise fails outright if the
        # requested level exceeds every feature's actual lineage depth,
        # which classifiers commonly don't reach for reads unrelated to
        # their training set.
        taxonomy_tsv_path = self.output_folder / 'biom_table' / 'taxonomy.tsv'
        max_depth = taxonomy.max_lineage_depth(taxonomy_tsv_path) if taxonomy_tsv_path.exists() else 6
        collapse_level = min(6, max_depth) if max_depth > 0 else 1
        print(f'Collapsing feature table to taxonomic level {collapse_level}...')
        collapsed_qza = self.output_folder / 'table-genus.qza'
        relative_qza = self.output_folder / 'table-genus-relative.qza'
        qiime_wrapper.taxa_collapse(table_qza, taxonomy_qza, collapse_level, collapsed_qza)
        qiime_wrapper.relative_frequency(collapsed_qza, relative_qza)

        print('Training sample classifiers for eligible metadata columns...')
        for column in eligible_columns:
            class_sizes = metadata_utils.class_sizes(self.metadata_file, column, final_sample_ids)
            smallest_class = min(class_sizes.values())
            if smallest_class < 2:
                print(f'\tSkipping classifier for "{column}": at least one class has fewer than 2 samples.')
                continue
            classifier_dir = self.output_folder / f'sample-classifier-{column}'
            try:
                qiime_wrapper.classify_samples(table_qza, self.metadata_file, column, classifier_dir,
                                                cv=min(5, smallest_class))
            except subprocess.CalledProcessError:
                print(f'\tSkipping classifier for "{column}": training failed (see the QIIME2 error above).')

    def _write_run_metadata(self):
        """Writes run_metadata.json: who ran this, when, with what exact
        command/parameters/inputs, against which QIIME2/plugin versions --
        for QA/audit purposes and for the PDF report's provenance pages.
        Written even with --skip-report, so it's available without one.
        """
        run_metadata = provenance.build_run_metadata(
            qiime2_its_version=__version__,
            command_line=provenance.format_command_line(sys.argv),
            start_time=self.start_time,
            end_time=datetime.now().astimezone(),
            username=getpass.getuser(),
            hostname=socket.gethostname(),
            platform_string=platform.platform(),
            conda_env=self.qiime2_env,
            qiime_info_text=qiime_wrapper.qiime_info(),
            bbduk_version_text=size_filter.bbduk_version(),
            input_folder=self.input_folder,
            metadata_file=self.metadata_file,
            classifier_file=self.qiime2_classifier,
            output_folder=self.output_folder,
            sample_dict=self.sample_dict,
            parameters=self.cli_args,
        )
        provenance.write_run_metadata(self.output_folder / 'run_metadata.json', run_metadata)

    def _import(self, fastq_folder, output_qza):
        if self.paired:
            qiime_wrapper.import_fastq_pe(fastq_folder, output_qza)
        else:
            qiime_wrapper.import_fastq_se(fastq_folder, output_qza)

    def _run_with_roundtrip(self, input_folder, demux_qza):
        """ITS extraction and/or size filtering: both require dropping back to
        fastq (ITSxpress runs on a qza; BBDuk size filtering does not), so this
        path imports once, optionally trims via ITSxpress, exports, optionally
        cleans/size-filters, then re-imports the final reads as `demux_qza`.
        """
        print('Importing data into QIIME2...')
        raw_qza = self.output_folder / 'raw-demux-seqs.qza'
        self._import(input_folder, raw_qza)
        current_qza = raw_qza

        if self.its1 or self.its2:
            print('Extracting ITS region with ITSxpress...')
            its_qza = self.output_folder / 'its-demux-seqs.qza'
            if self.paired:
                itsxpress_wrapper.trim_pair_unmerged(current_qza, its_qza, self.region, self.taxa,
                                                      threads=self.cpu)
            else:
                itsxpress_wrapper.trim_single(current_qza, its_qza, self.region, self.taxa, threads=self.cpu)
            current_qza = its_qza

        print('Exporting reads for post-processing...')
        export_folder = self.output_folder / 'exported_reads'
        qiime_wrapper.export(current_qza, export_folder)
        exported_fastq = fastq_utils.list_fastq(export_folder)
        # `qiime tools export` also writes MANIFEST/metadata.yml, which the
        # Casava-only re-import format below rejects.
        fastq_utils.strip_non_fastq_files(export_folder, exported_fastq)

        if self.its1 or self.its2:
            print('Checking for empty entries...')
            if self.single:
                fastq_utils.remove_empties_se_parallel(exported_fastq, self.parallel)
            else:
                exported_sample_dict = fastq_utils.parse_fastq_list(exported_fastq)
                fastq_utils.remove_empties_pe_parallel(exported_sample_dict, self.parallel)

        reimport_folder = export_folder
        if self.min_len > 0 or self.max_len > 0:
            print('Filtering reads based on size...')
            size_folder = self.output_folder / 'size_filtered'
            size_folder.mkdir(parents=True, exist_ok=True)
            if self.single:
                size_filter.size_select_se_parallel(exported_fastq, size_folder, self.min_len, self.max_len,
                                                     self.cpu, self.parallel)
            else:
                exported_sample_dict = fastq_utils.parse_fastq_list(exported_fastq)
                size_filter.size_select_pe_parallel(exported_sample_dict, size_folder, self.min_len,
                                                     self.max_len, self.cpu, self.parallel)
            reimport_folder = size_folder

        print('Re-importing processed reads into QIIME2...')
        self._import(reimport_folder, demux_qza)

    def checks(self):
        if not (self.paired or self.single):
            raise ValueError('You must state if reads are single-end or paired-end ("-se" or "-pe").')

        if not self.fastq_list:
            raise ValueError('No fastq files found in the provided input folder.')

        fastq_utils.validate_casava_filenames(self.fastq_list)

        empty_samples = sorted({Path(fq).name.split('_')[0] for fq in self.fastq_list
                                 if fastq_utils.is_empty_fastq(fq)})
        if empty_samples:
            raise ValueError(
                'The following sample(s) have an empty (zero-read) fastq file: {}. ITSxpress and '
                'DADA2 cannot process an empty sample -- remove it from the input folder and from '
                'the metadata file before running.'.format(', '.join(empty_samples)))

        env_checks.check_qiime2_env_active()

        if self.its1 and self.its2:
            raise ValueError('You cannot choose both ITS1 and ITS2 for the same analysis.')

        if (self.its1 or self.its2) and self.taxa not in TAXA_CODES:
            raise ValueError(f'--taxa must be one of: {", ".join(sorted(TAXA_CODES))}.')

        if self.min_len > 0 or self.max_len > 0:
            env_checks.check_executable(
                'bbduk.sh', 'Install BBTools/BBMap into your QIIME2 environment, e.g. '
                            '"conda install -c bioconda -c conda-forge bbmap".')


def build_parser():
    parser = argparse.ArgumentParser(description='Run QIIME2 on ITS amplicon data using ITSxpress and DADA2')
    parser.add_argument('-q', '--qiime2', metavar='rachis-qiime2-2026.7', required=True, type=str,
                         help='Name of your QIIME2 conda environment. Mandatory.')
    parser.add_argument('-i', '--input', metavar='/input_folder/', required=True, type=str,
                         help='Input folder where the fastq reads are located. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/', required=True, type=str,
                         help='Output folder for QIIME2 files. Mandatory.')
    parser.add_argument('-m', '--metadata', metavar='qiime2_metadata.tsv', required=True, type=str,
                         help='Validated QIIME2 metadata file (samples description). Mandatory.')
    parser.add_argument('-c', '--classifier', metavar='unite_classifier_qiime2.qza', required=True, type=str,
                         help='Classifier for QIIME2. See "qiime2-its-train-unite" to compile one. Mandatory.')
    parser.add_argument('-t', '--threads', metavar='4', default=4, type=int,
                         help='Number of CPUs. Default is 4.')
    parser.add_argument('-p', '--parallel-processes', metavar='1', default=1, type=int,
                         help='Samples to process in parallel. For example, with 16 threads, using 4 '
                              'parallel processes will run 4 samples in parallel using 4 threads each. '
                              'Default is 1.')
    parser.add_argument('-rc', '--reverse_complement', action='store_true',
                         help='Use this flag if your reads are in reverse complement, for example if you '
                              'sequenced from 5.8S to 18S. Optional.')
    parser.add_argument('--min-len', type=int, default=0,
                         help='Minimum read length to keep. Default is 0 (no minimum).')
    parser.add_argument('--max-len', type=int, default=0,
                         help='Maximum read length to keep. Default is 0 (no maximum).')

    read_type = parser.add_mutually_exclusive_group(required=True)
    read_type.add_argument('-se', action='store_true', help='Reads are single-end (one fastq file per sample).')
    read_type.add_argument('-pe', action='store_true', help='Reads are paired-end (two fastq files per sample).')

    its_region = parser.add_mutually_exclusive_group()
    its_region.add_argument('--extract-its1', action='store_true',
                             help='Extract the ITS1 region with ITSxpress. Cannot be used with --extract-its2.')
    its_region.add_argument('--extract-its2', action='store_true',
                             help='Extract the ITS2 region with ITSxpress. Cannot be used with --extract-its1.')

    parser.add_argument('--taxa', metavar='Fungi', default='Fungi', choices=sorted(TAXA_CODES),
                         help='Taxon of interest for ITSxpress. One of: ' + ', '.join(sorted(TAXA_CODES)))

    dada2 = parser.add_argument_group(
        'DADA2 denoising',
        'QIIME2\'s DADA2 defaults are tuned for Illumina data. For noisier single-end platforms '
        '(e.g. IonTorrent), consider raising --max-ee and passing --allow-one-off: their higher '
        'per-base error rate means more real reads get discarded by the default max-expected-errors '
        'threshold, and homopolymer-driven indels get miscalled as one-off bimeras by the default '
        'chimera search.')
    dada2.add_argument('--max-ee', metavar='2.0', type=float, default=2.0,
                        help='Reads (forward reads, for paired-end) with more expected errors than this '
                             'are discarded. Default is 2.0.')
    dada2.add_argument('--max-ee-r', metavar='2.0', type=float, default=None,
                        help='Max expected errors for the reverse read (paired-end only). Defaults to '
                             '--max-ee.')
    dada2.add_argument('--trunc-q', metavar='2', type=int, default=2,
                        help='Truncate reads at the first quality score at or below this value. Default is 2.')
    dada2.add_argument('--pooling-method', choices=['independent', 'pseudo'], default='independent',
                        help='Sample pooling strategy for denoising. Default is independent.')
    dada2.add_argument('--chimera-method', choices=['consensus', 'none'], default='consensus',
                        help='Chimera removal method. Default is consensus.')
    dada2.add_argument('--min-fold-parent-over-abundance', metavar='1.0', type=float, default=1.0,
                        help='Minimum abundance fold-change a chimera\'s parent must have over the '
                             'candidate chimera. Default is 1.0.')
    dada2.add_argument('--allow-one-off', action='store_true',
                        help='Also flag one-mismatch/indel-from-exact bimeras as chimeric. Off by default; '
                             'consider enabling for homopolymer-heavy single-end data (e.g. IonTorrent).')
    dada2.add_argument('--n-reads-learn', metavar='1000000', type=int, default=1000000,
                        help='Minimum number of reads used to train the DADA2 error model. Default is '
                             '1000000 (fewer reads speed up small runs at some cost to error-model quality).')

    diversity = parser.add_argument_group('Diversity analysis')
    diversity.add_argument('--sampling-depth', metavar='1000', type=int, default=1000,
                            help='Rarefaction depth for core-metrics-phylogenetic. Samples with fewer '
                                 'reads than this are excluded. Default is 1000 -- lower this for small/'
                                 'pilot datasets.')
    diversity.add_argument('--max-rarefaction-depth', metavar='4000', type=int, default=4000,
                            help='Maximum depth for the alpha-rarefaction plot. Default is 4000.')

    advanced = parser.add_argument_group(
        'Advanced stats and report',
        'Diversity group-significance tests, a genus-level composition table, best-effort sample '
        'classification per eligible metadata column, and a PDF summary report. Run by default '
        'after the core pipeline.')
    advanced.add_argument('--skip-advanced-stats', action='store_true',
                           help='Skip group-significance tests, taxonomy collapse, and sample '
                                'classification.')
    advanced.add_argument('--skip-report', action='store_true',
                           help='Skip building the PDF summary report.')
    advanced.add_argument('--report-metadata-column', metavar='COLUMN', default=None, type=str,
                           help='Categorical metadata column to group the PDF report\'s alpha/beta '
                                'diversity plots by. Defaults to the first eligible (>=2 distinct '
                                'values, >=2 samples each) categorical column. Group-significance '
                                'tests and sample classification still run against every eligible '
                                'column regardless of this choice -- it only affects the report.')

    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    Pipeline(args)


if __name__ == '__main__':
    main()
