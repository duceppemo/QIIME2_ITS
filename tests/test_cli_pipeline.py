import gzip
from pathlib import Path

import pytest

from qiime2_its.cli.pipeline import Pipeline, _safe_filename_component, build_parser

REQUIRED = ['-q', 'rachis-qiime2-2026.7', '-i', 'in/', '-o', 'out/', '-m', 'meta.tsv', '-c', 'clf.qza']


def test_requires_se_or_pe():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED)


def test_se_and_pe_are_mutually_exclusive():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED + ['-se', '-pe'])


def test_se_alone_parses():
    parser = build_parser()
    args = parser.parse_args(REQUIRED + ['-se'])
    assert args.se is True
    assert args.pe is False


def test_extract_its1_and_its2_are_mutually_exclusive():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED + ['-se', '--extract-its1', '--extract-its2'])


def test_defaults():
    parser = build_parser()
    args = parser.parse_args(REQUIRED + ['-se'])
    assert args.threads == 4
    assert args.parallel_processes == 1
    assert args.min_len == 0
    assert args.max_len == 0
    assert args.taxa == 'Fungi'
    assert args.reverse_complement is False
    assert args.max_ee == 2.0
    assert args.max_ee_r is None
    assert args.allow_one_off is False
    assert args.pooling_method == 'independent'
    assert args.chimera_method == 'consensus'
    assert args.sampling_depth == 1000
    assert args.max_rarefaction_depth == 4000
    assert args.skip_advanced_stats is False
    assert args.skip_report is False
    assert args.report_metadata_column is None


def test_advanced_stats_and_report_flags_parse():
    parser = build_parser()
    args = parser.parse_args(REQUIRED + ['-se', '--skip-advanced-stats', '--skip-report',
                                          '--report-metadata-column', 'site'])
    assert args.skip_advanced_stats is True
    assert args.skip_report is True
    assert args.report_metadata_column == 'site'


def test_dada2_overrides_parse():
    parser = build_parser()
    args = parser.parse_args(REQUIRED + ['-se', '--max-ee', '5', '--allow-one-off',
                                          '--pooling-method', 'pseudo', '--sampling-depth', '50'])
    assert args.max_ee == 5.0
    assert args.allow_one_off is True
    assert args.pooling_method == 'pseudo'
    assert args.sampling_depth == 50


def test_missing_required_argument_exits():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(['-se'])


def test_invalid_taxa_choice_exits():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED + ['-se', '--taxa', 'NotARealTaxon'])


def _make_pipeline(mocker, tmp_path, **overrides):
    """Build a Pipeline with run() stubbed out, so checks() can be exercised directly."""
    parser = build_parser()
    argv = REQUIRED + ['-se']
    args = parser.parse_args(argv)
    args.input = str(tmp_path)
    args.output = str(tmp_path.parent / (tmp_path.name + '_out'))
    # Real files: checks() verifies both exist and that every fastq sample
    # has a metadata row.
    metadata = tmp_path / 'meta.tsv'
    metadata.write_text('sample-id\tsite\nsample\tA\nsiteC-rep3\tA\nsiteD-rep9\tB\n')
    classifier = tmp_path / 'clf.qza'
    classifier.write_bytes(b'')
    args.metadata = str(metadata)
    args.classifier = str(classifier)
    for key, value in overrides.items():
        setattr(args, key, value)
    mocker.patch.object(Pipeline, 'run', lambda self: None)
    pipeline = Pipeline(args)
    fq = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
    with gzip.open(fq, 'wt') as f:
        f.write('@r1\nACGT\n+\nIIII\n')
    pipeline.fastq_list = [fq]
    return pipeline


class TestSafeFilenameComponent:
    def test_leaves_ordinary_column_names_unchanged(self):
        assert _safe_filename_component('host-plant') == 'host-plant'

    def test_replaces_path_separators(self):
        """Regression test: QIIME2 doesn't restrict metadata column names
        from containing '/' (e.g. a real column named "site/plot"), and
        _run_advanced_stats builds output paths directly from the column
        name -- unsanitized, '/' would be interpreted as a subdirectory
        separator instead of a literal character in the filename."""
        assert _safe_filename_component('site/plot') == 'site_plot'
        assert _safe_filename_component('../../etc/passwd') == '.._.._etc_passwd'

    def test_replaces_other_filesystem_significant_characters(self):
        assert _safe_filename_component('treatment (mg/L)') == 'treatment__mg_L_'


class TestPipelineChecks:
    def test_min_len_requires_bbduk_on_path(self, mocker, tmp_path, monkeypatch):
        """Regression test: --min-len/--max-len depend on bbduk.sh, which the
        QIIME2 environment does not install by default -- this must fail
        early with a clear message, not deep inside a subprocess call."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        mocker.patch('qiime2_its.env_checks.shutil.which', return_value=None)
        pipeline = _make_pipeline(mocker, tmp_path, min_len=100)

        with pytest.raises(EnvironmentError, match='bbduk.sh'):
            pipeline.checks()

    def test_min_len_rejects_an_equals_sign_in_the_output_path(self, mocker, tmp_path, monkeypatch):
        """Regression test: BBDuk parses its own arguments as key=value
        pairs -- confirmed empirically against the real bbduk.sh (39.80) --
        so an output path containing "=" (e.g. a date-stamped folder like
        "run=2026-09-18") truncates the out= argument at the "=" and BBDuk
        dies with a cryptic "Can't read file" error deep inside a
        ThreadPoolExecutor worker instead of a clear message here."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        mocker.patch('qiime2_its.env_checks.shutil.which', return_value='/usr/bin/bbduk.sh')
        pipeline = _make_pipeline(mocker, tmp_path, min_len=100, output=str(tmp_path.parent / 'run=2' / 'out'))

        with pytest.raises(ValueError, match='run=2'):
            pipeline.checks()

    def test_no_bbduk_check_when_no_size_filtering_requested(self, mocker, tmp_path, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        mocker.patch('qiime2_its.env_checks.shutil.which', return_value=None)
        pipeline = _make_pipeline(mocker, tmp_path)  # min_len/max_len default to 0

        pipeline.checks()  # should not raise despite bbduk.sh being "missing"

    def test_warns_when_declared_env_does_not_match_the_active_one(self, mocker, tmp_path, monkeypatch, capsys):
        """Regression test: -q/--qiime2 was accepted and recorded (for
        provenance) but never actually cross-checked against the real active
        environment -- so declaring one env while a different one was
        actually active went unnoticed, same class of bug already fixed for
        train_unite.py's -q."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path, qiime2='some-other-env')

        pipeline.checks()

        assert 'some-other-env' in capsys.readouterr().out

    def test_no_warning_when_declared_env_matches_the_active_one(self, mocker, tmp_path, monkeypatch, capsys):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path)  # qiime2 defaults to rachis-qiime2-2026.7 (REQUIRED)

        pipeline.checks()

        assert 'Warning' not in capsys.readouterr().out

    def test_empty_sample_rejected_with_clear_message(self, mocker, tmp_path, monkeypatch):
        """Regression test: a zero-read sample previously wasn't caught until
        several steps into the pipeline, where ITSxpress's HMM search fails
        with a cryptic external-tool stack trace instead of a clear message
        identifying which sample is the problem."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path)
        empty_fq = tmp_path / 'siteC-rep3_S7_L001_R1_001.fastq.gz'
        with gzip.open(empty_fq, 'wt'):
            pass
        pipeline.fastq_list.append(empty_fq)

        with pytest.raises(ValueError, match='siteC-rep3'):
            pipeline.checks()

    def test_unpaired_sample_rejected_with_clear_message(self, mocker, tmp_path, monkeypatch):
        """Regression test: a paired-end sample missing its R2 mate used to
        go undetected here and crash later with a bare IndexError inside a
        ThreadPoolExecutor worker (remove_empties_pe_parallel /
        size_select_pe_parallel both index reads[1] with no length check)."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path, se=False, pe=True)
        # Replace the fixture's generic R1-only "sample" with a distinctively
        # named one, so the assertion below can verify the *offending sample
        # id* made it into the message -- not just that the word "sample"
        # (present in the message's own generic wording, "sample(s)") shows
        # up somewhere, which would pass even if the id list were dropped
        # from the message entirely.
        r1 = tmp_path / 'siteD-rep9_S9_L001_R1_001.fastq.gz'
        with gzip.open(r1, 'wt') as f:
            f.write('@r1\nACGT\n+\nIIII\n')
        pipeline.fastq_list = [r1]

        with pytest.raises(ValueError, match='siteD-rep9'):
            pipeline.checks()

    def test_complete_pairs_are_not_rejected(self, mocker, tmp_path, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path, se=False, pe=True)
        r2 = tmp_path / 'sample_bc_L001_R2_001.fastq.gz'
        with gzip.open(r2, 'wt') as f:
            f.write('@r1\nACGT\n+\nIIII\n')
        pipeline.fastq_list.append(r2)

        pipeline.checks()  # should not raise: every sample has both mates

    def test_pairing_check_only_applies_in_paired_end_mode(self, mocker, tmp_path, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path)  # se=True (default), single R1-only fastq

        pipeline.checks()  # should not raise: single-end mode doesn't require R2s


    @pytest.mark.parametrize('argument, label', [('metadata', 'metadata'), ('classifier', 'classifier')])
    def test_missing_metadata_or_classifier_file_rejected_up_front(self, mocker, tmp_path, monkeypatch,
                                                                   argument, label):
        """Neither file is read until after DADA2 (metadata) or after
        phylogeny + diversity (classifier) -- a typo in the path used to
        surface only hours into a real run."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path, **{argument: str(tmp_path / 'typo')})

        with pytest.raises(ValueError, match=f'{label} file does not exist'):
            pipeline.checks()

    def test_output_folder_inside_the_input_folder_rejected(self, mocker, tmp_path, monkeypatch):
        """QIIME2's Casava importer refuses an input folder containing any
        subfolder, and a nested output folder's exported_reads/ would be
        picked up as input on a re-run."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path, output=str(tmp_path / 'qiime2_out'))

        with pytest.raises(ValueError, match='must not be the input folder or inside'):
            pipeline.checks()

    def test_output_folder_equal_to_the_input_folder_rejected(self, mocker, tmp_path, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path, output=str(tmp_path))

        with pytest.raises(ValueError, match='must not be the input folder or inside'):
            pipeline.checks()

    def test_symlinked_fastq_files_are_not_mistaken_for_subfolder_files(self, mocker, tmp_path, monkeypatch):
        """Regression test (0.3.1): the subfolder check resolved each *file*,
        following its symlink, so a flat input folder of symlinks to fastq
        files stored elsewhere -- a common layout QIIME2 imports fine, and
        what the error message itself recommends -- was rejected."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path)
        storage = tmp_path.parent / (tmp_path.name + '_storage')
        storage.mkdir()
        real = storage / 'siteC-rep3_S7_L001_R1_001.fastq.gz'
        with gzip.open(real, 'wt') as f:
            f.write('@r1\nACGT\n+\nIIII\n')
        link = tmp_path / real.name
        link.symlink_to(real)
        pipeline.fastq_list.append(link)

        pipeline.checks()  # should not raise

    def test_columns_colliding_once_sanitized_for_file_names_rejected(self, mocker, tmp_path, monkeypatch):
        """ "Host Plant" and "Host_Plant" would both write
        beta-group-significance-Host_Plant-*.qzv, overwriting each other."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path)
        Path(pipeline.metadata_file).write_text('sample-id\tHost Plant\tHost_Plant\nsample\tA\tB\n')

        with pytest.raises(ValueError, match='"Host Plant" / "Host_Plant"'):
            pipeline.checks()

    def test_fastq_in_a_subfolder_rejected_unless_reverse_complementing(self, mocker, tmp_path, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path)
        nested = tmp_path / 'run1' / 'siteC-rep3_S7_L001_R1_001.fastq.gz'
        nested.parent.mkdir()
        with gzip.open(nested, 'wt') as f:
            f.write('@r1\nACGT\n+\nIIII\n')
        pipeline.fastq_list.append(nested)

        with pytest.raises(ValueError, match='run1'):
            pipeline.checks()

        pipeline.reverse_complement = True
        pipeline.checks()  # -rc flattens every file into rc_reads/ first

    def test_sample_missing_from_metadata_rejected_up_front(self, mocker, tmp_path, monkeypatch):
        """QIIME2 refuses a feature table with sample IDs absent from the
        metadata, but only at `feature-table summarize` -- after DADA2."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path)
        extra = tmp_path / 'notInMetadata_S3_L001_R1_001.fastq.gz'
        with gzip.open(extra, 'wt') as f:
            f.write('@r1\nACGT\n+\nIIII\n')
        pipeline.fastq_list.append(extra)

        with pytest.raises(ValueError, match='notInMetadata'):
            pipeline.checks()


class TestAdvancedStatsAreBestEffort:
    def test_a_failing_beta_group_significance_column_does_not_abort_the_run(self, mocker, tmp_path, monkeypatch,
                                                                                 capsys):
        """Eligibility is computed from the metadata file; QIIME2 decides
        for itself (e.g. after dropping samples below the sampling depth).
        One refused column used to kill a run that was already past DADA2."""
        import subprocess
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        pipeline = _make_pipeline(mocker, tmp_path)
        pipeline.output_folder.mkdir(parents=True)
        core = pipeline.output_folder / 'core-metrics-results'
        core.mkdir()
        (core / 'bray_curtis_distance_matrix.qza').write_bytes(b'')
        Path(pipeline.metadata_file).write_text(
            'sample-id\tsite\ns1\tA\ns2\tA\ns3\tB\ns4\tB\n')
        pipeline.sample_dict = {f's{i}': [] for i in range(1, 5)}
        wrapper = mocker.patch('qiime2_its.cli.pipeline.qiime_wrapper')
        wrapper.beta_group_significance.side_effect = subprocess.CalledProcessError(1, 'qiime')

        pipeline._run_advanced_stats('table.qza', 'taxonomy.qza')

        assert 'Skipping bray_curtis group significance for "site"' in capsys.readouterr().out
        wrapper.taxa_collapse.assert_called_once()  # carried on past the failure


class TestThreadsAreHonored:
    def test_phylogeny_and_classification_use_the_requested_thread_count(self, mocker, tmp_path, monkeypatch):
        """Regression test: both steps used to take every core on the
        machine ('auto' / n_jobs=0) regardless of -t/--threads."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        parser = build_parser()
        fq = tmp_path / 'in' / 'sample_bc_L001_R1_001.fastq.gz'
        fq.parent.mkdir()
        with gzip.open(fq, 'wt') as f:
            f.write('@r1\nACGT\n+\nIIII\n')
        metadata = tmp_path / 'meta.tsv'
        metadata.write_text('sample-id\tsite\nsample\tA\n')
        classifier = tmp_path / 'clf.qza'
        classifier.write_bytes(b'')
        args = parser.parse_args(['-q', 'rachis-qiime2-2026.7', '-i', str(fq.parent), '-o', str(tmp_path / 'out'),
                                  '-m', str(metadata), '-c', str(classifier), '-se', '-t', '2',
                                  '--skip-advanced-stats', '--skip-report'])
        mocker.patch('qiime2_its.cli.pipeline.env_checks.clamp_cpu', return_value=2)
        wrapper = mocker.patch('qiime2_its.cli.pipeline.qiime_wrapper')
        mocker.patch('qiime2_its.cli.pipeline.biom_utils')
        mocker.patch.object(Pipeline, '_write_run_metadata')

        Pipeline(args)

        assert wrapper.phylogeny.call_args.kwargs['n_threads'] == 2
        assert wrapper.classify.call_args.kwargs['n_jobs'] == 2


class TestWriteRunMetadata:
    def test_records_the_actual_active_env_not_the_declared_one(self, mocker, tmp_path, monkeypatch):
        """Regression test: checks() already warns when -q/--qiime2 doesn't
        match the real active environment, but _write_run_metadata() used to
        record self.qiime2_env (the possibly-stale declared value) in
        run_metadata.json regardless -- defeating the point of a
        provenance/audit record."""
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2027.1')
        pipeline = _make_pipeline(mocker, tmp_path, qiime2='rachis-qiime2-2026.7')
        pipeline.checks()  # sets self.active_env

        mocker.patch('qiime2_its.cli.pipeline.qiime_wrapper.qiime_info', return_value='')
        mocker.patch('qiime2_its.cli.pipeline.size_filter.bbduk_version', return_value=None)
        mocker.patch('qiime2_its.cli.pipeline.provenance.write_run_metadata')
        build_metadata = mocker.patch('qiime2_its.cli.pipeline.provenance.build_run_metadata')

        pipeline._write_run_metadata()

        assert build_metadata.call_args.kwargs['conda_env'] == 'rachis-qiime2-2027.1'
