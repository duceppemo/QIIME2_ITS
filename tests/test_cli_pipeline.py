import pytest

from qiime2_its.cli.pipeline import Pipeline, build_parser

REQUIRED = ['-q', 'rachis-qiime2-2026.7', '-i', 'in/', '-o', 'out/', '-m', 'meta.tsv', '-c', 'clf.qza']


def test_requires_se_or_pe():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED)


def test_se_and_pe_are_mutually_exclusive():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED + ['-se', '-pe'])


def test_se_alone_parses(capsys):
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
    args.output = str(tmp_path / 'out')
    for key, value in overrides.items():
        setattr(args, key, value)
    mocker.patch.object(Pipeline, 'run', lambda self: None)
    pipeline = Pipeline(args)
    fq = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
    fq.write_text('x')
    pipeline.fastq_list = [fq]
    return pipeline


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

    def test_no_bbduk_check_when_no_size_filtering_requested(self, mocker, tmp_path, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        mocker.patch('qiime2_its.env_checks.shutil.which', return_value=None)
        pipeline = _make_pipeline(mocker, tmp_path)  # min_len/max_len default to 0

        pipeline.checks()  # should not raise despite bbduk.sh being "missing"
