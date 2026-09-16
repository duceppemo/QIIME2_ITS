import pytest

from qiime2_its import size_filter


@pytest.fixture
def mock_run(mocker):
    return mocker.patch('qiime2_its.size_filter.subprocess.run')


def test_min_and_max_len_both_included(mock_run):
    size_filter.size_select_se('r1.fastq.gz', 'out', min_len=100, max_len=300, threads=4)
    cmd = mock_run.call_args.args[0]
    assert 'minlength=100' in cmd
    assert 'maxlength=300' in cmd


def test_max_len_zero_omits_maxlength_flag(mock_run):
    """Regression test: the original code deleted --maxlength by a hardcoded
    list index, which broke if the argument order ever changed."""
    size_filter.size_select_se('r1.fastq.gz', 'out', min_len=100, max_len=0, threads=4)
    cmd = mock_run.call_args.args[0]
    assert 'minlength=100' in cmd
    assert not any(arg.startswith('maxlength=') for arg in cmd)


def test_min_len_zero_omits_minlength_flag(mock_run):
    size_filter.size_select_se('r1.fastq.gz', 'out', min_len=0, max_len=300, threads=4)
    cmd = mock_run.call_args.args[0]
    assert not any(arg.startswith('minlength=') for arg in cmd)
    assert 'maxlength=300' in cmd


def test_size_select_pe_includes_both_mates(mock_run):
    size_filter.size_select_pe('r1.fastq.gz', 'r2.fastq.gz', 'out', min_len=50, max_len=0, threads=2)
    cmd = mock_run.call_args.args[0]
    assert 'in=r1.fastq.gz' in cmd
    assert 'in2=r2.fastq.gz' in cmd
    assert any(arg.startswith('out2=') for arg in cmd)


def test_bbduk_version_returns_stderr(mock_run):
    mock_run.return_value.stderr = 'BBTools version 39.80\n'
    assert size_filter.bbduk_version() == 'BBTools version 39.80\n'
    mock_run.assert_called_once_with(['bbduk.sh', '--version'], capture_output=True, text=True)


def test_bbduk_version_none_when_not_installed(mocker):
    mocker.patch('qiime2_its.size_filter.subprocess.run', side_effect=FileNotFoundError)
    assert size_filter.bbduk_version() is None
