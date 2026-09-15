import pytest

from qiime2_its import itsxpress_wrapper


@pytest.fixture
def mock_run(mocker):
    return mocker.patch('qiime2_its.itsxpress_wrapper.subprocess.run')


def test_taxa_codes_cover_fungi_and_all():
    assert itsxpress_wrapper.TAXA_CODES['Fungi'] == 'F'
    assert itsxpress_wrapper.TAXA_CODES['All'] == 'ALL'


def test_trim_single_command(mock_run):
    itsxpress_wrapper.trim_single('raw.qza', 'trimmed.qza', 'ITS1', 'Fungi', threads=4, cluster_id=0.99)
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'itsxpress', 'trim-single',
                    '--i-per-sample-sequences', 'raw.qza',
                    '--p-region', 'ITS1', '--p-taxa', 'F',
                    '--p-threads', '4', '--p-cluster-id', '0.99',
                    '--o-trimmed', 'trimmed.qza']
    assert mock_run.call_args.kwargs.get('check') is True


def test_trim_pair_unmerged_command(mock_run):
    itsxpress_wrapper.trim_pair_unmerged('raw.qza', 'trimmed.qza', 'ITS2', 'Oomycota', threads=8)
    cmd = mock_run.call_args.args[0]
    assert cmd[:3] == ['qiime', 'itsxpress', 'trim-pair-output-unmerged']
    assert '--p-taxa' in cmd and cmd[cmd.index('--p-taxa') + 1] == 'O'
    assert '--p-threads' in cmd and cmd[cmd.index('--p-threads') + 1] == '8'
