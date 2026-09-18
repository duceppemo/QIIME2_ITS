import pytest

from qiime2_its.cli import train_ncbi


@pytest.fixture(autouse=True)
def mock_env_check(mocker):
    return mocker.patch('qiime2_its.cli.train_ncbi.env_checks.check_qiime2_env_active')


@pytest.fixture
def mock_pipeline(mocker, tmp_path):
    """Mocks every external call run() makes, so it can be exercised end to
    end without hitting NCBI, the real taxdump, or a real QIIME2 install."""
    mocker.patch('qiime2_its.cli.train_ncbi.downloader.download')
    mocker.patch('qiime2_its.cli.train_ncbi.downloader.extract_targz')
    mocker.patch('qiime2_its.cli.train_ncbi.download_sequences',
                 side_effect=lambda query, seq_file, email, api_key: seq_file.write_text('>ACC1\nACGT\n'))
    mocker.patch('qiime2_its.cli.train_ncbi.taxonomy.extract_accessions_from_fasta', return_value={'ACC1': '1001'})
    mocker.patch('qiime2_its.cli.train_ncbi.taxonomy.parse_accession2taxid', return_value={})
    mocker.patch('qiime2_its.cli.train_ncbi.taxonomy.accessions_to_taxids', return_value=[])
    mocker.patch('qiime2_its.cli.train_ncbi.taxonomy.parse_id_table', return_value={'1001': 'ACC1'})
    mocker.patch('qiime2_its.cli.train_ncbi.taxonomy.write_taxonomy_file')
    return {
        'import_seq': mocker.patch('qiime2_its.cli.train_ncbi.qiime_wrapper.import_sequences'),
        'import_taxo': mocker.patch('qiime2_its.cli.train_ncbi.qiime_wrapper.import_taxonomy'),
        'train': mocker.patch('qiime2_its.cli.train_ncbi.qiime_wrapper.train_naive_bayes_classifier'),
    }


class TestRun:
    def test_raises_if_qiime2_env_not_active(self, tmp_path, mock_env_check):
        """Regression test: run() used to never check this at all, unlike
        pipeline.py/train_unite.py."""
        mock_env_check.side_effect = EnvironmentError('You must activate your QIIME2 conda environment...')
        with pytest.raises(EnvironmentError, match='QIIME2 conda environment'):
            train_ncbi.run('some query', tmp_path / 'out', 4, 'a@b.com', None, None, None, None)

    def test_raises_on_empty_query(self, tmp_path):
        with pytest.raises(ValueError, match='query'):
            train_ncbi.run('', tmp_path / 'out', 4, 'a@b.com', None, None, None, None)

    def test_full_run_orchestrates_expected_calls(self, tmp_path, mock_pipeline):
        output_folder = tmp_path / 'out'

        train_ncbi.run('some query', output_folder, 4, 'a@b.com', None,
                        tmp_path / 'taxdump.tar.gz', tmp_path / 'acc2taxid.gz', tmp_path / 'dead.gz')

        mock_pipeline['import_seq'].assert_called_once()
        mock_pipeline['import_taxo'].assert_called_once()
        mock_pipeline['train'].assert_called_once()
        classifier_path = mock_pipeline['train'].call_args.args[2]
        assert str(classifier_path).endswith('seq_ncbi.qza')

    def test_reuses_existing_seq_fasta_with_a_warning(self, tmp_path, mock_pipeline, mocker, capsys):
        """Regression test: an existing seq.fasta was silently reused with
        no indication to the user -- dangerous for a crashed/partial prior
        download, or a stale file from a different query reusing the same
        output folder."""
        output_folder = tmp_path / 'out'
        output_folder.mkdir()
        stale = output_folder / 'seq.fasta'
        stale.write_text('>OLD\nTTTT\n')
        # mock_pipeline already patches download_sequences to write a fresh
        # file; re-patch it here bare so we can assert it's never called.
        mock_download = mocker.patch('qiime2_its.cli.train_ncbi.download_sequences')

        train_ncbi.run('some query', output_folder, 4, 'a@b.com', None,
                        tmp_path / 'taxdump.tar.gz', tmp_path / 'acc2taxid.gz', tmp_path / 'dead.gz')

        mock_download.assert_not_called()
        assert 'seq.fasta' in capsys.readouterr().out
        assert stale.read_text() == '>OLD\nTTTT\n'  # untouched, not silently overwritten either
