from pathlib import Path

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


class TestQueryFingerprint:
    def test_same_search_string_gives_the_same_fingerprint(self):
        assert train_ncbi._query_fingerprint('txid4762[Organism:exp]') == \
            train_ncbi._query_fingerprint('txid4762[Organism:exp]')

    def test_different_search_strings_give_different_fingerprints(self):
        assert train_ncbi._query_fingerprint('query A') != train_ncbi._query_fingerprint('query B')

    def test_accession_list_is_fingerprinted_by_content_not_path(self, tmp_path):
        """Editing the accession-list file (same path, different content)
        must count as a different query -- otherwise re-pointing -q at an
        updated list in place would still incorrectly reuse the old
        seq.fasta."""
        acc_list = tmp_path / 'accessions.list'
        acc_list.write_text('ACC1\nACC2\n')
        first = train_ncbi._query_fingerprint(str(acc_list))

        acc_list.write_text('ACC1\nACC2\nACC3\n')
        second = train_ncbi._query_fingerprint(str(acc_list))

        assert first != second


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

        seq_file = output_folder / 'seq.fasta'
        qiime2_seq = output_folder / 'seq.fasta.qza'
        taxonomy_file = output_folder / 'taxonomy.txt'
        qiime2_taxo = output_folder / 'taxonomy.txt.qza'
        classifier_file = output_folder / 'seq_ncbi.qza'
        # Regression coverage: each import call must get its own file (the
        # downloaded seq.fasta vs. the written taxonomy.txt), not each
        # other's -- assert_called_once() alone can't catch that mixup.
        mock_pipeline['import_seq'].assert_called_once_with(seq_file, qiime2_seq)
        mock_pipeline['import_taxo'].assert_called_once_with(taxonomy_file, qiime2_taxo)
        mock_pipeline['train'].assert_called_once_with(qiime2_seq, qiime2_taxo, classifier_file)

    def test_live_and_dead_accession2taxid_are_merged_not_overwritten(self, tmp_path, mock_pipeline, mocker):
        """Regression-shaped coverage: dead_nucl.accession2taxid.gz entries
        must be merged into (not replace) the ones already found in
        nucl_gb.accession2taxid.gz -- accessions_to_taxids() needs to see
        both, not just whichever accession2taxid file was parsed last."""
        def _parse_accession2taxid(path, acc_dict):
            # Path(path).name, not str(path): tmp_path's own directory is
            # derived from this test's name (which contains "dead"), so a
            # substring check against the full path would false-match the
            # *live* acc2taxid.gz too.
            return {'ACC2': '2002'} if Path(path).name.startswith('dead') else {'ACC1': '1001'}

        mocker.patch('qiime2_its.cli.train_ncbi.taxonomy.parse_accession2taxid',
                      side_effect=_parse_accession2taxid)
        accessions_to_taxids = mocker.patch('qiime2_its.cli.train_ncbi.taxonomy.accessions_to_taxids',
                                             return_value=[])
        output_folder = tmp_path / 'out'

        train_ncbi.run('some query', output_folder, 4, 'a@b.com', None,
                        tmp_path / 'taxdump.tar.gz', tmp_path / 'acc2taxid.gz', tmp_path / 'dead.gz')

        acc2taxid_dict = accessions_to_taxids.call_args.args[1]
        assert acc2taxid_dict == {'ACC1': '1001', 'ACC2': '2002'}

    def test_reuses_existing_seq_fasta_when_the_query_hash_matches(self, tmp_path, mock_pipeline, mocker, capsys):
        """A prior run for the *same* query recorded its fingerprint
        alongside seq.fasta -- reuse it without re-downloading."""
        output_folder = tmp_path / 'out'
        output_folder.mkdir()
        stale = output_folder / 'seq.fasta'
        stale.write_text('>OLD\nTTTT\n')
        (output_folder / 'seq.fasta.query_hash').write_text(
            train_ncbi._query_fingerprint('some query') + '\n')
        mock_download = mocker.patch('qiime2_its.cli.train_ncbi.download_sequences')

        train_ncbi.run('some query', output_folder, 4, 'a@b.com', None,
                        tmp_path / 'taxdump.tar.gz', tmp_path / 'acc2taxid.gz', tmp_path / 'dead.gz')

        mock_download.assert_not_called()
        assert 'matches the current query' in capsys.readouterr().out
        assert stale.read_text() == '>OLD\nTTTT\n'  # untouched, not silently overwritten either

    def test_redownloads_when_seq_fasta_exists_but_hash_does_not_match(self, tmp_path, mock_pipeline, capsys):
        """Regression test: an existing seq.fasta used to be reused
        unconditionally -- dangerous for a stale file from a *different*
        query reusing the same output folder. No recorded hash at all (an
        output folder from before this existed, or a crashed run that never
        got to write one) is treated the same way: don't trust it."""
        output_folder = tmp_path / 'out'
        output_folder.mkdir()
        stale = output_folder / 'seq.fasta'
        stale.write_text('>OLD\nTTTT\n')
        (output_folder / 'seq.fasta.query_hash').write_text(train_ncbi._query_fingerprint('a different query') + '\n')

        train_ncbi.run('some query', output_folder, 4, 'a@b.com', None,
                        tmp_path / 'taxdump.tar.gz', tmp_path / 'acc2taxid.gz', tmp_path / 'dead.gz')

        mock_pipeline['import_seq'].assert_called_once()  # ran to completion on the fresh download
        assert 'different query' in capsys.readouterr().out
        assert stale.read_text() == '>ACC1\nACGT\n'  # overwritten by mock_pipeline's download_sequences

    def test_records_query_hash_after_a_fresh_download(self, tmp_path, mock_pipeline):
        output_folder = tmp_path / 'out'

        train_ncbi.run('some query', output_folder, 4, 'a@b.com', None,
                        tmp_path / 'taxdump.tar.gz', tmp_path / 'acc2taxid.gz', tmp_path / 'dead.gz')

        recorded = (output_folder / 'seq.fasta.query_hash').read_text().strip()
        assert recorded == train_ncbi._query_fingerprint('some query')
