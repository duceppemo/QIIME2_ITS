import pytest

from qiime2_its.cli import train_fasta


@pytest.fixture(autouse=True)
def mock_env_check(mocker):
    return mocker.patch('qiime2_its.cli.train_fasta.env_checks.check_qiime2_env_active')


def test_raises_if_fasta_missing(tmp_path):
    id_table = tmp_path / 'ids.tsv'
    id_table.write_text('ACC1\t1001\n')
    with pytest.raises(ValueError, match='query'):
        train_fasta.run(tmp_path / 'missing.fasta', id_table, tmp_path / 'out', taxdump=None)


def test_raises_if_id_table_missing(tmp_path):
    fasta = tmp_path / 'seqs.fasta'
    fasta.write_text('>ACC1\nACGT\n')
    with pytest.raises(ValueError, match='id-table'):
        train_fasta.run(fasta, tmp_path / 'missing.tsv', tmp_path / 'out', taxdump=None)


def test_raises_if_qiime2_env_not_active(tmp_path, mock_env_check):
    """Regression test: run() used to never check this at all, unlike
    pipeline.py/train_unite.py -- a raw, confusing subprocess failure deep
    inside qiime_wrapper instead of a clear, immediate message."""
    mock_env_check.side_effect = EnvironmentError('You must activate your QIIME2 conda environment...')
    fasta = tmp_path / 'seqs.fasta'
    fasta.write_text('>ACC1\nACGT\n')
    id_table = tmp_path / 'ids.tsv'
    id_table.write_text('ACC1\t1001\n')
    with pytest.raises(EnvironmentError, match='QIIME2 conda environment'):
        train_fasta.run(fasta, id_table, tmp_path / 'out', taxdump=None)


def test_raises_up_front_when_a_fasta_id_is_missing_from_the_id_table(tmp_path, mocker):
    """Every sequence needs a taxonomy line. QIIME2 doesn't complain about a
    mismatch -- it quietly trains on the sequences that do have one -- so
    without this check a typo in the table silently shrinks the classifier."""
    fasta = tmp_path / 'seqs.fasta'
    fasta.write_text('>ACC1 desc\nACGT\n>ACC2.1\nACGT\n')
    id_table = tmp_path / 'ids.tsv'
    id_table.write_text('ACC1\t1001\nACC2\t1002\n')  # ACC2, not ACC2.1
    mock_download = mocker.patch('qiime2_its.cli.train_fasta.downloader.download')

    with pytest.raises(ValueError, match=r'ACC2\.1'):
        train_fasta.run(fasta, id_table, tmp_path / 'out', taxdump=None)

    mock_download.assert_not_called()


def test_raises_on_a_fasta_with_no_sequences(tmp_path):
    fasta = tmp_path / 'seqs.fasta'
    fasta.write_text('')
    id_table = tmp_path / 'ids.tsv'
    id_table.write_text('ACC1\t1001\n')
    with pytest.raises(ValueError, match='No sequences'):
        train_fasta.run(fasta, id_table, tmp_path / 'out', taxdump=None)


def test_full_run_orchestrates_expected_calls(tmp_path, mocker):
    fasta = tmp_path / 'seqs.fasta'
    fasta.write_text('>ACC1\nACGT\n')
    id_table = tmp_path / 'ids.tsv'
    id_table.write_text('ACC1\t1001\n')
    output_folder = tmp_path / 'out'

    taxdump_path = tmp_path / 'taxdump.tar.gz'
    taxdump_path.write_bytes(b'fake')

    mock_extract = mocker.patch('qiime2_its.cli.train_fasta.downloader.extract_targz')
    mock_write_taxo = mocker.patch('qiime2_its.cli.train_fasta.taxonomy.write_taxonomy_file')
    mock_import_seq = mocker.patch('qiime2_its.cli.train_fasta.qiime_wrapper.import_sequences')
    mock_import_taxo = mocker.patch('qiime2_its.cli.train_fasta.qiime_wrapper.import_taxonomy')
    mock_train = mocker.patch('qiime2_its.cli.train_fasta.qiime_wrapper.train_naive_bayes_classifier')

    train_fasta.run(fasta, id_table, output_folder, taxdump=taxdump_path)

    mock_extract.assert_called_once_with(taxdump_path, output_folder)

    taxonomy_file = output_folder / 'taxonomy.txt'
    mock_write_taxo.assert_called_once_with(
        {'ACC1': '1001'}, taxonomy_file,
        output_folder / 'nodes.dmp', output_folder / 'names.dmp', output_folder / 'merged.dmp')

    # Regression coverage: import_sequences/import_taxonomy/
    # train_naive_bayes_classifier must each get the right file for the
    # right argument (the raw fasta vs. the written taxonomy.txt, and their
    # imported .qza counterparts) -- assert_called_once() alone doesn't
    # catch e.g. import_taxonomy being handed the sequence fasta instead of
    # taxonomy.txt.
    qiime2_seq = output_folder / 'seqs.qza'
    qiime2_taxo = output_folder / 'taxonomy.txt.qza'
    mock_import_seq.assert_called_once_with(fasta, qiime2_seq)
    mock_import_taxo.assert_called_once_with(taxonomy_file, qiime2_taxo)
    mock_train.assert_called_once_with(
        qiime2_seq, qiime2_taxo, output_folder / 'naive-bayes_classifier.qza', verbose=True)
