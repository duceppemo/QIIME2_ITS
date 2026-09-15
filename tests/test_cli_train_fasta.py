import pytest

from qiime2_its.cli import train_fasta


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
    mock_write_taxo.assert_called_once()
    mock_import_seq.assert_called_once()
    mock_import_taxo.assert_called_once()
    mock_train.assert_called_once()
    assert mock_train.call_args.kwargs.get('verbose') is True
