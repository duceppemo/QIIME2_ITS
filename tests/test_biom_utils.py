import pytest

from qiime2_its import biom_utils


@pytest.fixture
def mock_run(mocker):
    return mocker.patch('qiime2_its.biom_utils.subprocess.run')


def test_rewrite_taxonomy_header(tmp_path):
    taxo = tmp_path / 'taxonomy.tsv'
    taxo.write_text('Feature ID\tTaxon\tConfidence\n'
                     'ASV1\tk__Fungi;...\t0.99\n'
                     'ASV2\tk__Fungi;...\t0.87\n')

    biom_utils.rewrite_taxonomy_header(taxo)

    lines = taxo.read_text().splitlines()
    assert lines[0] == '#OTUID\ttaxonomy\tconfidence'
    assert lines[1] == 'ASV1\tk__Fungi;...\t0.99'
    assert lines[2] == 'ASV2\tk__Fungi;...\t0.87'
    assert len(lines) == 3


def test_add_metadata_command(mock_run):
    biom_utils.add_metadata('table.biom', 'taxonomy.tsv', 'table-with-taxo.biom')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['biom', 'add-metadata', '--sc-separated', 'taxonomy',
                    '-i', 'table.biom', '--observation-metadata-fp', 'taxonomy.tsv',
                    '-o', 'table-with-taxo.biom']


def test_convert_to_tsv_command(mock_run):
    biom_utils.convert_to_tsv('table.biom', 'taxonomy.tsv', 'table.biom.tsv')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['biom', 'convert', '--to-tsv', '--header-key', 'taxonomy',
                    '-i', 'table.biom', '--observation-metadata-fp', 'taxonomy.tsv',
                    '-o', 'table.biom.tsv']
