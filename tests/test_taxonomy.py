import gzip

import pytest

from qiime2_its import taxonomy


class TestParseIdTable:
    def test_parses_two_column_tsv_into_taxid_keyed_dict(self, tmp_path):
        table = tmp_path / 'id_table.tsv'
        table.write_text('ACC1\t1001\nACC2\t1002\n')
        assert taxonomy.parse_id_table(table) == {'1001': 'ACC1', '1002': 'ACC2'}

    def test_skips_blank_lines(self, tmp_path):
        table = tmp_path / 'id_table.tsv'
        table.write_text('ACC1\t1001\n\nACC2\t1002\n')
        assert taxonomy.parse_id_table(table) == {'1001': 'ACC1', '1002': 'ACC2'}

    def test_rejects_wrong_column_count(self, tmp_path):
        table = tmp_path / 'id_table.tsv'
        table.write_text('ACC1\t1001\textra\n')
        with pytest.raises(ValueError):
            taxonomy.parse_id_table(table)

    def test_handles_gzipped_input(self, tmp_path):
        table = tmp_path / 'id_table.tsv.gz'
        with gzip.open(table, 'wt') as f:
            f.write('ACC1\t1001\n')
        assert taxonomy.parse_id_table(table) == {'1001': 'ACC1'}


class TestExtractAccessionsFromFasta:
    def test_extracts_accession_before_first_whitespace(self, tmp_path):
        fasta = tmp_path / 'seqs.fasta'
        fasta.write_text('>ACC1.1 some description\nACGT\n>ACC2.2\nTTTT\n')
        acc_file = tmp_path / 'acc.list'

        acc_dict = taxonomy.extract_accessions_from_fasta(fasta, acc_file)

        assert acc_dict == {'ACC1.1': '', 'ACC2.2': ''}
        assert acc_file.read_text() == 'ACC1.1\nACC2.2\n'


class TestParseAccession2Taxid:
    def test_matches_on_unversioned_accession(self, tmp_path):
        acc2taxid = tmp_path / 'nucl.accession2taxid'
        acc2taxid.write_text('accession\taccession.version\ttaxid\tgi\n'
                              'ACC1\tACC1.1\t1001\t0\n'
                              'ACC2\tACC2.3\t1002\t0\n'
                              'UNRELATED\tUNRELATED.1\t9999\t0\n')

        result = taxonomy.parse_accession2taxid(acc2taxid, ['ACC1.1', 'ACC2.9'])

        assert result == {'ACC1': '1001', 'ACC2': '1002'}


class TestAccessionsToTaxids:
    def test_falls_back_to_unclassified_for_missing(self, tmp_path):
        out_file = tmp_path / 'taxid.list'
        missing = taxonomy.accessions_to_taxids(
            acc_dict={'ACC1.1': '', 'ACC2.1': ''},
            acc2taxid_dict={'ACC1': '1001'},
            output_file=out_file,
        )
        assert missing == ['ACC2.1']
        assert out_file.read_text() == 'ACC1.1\t1001\nACC2.1\t12908\n'


class TestTaxdumpParsing:
    @pytest.fixture
    def taxdump_files(self, tmp_path):
        # Mirrors the real NCBI taxdump shape for Fungi: a 'clade' node
        # (Opisthokonta-like, not a Linnean rank) sits between kingdom and
        # phylum and must NOT be treated as a class.
        nodes = tmp_path / 'nodes.dmp'
        nodes.write_text(
            '1\t|\t1\t|\tno rank\t|\n'
            '4751\t|\t1\t|\tkingdom\t|\n'
            '33154\t|\t4751\t|\tclade\t|\n'
            '4890\t|\t33154\t|\tphylum\t|\n'
            '147545\t|\t4890\t|\tclass\t|\n'
        )
        names = tmp_path / 'names.dmp'
        names.write_text(
            '4751\t|\tFungi\t|\t\t|\tscientific name\t|\n'
            '33154\t|\tOpisthokonta\t|\t\t|\tscientific name\t|\n'
            '4890\t|\tAscomycota\t|\t\t|\tscientific name\t|\n'
            '147545\t|\tSaccharomycetes\t|\t\t|\tscientific name\t|\n'
        )
        merged = tmp_path / 'merged.dmp'
        merged.write_text('')
        return nodes, names, merged

    def test_parse_nodes_dmp(self, taxdump_files):
        nodes_file, _, _ = taxdump_files
        node_dict = taxonomy.parse_nodes_dmp(nodes_file)
        assert node_dict['4890'] == ('33154', 'phylum')

    def test_parse_names_dmp_only_keeps_scientific_name(self, taxdump_files):
        _, names_file, _ = taxdump_files
        names_dict = taxonomy.parse_names_dmp(names_file)
        assert names_dict == {'4751': 'Fungi', '33154': 'Opisthokonta', '4890': 'Ascomycota',
                               '147545': 'Saccharomycetes'}

    def test_apply_merged_taxids_remaps_old_to_new(self, tmp_path):
        merged = tmp_path / 'merged.dmp'
        merged.write_text('111\t|\t222\t|\n')
        result = taxonomy.apply_merged_taxids({'111': 'ACC1'}, merged)
        assert result == {'222': 'ACC1'}
        assert '111' not in result

    def test_lineage_string_skips_the_intervening_clade(self, taxdump_files):
        """Regression test: a real NCBI lineage has a 'clade' node
        (Opisthokonta) between kingdom and phylum. It must be skipped, not
        mistaken for class -- 'clade' is not a Linnean rank."""
        nodes_file, names_file, _ = taxdump_files
        node_dict = taxonomy.parse_nodes_dmp(nodes_file)
        names_dict = taxonomy.parse_names_dmp(names_file)
        lineage = taxonomy.lineage_string('4890', node_dict, names_dict)
        assert lineage == 'k__Fungi;p__Ascomycota;c__unidentified;o__unidentified;' \
                           'f__unidentified;g__unidentified;s__unidentified'
        assert 'Opisthokonta' not in lineage

    def test_lineage_string_maps_class_rank_correctly(self, taxdump_files):
        """Regression test: the original rank map used NCBI rank 'clade' for
        the class slot and 'superkingdom' for kingdom -- neither matches real
        NCBI taxonomy, which uses 'class' and 'kingdom'."""
        nodes_file, names_file, _ = taxdump_files
        node_dict = taxonomy.parse_nodes_dmp(nodes_file)
        names_dict = taxonomy.parse_names_dmp(names_file)
        lineage = taxonomy.lineage_string('147545', node_dict, names_dict)
        assert lineage == 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__unidentified;' \
                           'f__unidentified;g__unidentified;s__unidentified'

    def test_write_taxonomy_file(self, tmp_path, taxdump_files):
        nodes_file, names_file, merged_file = taxdump_files
        taxonomy_file = tmp_path / 'taxonomy.txt'
        taxonomy.write_taxonomy_file({'4890': 'ACC1'}, taxonomy_file, nodes_file, names_file, merged_file)
        line = taxonomy_file.read_text().strip()
        assert line.startswith('ACC1\tk__Fungi;p__Ascomycota;')


class TestMaxLineageDepth:
    def test_returns_the_deepest_lineage_present(self, tmp_path):
        """Regression test: `qiime taxa collapse --p-level N` fails outright
        if N exceeds the deepest lineage any feature actually has -- common
        when a classifier can't resolve reads unrelated to its training set
        past family/order, well short of the genus level this pipeline
        otherwise requests by default."""
        path = tmp_path / 'taxonomy.tsv'
        path.write_text(
            '#OTUID\ttaxonomy\tconfidence\n'
            'f1\tk__Fungi;p__Ascomycota\t1.0\n'
            'f2\tk__Fungi;p__Ascomycota;c__Sordariomycetes;o__Sordariales\t0.8\n'
            'f3\tk__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae\t0.79\n'
        )
        assert taxonomy.max_lineage_depth(path) == 5

    def test_returns_zero_for_header_only_file(self, tmp_path):
        path = tmp_path / 'taxonomy.tsv'
        path.write_text('#OTUID\ttaxonomy\tconfidence\n')
        assert taxonomy.max_lineage_depth(path) == 0
