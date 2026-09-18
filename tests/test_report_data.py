import json
import zipfile

import pytest

from qiime2_its import report_data


class TestParseSampleFrequencies:
    def test_parses_frequencies_including_zero(self, tmp_path):
        """Regression test: a near-empty input sample can survive DADA2 as a
        zero-read row (retain-all-samples defaults to True). Downstream
        eligibility checks need to see that zero, not skip the row."""
        path = tmp_path / 'sample-frequencies.tsv'
        path.write_text(
            'Sample ID\tFrequency\tNo. of Associated Features\n'
            '#q2:types\tcategorical\tcategorical\n'
            'siteA-rep1\t19.0\t3\n'
            'siteC-rep2\t0.0\t0\n'
        )
        frequencies = report_data.parse_sample_frequencies(path)
        assert frequencies == {'siteA-rep1': 19.0, 'siteC-rep2': 0.0}

    def test_parses_thousands_separated_frequencies(self, tmp_path):
        """Regression test: QIIME2 writes this file with a comma thousands
        separator once a sample's frequency reaches four digits (e.g.
        "110,406.0") -- found running against a real 53-sample production
        dataset whose samples had ~100k+ reads each; every fixture used
        before that had frequencies too small to trigger it."""
        path = tmp_path / 'sample-frequencies.tsv'
        path.write_text(
            'Sample ID\tFrequency\tNo. of Associated Features\n'
            '#q2:types\tcategorical\tcategorical\n'
            'K.BeL.1.1\t110,406.0\t425\n'
            'P.Leg.1.1\t389,161.0\t494\n'
        )
        frequencies = report_data.parse_sample_frequencies(path)
        assert frequencies == {'K.BeL.1.1': 110406.0, 'P.Leg.1.1': 389161.0}


class TestParseDada2Stats:
    def test_parses_and_coerces_numeric_columns(self, tmp_path):
        path = tmp_path / 'stats.tsv'
        path.write_text(
            'sample-id\tinput\tnon-chimeric\tpercentage of input non-chimeric\n'
            '#q2:types\tnumeric\tnumeric\tnumeric\n'
            'sampleA\t100\t80\t80.0\n'
            'sampleB\t50\t0\t0\n'
        )
        df = report_data.parse_dada2_stats(path)
        assert list(df.index) == ['sampleA', 'sampleB']
        assert df.loc['sampleA', 'input'] == 100
        assert df.loc['sampleB', 'non-chimeric'] == 0


class TestParseOrdination:
    def test_parses_site_coordinates_and_proportion_explained(self, tmp_path):
        path = tmp_path / 'ordination.txt'
        path.write_text(
            'Eigvals\t2\n'
            '0.5\t0.03\n'
            '\n'
            'Proportion explained\t2\n'
            '0.93\t0.06\n'
            '\n'
            'Species\t0\t0\n'
            '\n'
            'Site\t2\t2\n'
            'sampleA\t0.45\t-0.08\n'
            'sampleB\t-0.54\t-0.06\n'
            '\n'
            'Biplot\t0\t0\n'
            '\n'
            'Site constraints\t0\t0\n'
        )
        sample_coords, proportion_explained = report_data.parse_ordination(path)
        assert sample_coords == {'sampleA': (0.45, -0.08), 'sampleB': (-0.54, -0.06)}
        assert proportion_explained == (0.93, 0.06)


def _write_alpha_group_significance_qzv(path, column_stats):
    """column_stats: {column_name: (h, p, {group_label: [values]})}."""
    with zipfile.ZipFile(path, 'w') as zf:
        for column, (h, p, groups) in column_stats.items():
            index = list(groups.keys())
            data = list(groups.values())
            group_data = json.dumps({'name': None, 'index': index, 'data': data})
            content = (
                f"load_data('{column}',{group_data},"
                f'{{"initial": 4, "filtered": 4}},'
                f'{{"H": {h}, "p": {p}}},'
                f"'<table></table>','kruskal-wallis-pairwise-{column}.csv', 'shannon_entropy');"
            )
            zf.writestr(f'uuid1234/data/column-{column}.jsonp', content)


class TestParseAlphaGroupSignificance:
    def test_extracts_h_and_p_per_column(self, tmp_path):
        qzv = tmp_path / 'alpha.qzv'
        _write_alpha_group_significance_qzv(qzv, {
            'site': (0.9166666666666659, 0.6323366621862501,
                     {'siteA (n=2)': [0.97, 1.52], 'siteB (n=1)': [0.97]}),
            'host-plant': (1.5, 0.22, {'Pinus (n=1)': [1.1]}),
        })
        result = report_data.parse_alpha_group_significance(qzv)
        assert result['site']['h_statistic'] == 0.9166666666666659
        assert result['site']['p_value'] == 0.6323366621862501
        assert result['host-plant']['h_statistic'] == 1.5

    def test_extracts_raw_group_values_with_sample_counts_stripped(self, tmp_path):
        qzv = tmp_path / 'alpha.qzv'
        _write_alpha_group_significance_qzv(qzv, {
            'site': (0.9, 0.6, {'siteA (n=2)': [0.97, 1.52], 'siteB (n=1)': [0.97]}),
        })
        result = report_data.parse_alpha_group_significance(qzv)
        assert result['site']['groups'] == {'siteA': [0.97, 1.52], 'siteB': [0.97]}


def _write_beta_group_significance_qzv(path, overview_rows):
    html = '<html><body><table>' + ''.join(
        f'<th>{k}</th><td>{v}</td>' for k, v in overview_rows.items()
    ) + '</table></body></html>'
    with zipfile.ZipFile(path, 'w') as zf:
        zf.writestr('uuid5678/data/index.html', html)


class TestParseBetaGroupSignificance:
    def test_extracts_overview_table(self, tmp_path):
        qzv = tmp_path / 'beta.qzv'
        _write_beta_group_significance_qzv(qzv, {
            'method name': 'PERMANOVA',
            'test statistic name': 'pseudo-F',
            'sample size': '6',
            'number of groups': '3',
            'test statistic': '2.345',
            'p-value': '0.012',
            'number of permutations': '999',
        })
        result = report_data.parse_beta_group_significance(qzv)
        assert result == {
            'method_name': 'PERMANOVA',
            'test_statistic_name': 'pseudo-F',
            'sample_size': 6,
            'number_of_groups': 3,
            'test_statistic': 2.345,
            'p_value': 0.012,
        }


def _write_rarefaction_qzv(path, metric, rows):
    """rows: {sample_id: {col_name: value_str}}."""
    fieldnames = sorted({col for row in rows.values() for col in row})
    lines = ['sample-id,' + ','.join(fieldnames)]
    for sample_id, row in rows.items():
        lines.append(sample_id + ',' + ','.join(row.get(c, '') for c in fieldnames))
    with zipfile.ZipFile(path, 'w') as zf:
        zf.writestr(f'uuidabcd/data/{metric}.csv', '\n'.join(lines))


def _write_classifier_accuracy_qzv(path, tsv_content):
    with zipfile.ZipFile(path, 'w') as zf:
        zf.writestr('uuidef01/data/predictive_accuracy.tsv', tsv_content)


class TestParseClassifierAccuracy:
    def test_extracts_summary_rows_ignoring_confusion_matrix(self, tmp_path):
        """Regression test: the same numbers also appear in an HTML
        confusion-matrix table with a variable number of columns per class;
        a fixed <th>/<td> pair extraction (as used for beta-group-
        significance) would grab an empty leading cell instead. The
        plain-text export's summary rows always end with the real value."""
        qzv = tmp_path / 'accuracy.qzv'
        _write_classifier_accuracy_qzv(qzv,
            '\t1\t2\tOverall Accuracy\n'
            '1\t1.0\t0.0\t\n'
            '2\t1.0\t0.0\t\n'
            'Overall Accuracy\t\t\t0.3333333333333333\n'
            'Baseline Accuracy\t\t\t0.6666666666666666\n'
            'Accuracy Ratio\t\t\t0.5\n')
        result = report_data.parse_classifier_accuracy(qzv)
        assert result == {
            'overall_accuracy': 0.3333333333333333,
            'baseline_accuracy': 0.6666666666666666,
            'accuracy_ratio': 0.5,
        }


class TestParseRarefactionCurve:
    def test_averages_iterations_per_depth_sorted_by_depth(self, tmp_path):
        qzv = tmp_path / 'raref.qzv'
        _write_rarefaction_qzv(qzv, 'observed_features', {
            'sampleA': {
                'depth-1_iter-1': '1.0', 'depth-1_iter-2': '1.0',
                'depth-5_iter-1': '3.0', 'depth-5_iter-2': '5.0',
            },
        })
        curves = report_data.parse_rarefaction_curve(qzv, 'observed_features')
        assert curves['sampleA'] == [(1, 1.0), (5, 4.0)]

    def test_skips_depths_beyond_a_samples_reads(self, tmp_path):
        """A shallow sample has blank cells at depths it can't reach -- those
        should be skipped, not treated as zero."""
        qzv = tmp_path / 'raref.qzv'
        _write_rarefaction_qzv(qzv, 'shannon', {
            'sampleA': {'depth-1_iter-1': '1.0'},
            'sampleB': {},  # blank at every depth column
        })
        curves = report_data.parse_rarefaction_curve(qzv, 'shannon')
        assert curves['sampleB'] == []


class TestBuildGenusAbundanceTable:
    def test_collapses_to_genus_and_normalizes_per_sample(self, tmp_path):
        path = tmp_path / 'table-with-taxonomy.biom.tsv'
        path.write_text(
            '# Constructed from biom file\n'
            '#OTU ID\tsampleA\tsampleB\ttaxonomy\n'
            'f1\t80.0\t0.0\tk__Fungi; p__Ascomycota; c__Sordariomycetes; g__Fusarium\n'
            'f2\t20.0\t50.0\tk__Fungi; p__Ascomycota; c__Sordariomycetes; g__Fusarium\n'
            'f3\t0.0\t50.0\tk__Fungi; p__Basidiomycota; g__Trichosporon\n'
        )
        table = report_data.build_genus_abundance_table(path)
        assert table.loc['Fusarium', 'sampleA'] == pytest.approx(1.0)
        assert table.loc['Fusarium', 'sampleB'] == pytest.approx(0.5)
        assert table.loc['Trichosporon', 'sampleB'] == pytest.approx(0.5)

    def test_unidentified_and_missing_genus_become_unclassified(self, tmp_path):
        path = tmp_path / 'table-with-taxonomy.biom.tsv'
        path.write_text(
            '# Constructed from biom file\n'
            '#OTU ID\tsampleA\ttaxonomy\n'
            'f1\t10.0\tk__Fungi; p__Ascomycota; g__unidentified\n'
            'f2\t10.0\tk__Fungi; p__Ascomycota\n'
        )
        table = report_data.build_genus_abundance_table(path)
        assert table.loc['Unclassified', 'sampleA'] == pytest.approx(1.0)

    def test_collapses_low_abundance_genera_into_other(self, tmp_path):
        rows = ''.join(
            f'f{i}\t1.0\tk__Fungi; g__Genus{i}\n' for i in range(15)
        )
        path = tmp_path / 'table-with-taxonomy.biom.tsv'
        path.write_text(
            '# Constructed from biom file\n#OTU ID\tsampleA\ttaxonomy\n' + rows
        )
        table = report_data.build_genus_abundance_table(path, top_n=10)
        assert len(table) == 11  # top 10 + Other
        assert 'Other' in table.index
        assert table['sampleA'].sum() == pytest.approx(1.0)


class TestParseRunMetadata:
    def test_parses_written_json(self, tmp_path):
        path = tmp_path / 'run_metadata.json'
        path.write_text('{"pipeline": {"qiime2_its_version": "0.2.0"}}')
        assert report_data.parse_run_metadata(path) == {'pipeline': {'qiime2_its_version': '0.2.0'}}

    def test_returns_none_when_missing(self, tmp_path):
        assert report_data.parse_run_metadata(tmp_path / 'does-not-exist.json') is None


class TestParseFastaSequenceLengths:
    def test_parses_single_line_sequences(self, tmp_path):
        path = tmp_path / 'dna-sequences.fasta'
        path.write_text('>seq1\nACGTACGT\n>seq2\nACGT\n')
        assert report_data.parse_fasta_sequence_lengths(path) == [8, 4]

    def test_sums_a_sequence_wrapped_across_multiple_lines(self, tmp_path):
        """A record's sequence isn't guaranteed to be a single line -- some
        FASTA writers wrap at a fixed width -- so lines must be accumulated
        between headers rather than assumed to be one sequence per line."""
        path = tmp_path / 'dna-sequences.fasta'
        path.write_text('>seq1\nACGT\nACGT\n>seq2\nAC\n')
        assert report_data.parse_fasta_sequence_lengths(path) == [8, 2]

    def test_empty_file_returns_empty_list(self, tmp_path):
        path = tmp_path / 'dna-sequences.fasta'
        path.write_text('')
        assert report_data.parse_fasta_sequence_lengths(path) == []


class TestParseTaxonomyConfidence:
    def test_parses_confidence_column(self, tmp_path):
        path = tmp_path / 'taxonomy.tsv'
        path.write_text(
            '#OTUID\ttaxonomy\tconfidence\n'
            'f1\tk__Fungi;p__Ascomycota\t0.999\n'
            'f2\tk__Fungi\t1.0\n'
        )
        assert report_data.parse_taxonomy_confidence(path) == [0.999, 1.0]
