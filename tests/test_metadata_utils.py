import pytest

from qiime2_its import metadata_utils


@pytest.fixture
def metadata_with_types(tmp_path):
    path = tmp_path / 'metadata.tsv'
    path.write_text(
        'sample-id\tsite\treplicate\tph\tcollection-date\n'
        '#q2:types\tcategorical\tcategorical\tnumeric\tcategorical\n'
        'sampleA\tsiteA\t1\t5.8\t2026-06-02\n'
        'sampleB\tsiteA\t2\t5.9\t2026-06-02\n'
        'sampleC\tsiteB\t1\t6.4\t2026-06-03\n'
    )
    return path


@pytest.fixture
def metadata_without_types(tmp_path):
    path = tmp_path / 'metadata.tsv'
    path.write_text(
        'sample-id\tcondition\televation-m\n'
        'sampleA\tA\t1340\n'
        'sampleB\tB\t980\n'
    )
    return path


class TestParseMetadataColumns:
    def test_honors_explicit_types_row(self, metadata_with_types):
        columns = metadata_utils.parse_metadata_columns(metadata_with_types)
        assert columns == {
            'site': 'categorical',
            'replicate': 'categorical',
            'ph': 'numeric',
            'collection-date': 'categorical',
        }

    def test_infers_types_when_no_types_row(self, metadata_without_types):
        columns = metadata_utils.parse_metadata_columns(metadata_without_types)
        assert columns == {'condition': 'categorical', 'elevation-m': 'numeric'}

    def test_missing_values_do_not_block_numeric_inference(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        path.write_text('sample-id\televation-m\nsampleA\t1340\nsampleB\t\n')
        columns = metadata_utils.parse_metadata_columns(path)
        assert columns == {'elevation-m': 'numeric'}


class TestReadMetadataTable:
    def test_reads_rows_keyed_by_sample_id(self, metadata_with_types):
        table = metadata_utils.read_metadata_table(metadata_with_types)
        assert table['sampleA'] == {
            'site': 'siteA', 'replicate': '1', 'ph': '5.8', 'collection-date': '2026-06-02'}
        assert table['sampleC']['site'] == 'siteB'


class TestClassSizes:
    def test_counts_values_for_given_samples(self, metadata_with_types):
        counts = metadata_utils.class_sizes(metadata_with_types, 'site',
                                             ['sampleA', 'sampleB', 'sampleC'])
        assert counts == {'siteA': 2, 'siteB': 1}

    def test_ignores_sample_ids_not_in_metadata(self, metadata_with_types):
        counts = metadata_utils.class_sizes(metadata_with_types, 'site', ['sampleA', 'unknown'])
        assert counts == {'siteA': 1}


class TestEligibleCategoricalColumns:
    def test_excludes_numeric_columns(self, metadata_with_types):
        eligible = metadata_utils.eligible_categorical_columns(
            metadata_with_types, ['sampleA', 'sampleB', 'sampleC'])
        assert 'ph' not in eligible

    def test_excludes_columns_with_a_singleton_group(self, metadata_with_types):
        """'replicate' has values 1,2,1 across 3 samples -- '2' is a group of
        one, so only one group would have >=2 members: not eligible."""
        eligible = metadata_utils.eligible_categorical_columns(
            metadata_with_types, ['sampleA', 'sampleB', 'sampleC'])
        assert 'replicate' not in eligible

    def test_includes_column_with_two_balanced_groups(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        path.write_text(
            'sample-id\tsite\n'
            '#q2:types\tcategorical\n'
            'sampleA\tsiteA\nsampleB\tsiteA\nsampleC\tsiteB\nsampleD\tsiteB\n'
        )
        eligible = metadata_utils.eligible_categorical_columns(
            path, ['sampleA', 'sampleB', 'sampleC', 'sampleD'])
        assert eligible == ['site']

    def test_respects_min_per_group_override(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        path.write_text(
            'sample-id\tsite\n'
            '#q2:types\tcategorical\n'
            'sampleA\tsiteA\nsampleB\tsiteA\nsampleC\tsiteB\nsampleD\tsiteB\n'
        )
        eligible = metadata_utils.eligible_categorical_columns(
            path, ['sampleA', 'sampleB', 'sampleC', 'sampleD'], min_per_group=3)
        assert eligible == []

    def test_only_considers_given_sample_ids(self, metadata_with_types):
        """A near-empty/excluded sample shouldn't count toward group eligibility."""
        eligible = metadata_utils.eligible_categorical_columns(
            metadata_with_types, ['sampleA', 'sampleC'])  # sampleB excluded
        assert 'site' not in eligible  # siteA now has only 1 member (sampleA)


class TestHasAlphaGroupSignificanceColumn:
    def test_false_when_every_sample_has_a_unique_value(self, tmp_path):
        """Regression test: qiime diversity alpha-group-significance fails
        outright (for the whole metadata file, not a graceful per-column
        skip) if no categorical column has a repeated value -- exactly the
        case for a small 2-sample metadata file with one ID-like column."""
        path = tmp_path / 'metadata.tsv'
        path.write_text('sample-id\tcondition\n#q2:types\tcategorical\nsampleA\tA\nsampleB\tB\n')
        assert metadata_utils.has_alpha_group_significance_column(path, ['sampleA', 'sampleB']) is False

    def test_false_when_column_is_constant(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        path.write_text('sample-id\tcondition\n#q2:types\tcategorical\nsampleA\tA\nsampleB\tA\n')
        assert metadata_utils.has_alpha_group_significance_column(path, ['sampleA', 'sampleB']) is False

    def test_true_when_a_value_repeats_but_not_all(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        path.write_text(
            'sample-id\tsite\n#q2:types\tcategorical\n'
            'sampleA\tsiteA\nsampleB\tsiteA\nsampleC\tsiteB\n'
        )
        assert metadata_utils.has_alpha_group_significance_column(
            path, ['sampleA', 'sampleB', 'sampleC']) is True
