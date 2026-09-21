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


class TestReadMetadataRows:
    def test_comment_rows_are_not_samples(self, tmp_path):
        """QIIME2 ignores "#"-prefixed rows (anywhere, including before the
        header) and empty rows; they used to be read as samples here."""
        path = tmp_path / 'metadata.tsv'
        path.write_text(
            '# exported from the LIMS\n'
            'sample-id\tsite\n'
            '#q2:types\tcategorical\n'
            '# second batch below\n'
            '\t\n'
            'sampleA\tsiteA\n'
        )
        assert metadata_utils.read_metadata_table(path) == {'sampleA': {'site': 'siteA'}}
        assert metadata_utils.parse_metadata_columns(path) == {'site': 'categorical'}

    def test_legacy_hash_id_header_is_the_header_not_a_comment(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        path.write_text('#SampleID\tsite\nsampleA\tsiteA\n')
        assert metadata_utils.read_metadata_table(path) == {'sampleA': {'site': 'siteA'}}

    def test_cells_are_stripped_and_quotes_honored(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        path.write_text('sample-id\tsite\nsampleA \t"site A"\n')
        assert metadata_utils.read_metadata_table(path) == {'sampleA': {'site': 'site A'}}

    def test_raises_clearly_without_a_header(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        path.write_text('# only a comment\n')
        with pytest.raises(ValueError, match='no header'):
            metadata_utils.read_metadata_table(path)


class TestClassSizes:
    def test_missing_values_are_not_a_group(self, tmp_path):
        """Regression test: '' used to count as a class, so one real group
        plus two blanks read as two groups of >=2 -- "eligible" -- and
        beta-group-significance then failed the whole run on that column."""
        path = tmp_path / 'metadata.tsv'
        path.write_text(
            'sample-id\tsite\n#q2:types\tcategorical\n'
            'sampleA\tsiteA\nsampleB\tsiteA\nsampleC\t\nsampleD\t\n'
        )
        ids = ['sampleA', 'sampleB', 'sampleC', 'sampleD']
        assert metadata_utils.class_sizes(path, 'site', ids) == {'siteA': 2}
        assert metadata_utils.eligible_categorical_columns(path, ids) == []
        assert metadata_utils.has_alpha_group_significance_column(path, ids) is False

    def test_counts_values_for_given_samples(self, metadata_with_types):
        counts = metadata_utils.class_sizes(metadata_with_types, 'site',
                                             ['sampleA', 'sampleB', 'sampleC'])
        assert counts == {'siteA': 2, 'siteB': 1}

    def test_ignores_sample_ids_not_in_metadata(self, metadata_with_types):
        counts = metadata_utils.class_sizes(metadata_with_types, 'site', ['sampleA', 'unknown'])
        assert counts == {'siteA': 1}


class TestFinalSampleIds:
    def test_filters_out_zero_read_samples(self):
        sample_frequencies = {'sampleA': 19.0, 'sampleB': 0.0, 'sampleC': 7.0}
        assert metadata_utils.final_sample_ids(sample_frequencies, ['sampleA', 'sampleB', 'sampleC']) \
            == ['sampleA', 'sampleC']

    def test_falls_back_to_every_sample_when_frequencies_missing(self):
        assert metadata_utils.final_sample_ids({}, ['sampleA', 'sampleB']) == ['sampleA', 'sampleB']

    def test_all_zero_read_samples_returns_empty_not_the_fallback(self):
        """Regression test: a naive `filtered or fallback_sample_ids`
        implementation can't distinguish "frequencies were never loaded"
        from "every sample loaded at zero reads" -- both produce an empty
        filtered list, but only the first should fall back to every sample.
        A run where literally every sample ended up zero-read must report
        zero eligible samples, not silently un-exclude them all."""
        sample_frequencies = {'sampleA': 0.0, 'sampleB': 0.0}
        assert metadata_utils.final_sample_ids(sample_frequencies, ['sampleA', 'sampleB']) == []


class TestEligibleCategoricalColumns:
    def test_excludes_numeric_columns(self, metadata_with_types):
        """'ph' in this fixture is both numeric-typed AND has all-unique
        values (5.8/5.9/6.4), so this alone can't tell whether the type
        check or the group-count check is what excludes it -- see the
        dedicated test below for a fixture that isolates the type check."""
        eligible = metadata_utils.eligible_categorical_columns(
            metadata_with_types, ['sampleA', 'sampleB', 'sampleC'])
        assert 'ph' not in eligible

    def test_excludes_numeric_columns_even_with_repeated_values(self, tmp_path):
        """A numeric column with a repeated value (here, replicate numbers
        1,1,2,2 -- two groups of 2, which would otherwise satisfy the
        eligibility rule) isolates the type check: without it, this column
        would read as eligible on group counts alone."""
        path = tmp_path / 'metadata.tsv'
        path.write_text(
            'sample-id\treplicate-number\n#q2:types\tnumeric\n'
            'sampleA\t1\nsampleB\t1\nsampleC\t2\nsampleD\t2\n'
        )
        eligible = metadata_utils.eligible_categorical_columns(
            path, ['sampleA', 'sampleB', 'sampleC', 'sampleD'])
        assert eligible == []

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

    def test_only_considers_given_sample_ids(self, tmp_path):
        """A near-empty/excluded sample shouldn't count toward group
        eligibility -- verified by a fixture where excluding one sample
        actually flips the result: with all 4 samples, 'site' has two
        groups of 2 (eligible); excluding sampleD drops siteB to a
        singleton, and only one group is then >= min_per_group (not
        eligible). (A fixture where exclusion can't change the outcome
        either way doesn't prove sample_ids is actually being applied.)"""
        path = tmp_path / 'metadata.tsv'
        path.write_text(
            'sample-id\tsite\n#q2:types\tcategorical\n'
            'sampleA\tsiteA\nsampleB\tsiteA\nsampleC\tsiteB\nsampleD\tsiteB\n'
        )
        eligible_all = metadata_utils.eligible_categorical_columns(
            path, ['sampleA', 'sampleB', 'sampleC', 'sampleD'])
        assert 'site' in eligible_all

        eligible_excluding_sampleD = metadata_utils.eligible_categorical_columns(
            path, ['sampleA', 'sampleB', 'sampleC'])
        assert 'site' not in eligible_excluding_sampleD


class TestReadsTheMetadataFileOnce:
    """eligible_categorical_columns() and has_alpha_group_significance_column()
    used to re-read and re-parse the whole metadata TSV once per categorical
    column (class_sizes() -> read_metadata_table()), on top of
    parse_metadata_columns()'s own read -- 1 + N file reads per call for N
    categorical columns (4 for this 3-categorical-column fixture)."""

    @pytest.fixture
    def read_calls(self, monkeypatch):
        calls = []
        real_read = metadata_utils.read_metadata_rows

        def counting_read(path):
            calls.append(path)
            return real_read(path)

        monkeypatch.setattr(metadata_utils, 'read_metadata_rows', counting_read)
        return calls

    def test_eligible_categorical_columns(self, metadata_with_types, read_calls):
        metadata_utils.eligible_categorical_columns(metadata_with_types, ['sampleA', 'sampleB', 'sampleC'])
        assert len(read_calls) == 1

    def test_has_alpha_group_significance_column(self, metadata_with_types, read_calls):
        metadata_utils.has_alpha_group_significance_column(metadata_with_types, ['sampleA', 'sampleB', 'sampleC'])
        assert len(read_calls) == 1


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

    def test_excluded_sample_does_not_count_toward_the_repeat(self, tmp_path):
        """A sample outside sample_ids must not count toward the
        repeated-value check -- verified by a fixture where excluding the
        one sample that creates the repeat actually flips the result from
        True to False (sampleA/sampleB both 'A' when both count; excluding
        sampleB leaves every remaining value unique)."""
        path = tmp_path / 'metadata.tsv'
        path.write_text(
            'sample-id\tcondition\n#q2:types\tcategorical\n'
            'sampleA\tA\nsampleB\tA\nsampleC\tB\n'
        )
        assert metadata_utils.has_alpha_group_significance_column(
            path, ['sampleA', 'sampleB', 'sampleC']) is True

        assert metadata_utils.has_alpha_group_significance_column(
            path, ['sampleA', 'sampleC']) is False
