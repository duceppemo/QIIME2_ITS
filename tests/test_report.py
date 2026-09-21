"""Integration-style test for report.build_report(): builds a full PDF from a
synthetic (but realistically-shaped) output_folder, so the whole
report_data -> matplotlib -> fpdf2 pipeline is exercised in CI without
needing QIIME2 installed. Regression coverage for exactly the kind of bug
real-data testing caught: fpdf2's multi_cell() not resetting the cursor to
the left margin between calls the way cell(new_x='LMARGIN') does, which blew
up on the second line of the title page.
"""
import json
import zipfile

from qiime2_its import report


def _page_titles(mocker):
    """Spy on every _ReportPDF page-adding method (add_title_page,
    add_table_page, add_keyvalue_page, add_figure_page -- title is always
    their first argument after self) and return a list that accumulates
    each one's title, in call order, as build_report() runs. Lets
    TestBuildReport assert on which sections actually made it into the
    report instead of only on the resulting PDF's byte size, which stays
    roughly stable (or even grows) whether or not a given section's code
    path actually ran -- e.g. dropping the classifier-results table
    entirely, or never emitting the alpha boxplot pages, changes the PDF by
    only a few hundred bytes against a report that's otherwise many
    matplotlib figures, well within the noise `> 1000 bytes`-style
    assertions can't see."""
    titles = []
    for method_name in ('add_title_page', 'add_table_page', 'add_keyvalue_page', 'add_figure_page'):
        original = getattr(report._ReportPDF, method_name)

        def wrapper(self, title, *args, _original=original, **kwargs):
            titles.append(title)
            return _original(self, title, *args, **kwargs)

        mocker.patch.object(report._ReportPDF, method_name, wrapper)
    return titles


def _write_alpha_qzv(path, column_stats):
    with zipfile.ZipFile(path, 'w') as zf:
        for column, (h, p, groups) in column_stats.items():
            group_data = json.dumps({'name': None, 'index': list(groups), 'data': list(groups.values())})
            content = (f"load_data('{column}',{group_data},{{}},"
                       f'{{"H": {h}, "p": {p}}},\'<table></table>\',\'x.csv\', \'metric\');')
            zf.writestr(f'uuid1/data/column-{column}.jsonp', content)


def _write_beta_qzv(path, overview_rows):
    # Real q2templates output (pandas DataFrame.to_html()) puts each <td> on
    # its own indented line after its <th>, not directly adjacent -- see the
    # matching comment in test_report_data.py's _write_beta_group_significance_qzv.
    html = '<html><body><table>' + ''.join(
        f'<tr>\n      <th>{k}</th>\n      <td>{v}</td>\n    </tr>' for k, v in overview_rows.items()
    ) + '</table></body></html>'
    with zipfile.ZipFile(path, 'w') as zf:
        zf.writestr('uuid2/data/index.html', html)


def _write_rarefaction_qzv(path, metrics_rows):
    """metrics_rows: {metric: {sample_id: {col: value}}}."""
    with zipfile.ZipFile(path, 'w') as zf:
        for metric, rows in metrics_rows.items():
            fieldnames = sorted({c for row in rows.values() for c in row})
            lines = ['sample-id,' + ','.join(fieldnames)]
            for sample_id, row in rows.items():
                lines.append(sample_id + ',' + ','.join(row.get(c, '') for c in fieldnames))
            zf.writestr(f'uuid3/data/{metric}.csv', '\n'.join(lines))


def _write_classifier_qzv(path, accuracy_tsv):
    with zipfile.ZipFile(path, 'w') as zf:
        zf.writestr('uuid4/data/predictive_accuracy.tsv', accuracy_tsv)


def _build_synthetic_output_folder(tmp_path):
    out = tmp_path / 'output'
    (out / 'dada2_stats').mkdir(parents=True)
    (out / 'dada2_stats' / 'stats.tsv').write_text(
        'sample-id\tinput\tnon-chimeric\tpercentage of input non-chimeric\n'
        '#q2:types\tnumeric\tnumeric\tnumeric\n'
        'sampleA\t100\t80\t80.0\n'
        'sampleB\t90\t60\t66.7\n'
        'sampleC\t95\t70\t73.7\n'
        'sampleD\t85\t50\t58.8\n'
    )

    (out / 'sample_frequencies').mkdir(parents=True)
    (out / 'sample_frequencies' / 'metadata.tsv').write_text(
        'Sample ID\tFrequency\n#q2:types\tcategorical\n'
        'sampleA\t80.0\nsampleB\t60.0\nsampleC\t70.0\nsampleD\t50.0\n'
    )

    (out / 'biom_table').mkdir(parents=True)
    (out / 'biom_table' / 'table-with-taxonomy.biom.tsv').write_text(
        '# Constructed from biom file\n'
        '#OTU ID\tsampleA\tsampleB\tsampleC\tsampleD\ttaxonomy\n'
        'f1\t50.0\t30.0\t40.0\t20.0\tk__Fungi; p__Ascomycota; g__Fusarium\n'
        'f2\t30.0\t30.0\t30.0\t30.0\tk__Fungi; p__Basidiomycota; g__Trichosporon\n'
    )

    for metric in ('bray_curtis', 'unweighted_unifrac'):
        d = out / 'core-metrics-results' / f'{metric}_pcoa_export'
        d.mkdir(parents=True)
        (d / 'ordination.txt').write_text(
            'Eigvals\t2\n0.5\t0.03\n\nProportion explained\t2\n0.8\t0.1\n\nSpecies\t0\t0\n\n'
            'Site\t4\t2\n'
            'sampleA\t0.4\t-0.1\nsampleB\t-0.5\t-0.05\nsampleC\t0.1\t0.2\nsampleD\t0.0\t0.0\n\n'
            'Biplot\t0\t0\n\nSite constraints\t0\t0\n'
        )
        dm_dir = out / 'core-metrics-results' / f'{metric}_distance_export'
        dm_dir.mkdir(parents=True)
        (dm_dir / 'distance-matrix.tsv').write_text(
            '\tsampleA\tsampleB\tsampleC\tsampleD\n'
            'sampleA\t0.0\t0.3\t0.6\t0.7\n'
            'sampleB\t0.3\t0.0\t0.5\t0.6\n'
            'sampleC\t0.6\t0.5\t0.0\t0.2\n'
            'sampleD\t0.7\t0.6\t0.2\t0.0\n'
        )

    _write_alpha_qzv(out / 'alpha-group-significance-shannon.qzv', {
        'site': (0.9, 0.6, {'siteA (n=2)': [0.97, 1.5], 'siteB (n=2)': [1.2, 1.3]}),
    })
    _write_beta_qzv(out / 'beta-group-significance-site-bray_curtis.qzv', {
        'method name': 'PERMANOVA', 'test statistic name': 'pseudo-F',
        'sample size': '4', 'number of groups': '2', 'test statistic': '1.5', 'p-value': '0.2',
    })
    _write_beta_qzv(out / 'beta-group-significance-host-plant-unweighted_unifrac.qzv', {
        'method name': 'PERMANOVA', 'test statistic name': 'pseudo-F',
        'sample size': '4', 'number of groups': '2', 'test statistic': '0.4', 'p-value': '0.5',
    })
    _write_rarefaction_qzv(out / 'alpha-rarefaction.qzv', {
        'observed_features': {
            'sampleA': {'depth-1_iter-1': '1.0', 'depth-1_iter-2': '1.0'},
            'sampleB': {'depth-1_iter-1': '1.0', 'depth-1_iter-2': '2.0'},
        },
    })
    (out / 'sample-classifier-site').mkdir(parents=True)
    _write_classifier_qzv(out / 'sample-classifier-site' / 'accuracy_results.qzv',
                           'x\tA\tB\tOverall Accuracy\nA\t1.0\t0.0\t\nB\t0.0\t1.0\t\n'
                           'Overall Accuracy\t\t\t0.75\nBaseline Accuracy\t\t\t0.5\nAccuracy Ratio\t\t\t1.5\n')

    (out / 'run_metadata.json').write_text(json.dumps(_SAMPLE_RUN_METADATA))

    metadata_path = tmp_path / 'metadata.tsv'
    metadata_path.write_text(
        'sample-id\tsite\n#q2:types\tcategorical\n'
        'sampleA\tsiteA\nsampleB\tsiteA\nsampleC\tsiteB\nsampleD\tsiteB\n'
    )
    return out, metadata_path


class TestSplitBetaStem:
    """Regression coverage for the overlapping-text bug real report.pdf
    generation caught: the beta-group-significance table used to dump the
    whole beta-group-significance-{column}-{metric}.qzv stem into one cell,
    which overflowed its column and drew over the next one for anything but
    the shortest column names."""

    def test_splits_column_and_metric(self):
        assert report._split_beta_stem('beta-group-significance-site-bray_curtis') == \
            ('site', 'bray_curtis')

    def test_handles_hyphenated_column_names(self):
        # Metadata column names can themselves contain hyphens (e.g.
        # "host-plant"), so this can't be a naive rsplit('-', 1).
        assert report._split_beta_stem('beta-group-significance-host-plant-unweighted_unifrac') == \
            ('host-plant', 'unweighted_unifrac')

    def test_unrecognized_metric_suffix_falls_back_to_whole_remainder(self):
        column, metric = report._split_beta_stem('beta-group-significance-site-some_new_metric')
        assert metric == ''
        assert column == 'site-some_new_metric'


class TestSignificantColumns:
    def test_column_with_any_significant_metric_is_included(self):
        pvalues_by_column = {'site': [0.6, 0.01], 'host-plant': [0.5, 0.8]}
        assert report._significant_columns(pvalues_by_column) == ['site']

    def test_none_pvalues_are_ignored_not_treated_as_significant(self):
        """A test that couldn't be computed for some column/metric (None)
        must not be mistaken for a significant (low) p-value."""
        assert report._significant_columns({'site': [None, None]}) == []

    def test_no_significant_columns_returns_empty_list(self):
        assert report._significant_columns({'site': [0.6], 'host-plant': [0.8]}) == []

    def test_result_is_sorted(self):
        pvalues_by_column = {'zzz-column': [0.01], 'aaa-column': [0.02]}
        assert report._significant_columns(pvalues_by_column) == ['aaa-column', 'zzz-column']


class TestSignificanceNote:
    def test_significant_column_gets_the_significance_reason(self):
        note = report._significance_note('site', ['site'])
        # 'significant' alone is vacuous: it also appears in the
        # non-significant branch's own wording ("did not reach statistical
        # significance"), so it can't tell the two branches apart. The
        # phrase below only appears in the significant branch.
        assert 'came back statistically significant' in note
        assert 'default' not in note
        assert 'site' in note

    def test_non_significant_default_column_gets_the_default_reason(self):
        note = report._significance_note('site', [])
        assert 'default' in note
        assert 'came back statistically significant' not in note


class TestFitCellText:
    def test_short_text_passes_through_unchanged(self):
        pdf = report._ReportPDF()
        pdf.add_page()
        pdf.set_font('Helvetica', '', 8)
        assert pdf._fit_cell_text('site', 40) == 'site'

    def test_long_text_is_truncated_with_an_ellipsis_and_fits(self):
        pdf = report._ReportPDF()
        pdf.add_page()
        pdf.set_font('Helvetica', '', 8)
        long_text = 'beta-group-significance-collection-date-bray_curtis'
        fitted = pdf._fit_cell_text(long_text, 30)
        assert fitted != long_text
        assert fitted.endswith('...')
        assert pdf.get_string_width(fitted) <= 30 - 2 * pdf.c_margin


class TestColumnWidths:
    def test_short_content_does_not_fill_the_page(self):
        """Regression test: add_table_page used to always split the full
        page width evenly across columns, so a short two-column table (e.g.
        "parameter"/"value") stretched edge to edge with a wide, pointless
        gap between short values."""
        pdf = report._ReportPDF()
        pdf.add_page()
        equal_width = (pdf.w - 2 * pdf.l_margin) / 2
        widths = pdf._column_widths(['parameter', 'value'], [['max-ee', '4.0'], ['taxa', 'Fungi']])
        assert all(w < equal_width for w in widths)

    def test_column_width_never_exceeds_equal_division(self):
        """A column with genuinely long content is capped at what equal
        division across all columns would give it, so it can't crowd out
        the others or push the table past the page width."""
        pdf = report._ReportPDF()
        pdf.add_page()
        equal_width = (pdf.w - 2 * pdf.l_margin) / 2
        widths = pdf._column_widths(
            ['comparison', 'value'],
            [['a very very very long value that would otherwise stretch this column '
              'well past what an even two-column split would give it', '1.0']])
        assert widths[0] == equal_width

    def test_total_width_never_exceeds_available_page_width(self):
        pdf = report._ReportPDF()
        pdf.add_page()
        available = pdf.w - 2 * pdf.l_margin
        widths = pdf._column_widths(
            ['a', 'b', 'c'],
            [['short', 'a moderately long value here', 'x']])
        assert sum(widths) <= available


class TestPcoaFigure:
    def test_groups_get_distinct_markers_past_one_palette_pass(self):
        """The scatter series' actual marker (not just group_marker's own
        mapping) must reflect the group -- this exercises _pcoa_figure's own
        use of group_marker, not just the pure mapping function."""
        import matplotlib.pyplot as plt
        n = len(report._COLORBLIND_PALETTE) + 1
        sample_ids = [f'sample{i}' for i in range(n)]
        sample_coords = {sid: (float(i), float(i)) for i, sid in enumerate(sample_ids)}
        metadata_table = {sid: {'site': f'site{i:02d}'} for i, sid in enumerate(sample_ids)}
        fig = report._pcoa_figure(sample_coords, (0.5, 0.2), metadata_table, 'site', 'title')
        ax = fig.axes[0]
        markers = {tuple(coll.get_paths()[0].vertices.round(3).flatten()) for coll in ax.collections}
        # n groups, each its own scatter() call -> n distinct path-collections;
        # at least two different marker shapes once past one palette pass.
        assert len(ax.collections) == n
        assert len(markers) > 1
        plt.close(fig)

    def test_each_group_gets_its_own_mapped_color(self):
        """The metadata->color mapping otherwise has no direct coverage --
        e.g. every group being plotted in the same (say, first palette)
        color would pass every other test in this class."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        sample_coords = {'sampleA': (0.1, 0.1), 'sampleB': (-0.1, -0.1)}
        metadata_table = {'sampleA': {'site': 'siteA'}, 'sampleB': {'site': 'siteB'}}
        group_color = {'siteA': '#111111', 'siteB': '#222222'}
        fig = report._pcoa_figure(sample_coords, (0.5, 0.2), metadata_table, 'site', 'title',
                                   group_color=group_color)
        ax = fig.axes[0]
        colors_by_group = {coll.get_label(): tuple(coll.get_facecolor()[0]) for coll in ax.collections}
        assert colors_by_group['siteA'] == mcolors.to_rgba('#111111', alpha=0.7)
        assert colors_by_group['siteB'] == mcolors.to_rgba('#222222', alpha=0.7)
        plt.close(fig)

    def test_sample_without_a_group_value_is_labelled_missing_in_neutral_grey(self):
        """A sample with an empty value (or no metadata row) used to be
        plotted under an empty-string label, which matplotlib silently
        leaves out of the legend; it must not take a palette color either."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        sample_coords = {'sampleA': (0.1, 0.1), 'sampleB': (-0.1, -0.1), 'sampleC': (0.0, 0.2)}
        metadata_table = {'sampleA': {'site': 'siteA'}, 'sampleB': {'site': ''}}  # sampleC: no row at all
        fig = report._pcoa_figure(sample_coords, (0.5, 0.2), metadata_table, 'site', 'title')
        ax = fig.axes[0]
        legend_labels = [text.get_text() for text in ax.get_legend().get_texts()]
        assert legend_labels == ['(missing)', 'siteA']
        missing = next(coll for coll in ax.collections if coll.get_label() == '(missing)')
        assert len(missing.get_offsets()) == 2
        assert tuple(missing.get_facecolor()[0]) == mcolors.to_rgba('#666666', alpha=0.7)
        assert report._group_color_map(metadata_table, 'site') == {'siteA': report._COLORBLIND_PALETTE[0]}
        plt.close(fig)

    def test_axis_labels_use_the_matching_proportion_explained_value(self):
        """Regression-shaped coverage: PC1's label must show
        proportion_explained[0], PC2's proportion_explained[1] -- swapped
        would silently mislabel which axis explains how much variance."""
        import matplotlib.pyplot as plt
        sample_coords = {'sampleA': (0.1, 0.2), 'sampleB': (-0.1, -0.2)}
        metadata_table = {'sampleA': {'site': 'siteA'}, 'sampleB': {'site': 'siteB'}}
        fig = report._pcoa_figure(sample_coords, (0.83, 0.17), metadata_table, 'site', 'title')
        ax = fig.axes[0]
        assert '83.0%' in ax.get_xlabel()
        assert '17.0%' in ax.get_ylabel()
        plt.close(fig)


class TestDendrogramFigure:
    def test_builds_a_figure_with_one_axes(self):
        import matplotlib.pyplot as plt
        import pandas as pd
        distance_df = pd.DataFrame(
            [[0.0, 0.3, 0.6], [0.3, 0.0, 0.5], [0.6, 0.5, 0.0]],
            index=['sampleA', 'sampleB', 'sampleC'], columns=['sampleA', 'sampleB', 'sampleC'])
        metadata_table = {'sampleA': {'site': 'siteA'}, 'sampleB': {'site': 'siteA'},
                           'sampleC': {'site': 'siteB'}}
        fig = report._dendrogram_figure(distance_df, metadata_table, 'site', 'title')
        assert len(fig.axes) == 1
        plt.close(fig)

    def test_leaf_labels_include_every_sample(self):
        import matplotlib.pyplot as plt
        import pandas as pd
        distance_df = pd.DataFrame(
            [[0.0, 0.3, 0.6], [0.3, 0.0, 0.5], [0.6, 0.5, 0.0]],
            index=['sampleA', 'sampleB', 'sampleC'], columns=['sampleA', 'sampleB', 'sampleC'])
        fig = report._dendrogram_figure(distance_df, {}, None, 'title')
        labels = {t.get_text() for t in fig.axes[0].get_xticklabels()}
        assert labels == {'sampleA', 'sampleB', 'sampleC'}
        plt.close(fig)

    def test_switches_to_left_orientation_past_the_sample_threshold(self):
        import matplotlib.pyplot as plt
        import pandas as pd
        n = report._MANY_SAMPLES_THRESHOLD + 1
        ids = [f'sample{i}' for i in range(n)]
        import numpy as np
        rng = np.random.default_rng(0)
        values = rng.random((n, n))
        values = (values + values.T) / 2
        for i in range(n):
            values[i, i] = 0.0
        distance_df = pd.DataFrame(values, index=ids, columns=ids)
        fig = report._dendrogram_figure(distance_df, {}, None, 'title')
        labels = {t.get_text() for t in fig.axes[0].get_yticklabels()}
        assert labels == set(ids)
        plt.close(fig)

    def test_leaf_label_color_matches_the_samples_group(self):
        """The metadata->color mapping otherwise has no direct coverage --
        e.g. removing the tick_label.set_color(...) call entirely (leaving
        every label at matplotlib's default black) would pass every other
        test in this class."""
        import matplotlib.pyplot as plt
        import pandas as pd
        distance_df = pd.DataFrame(
            [[0.0, 0.3, 0.6], [0.3, 0.0, 0.5], [0.6, 0.5, 0.0]],
            index=['sampleA', 'sampleB', 'sampleC'], columns=['sampleA', 'sampleB', 'sampleC'])
        metadata_table = {'sampleA': {'site': 'siteA'}, 'sampleB': {'site': 'siteA'},
                           'sampleC': {'site': 'siteB'}}
        group_color = {'siteA': '#111111', 'siteB': '#222222'}
        fig = report._dendrogram_figure(distance_df, metadata_table, 'site', 'title', group_color=group_color)
        ax = fig.axes[0]
        colors_by_sample = {t.get_text(): t.get_color() for t in ax.get_xticklabels()}
        assert colors_by_sample['sampleA'] == '#111111'
        assert colors_by_sample['sampleB'] == '#111111'
        assert colors_by_sample['sampleC'] == '#222222'
        plt.close(fig)


class TestGroupColorMap:
    def test_deterministic_for_the_same_table(self):
        """build_report() relies on computing this once and passing the same
        dict to both the PCoA and dendrogram figures -- meaningless unless
        calling it twice on the same input gives back the same mapping."""
        table = {'sampleA': {'site': 'siteA'}, 'sampleB': {'site': 'siteB'}}
        assert report._group_color_map(table, 'site') == report._group_color_map(table, 'site')

    def test_groups_assigned_in_sorted_order(self):
        table = {'sampleA': {'site': 'siteB'}, 'sampleB': {'site': 'siteA'}}
        color_map = report._group_color_map(table, 'site')
        assert color_map['siteA'] == report._COLORBLIND_PALETTE[0]
        assert color_map['siteB'] == report._COLORBLIND_PALETTE[1]

    def test_no_report_column_returns_empty_map(self):
        assert report._group_color_map({'sampleA': {'site': 'siteA'}}, None) == {}


class TestGroupMarkerMap:
    def test_groups_within_one_palette_length_all_get_the_first_marker(self):
        table = {f'sample{i}': {'site': f'site{i}'} for i in range(len(report._COLORBLIND_PALETTE))}
        marker_map = report._group_marker_map(table, 'site')
        assert set(marker_map.values()) == {report._MARKER_CYCLE[0]}

    def test_a_group_past_one_full_palette_pass_gets_the_next_marker(self):
        """The whole point: once there are more groups than palette colors,
        a group that reuses an earlier group's color must get a different
        marker so the two are still visually distinguishable."""
        n = len(report._COLORBLIND_PALETTE) + 1
        table = {f'sample{i}': {'site': f'site{i:02d}'} for i in range(n)}
        marker_map = report._group_marker_map(table, 'site')
        colors = report._group_color_map(table, 'site')
        groups = sorted(f'site{i:02d}' for i in range(n))
        first_group, repeated_group = groups[0], groups[len(report._COLORBLIND_PALETTE)]
        assert colors[first_group] == colors[repeated_group]  # same color, as expected
        assert marker_map[first_group] != marker_map[repeated_group]  # but different marker

    def test_no_report_column_returns_empty_map(self):
        assert report._group_marker_map({'sampleA': {'site': 'siteA'}}, None) == {}


class TestWrapLegendLabel:
    def test_short_label_is_unchanged(self):
        assert report._wrap_legend_label('country') == 'country'

    def test_long_label_is_wrapped_onto_multiple_lines(self):
        wrapped = report._wrap_legend_label('self-heating coal mine waste material')
        assert '\n' in wrapped
        assert all(len(line) <= 18 for line in wrapped.split('\n'))

    def test_wrapping_does_not_drop_or_reorder_words(self):
        label = 'self-heating coal mine waste material'
        wrapped = report._wrap_legend_label(label)
        assert wrapped.replace('\n', ' ') == label


class TestReadableTextColor:
    def test_pale_color_is_darkened(self):
        # The old Okabe-Ito palette's yellow -- the concrete color that
        # prompted this: legible as a dot/bar fill, too pale to read as
        # small text. Not in the current palette, but still a valid case
        # for this function on any pale color, whatever palette is in use.
        darkened = report._readable_text_color('#F0E442')
        assert darkened != '#F0E442'
        r, g, b = (int(darkened[i:i + 2], 16) for i in (1, 3, 5))
        brightness = (299 * r + 587 * g + 114 * b) / 1000
        assert brightness <= 190

    def test_already_dark_color_is_unchanged(self):
        assert report._readable_text_color('#0072B2') == '#0072B2'


class TestFitWidthMm:
    def test_landscape_figure_uses_default_width(self):
        import matplotlib.pyplot as plt
        fig, _ax = plt.subplots(figsize=(6, 5))
        assert report._fit_width_mm(fig) == 180
        plt.close(fig)

    def test_tall_figure_is_narrowed_to_fit_one_page(self):
        import matplotlib.pyplot as plt
        fig, _ax = plt.subplots(figsize=(7, 11))
        width = report._fit_width_mm(fig)
        height = width * 11 / 7
        assert width < 180
        assert height <= 230
        plt.close(fig)


class TestGenusColors:
    def test_unclassified_and_other_get_distinct_fixed_neutral_greys(self):
        colors = report._genus_colors(['Unclassified', 'Fusarium', 'Penicillium', 'Other'])
        assert colors[0] == report._NEUTRAL_GREY_UNCLASSIFIED
        assert colors[3] == report._NEUTRAL_GREY_OTHER
        assert colors[0] != colors[3]  # must not look the same as each other

    def test_named_genera_cycle_the_colorblind_palette_without_using_grey(self):
        colors = report._genus_colors(['Unclassified', 'Fusarium', 'Penicillium', 'Other'])
        assert colors[1] == report._COLORBLIND_PALETTE[0]
        assert colors[2] == report._COLORBLIND_PALETTE[1]
        assert colors[1] not in (report._NEUTRAL_GREY_UNCLASSIFIED, report._NEUTRAL_GREY_OTHER)
        assert colors[2] not in (report._NEUTRAL_GREY_UNCLASSIFIED, report._NEUTRAL_GREY_OTHER)


class TestGenusBarplotFigure:
    def test_switches_to_horizontal_bars_past_the_sample_threshold(self):
        import pandas as pd
        many_samples = pd.DataFrame(
            {f'sample{i}': [0.5, 0.5] for i in range(report._MANY_SAMPLES_THRESHOLD + 1)},
            index=['Fusarium', 'Other'])
        fig = report._genus_barplot_figure(many_samples)
        ax = fig.axes[0]
        assert ax.get_xlabel() == 'Relative abundance'  # horizontal: abundance is the x-axis
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_stays_vertical_at_or_below_the_sample_threshold(self):
        import pandas as pd
        few_samples = pd.DataFrame(
            {f'sample{i}': [0.5, 0.5] for i in range(report._MANY_SAMPLES_THRESHOLD)},
            index=['Fusarium', 'Other'])
        fig = report._genus_barplot_figure(few_samples)
        ax = fig.axes[0]
        assert ax.get_ylabel() == 'Relative abundance'  # vertical: abundance is the y-axis
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestRarefactionFigure:
    def test_legend_splits_into_two_columns_past_the_sample_threshold(self):
        curves = {f'sample{i}': [(0, 1), (100, 2)] for i in range(report._MANY_SAMPLES_THRESHOLD + 1)}
        fig = report._rarefaction_figure(curves, 'observed_features')
        legend = fig.axes[0].get_legend()
        assert legend._ncols == 2
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_legend_stays_single_column_at_or_below_the_threshold(self):
        curves = {f'sample{i}': [(0, 1), (100, 2)] for i in range(report._MANY_SAMPLES_THRESHOLD)}
        fig = report._rarefaction_figure(curves, 'observed_features')
        legend = fig.axes[0].get_legend()
        assert legend._ncols == 1
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestSeqLengthHistogramFigure:
    def test_builds_a_figure_with_one_axes(self):
        import matplotlib.pyplot as plt
        fig = report._seq_length_histogram_figure([200, 210, 205, 198, 300])
        assert len(fig.axes) == 1
        assert fig.axes[0].get_xlabel() == 'Sequence length (bp)'
        plt.close(fig)


class TestConfidenceHistogramFigure:
    def test_builds_a_figure_with_one_axes(self):
        import matplotlib.pyplot as plt
        fig = report._confidence_histogram_figure([0.9, 0.95, 1.0, 0.6])
        assert len(fig.axes) == 1
        assert fig.axes[0].get_xlabel() == 'Classification confidence'
        plt.close(fig)


class TestAlphaBoxplotFigure:
    def test_boxes_use_the_shared_group_colors_not_their_position(self):
        """Regression test: boxes were colored by position among the groups
        present in the alpha results, so a group missing there (all its
        samples dropped at rarefaction) shifted every later group to a
        different color than its PCoA dots/dendrogram labels."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        alpha_results = {'shannon': {'site': {'groups': {'siteB': [1, 2], 'siteC': [3, 4]}}}}
        group_color = {'siteA': '#111111', 'siteB': '#222222', 'siteC': '#333333'}
        fig = report._alpha_boxplot_figure(alpha_results, 'site', group_color=group_color)
        face_colors = [tuple(patch.get_facecolor()) for patch in fig.axes[0].patches]
        assert face_colors == [mcolors.to_rgba('#222222', alpha=0.75), mcolors.to_rgba('#333333', alpha=0.75)]
        plt.close(fig)

    def test_four_metrics_form_a_2x2_grid_not_a_single_row(self):
        alpha_results = {
            metric: {'site': {'groups': {'siteA': [1, 2], 'siteB': [3, 4]}}}
            for metric in ('faith_pd', 'observed_features', 'shannon', 'evenness')
        }
        fig = report._alpha_boxplot_figure(alpha_results, 'site')
        visible_axes = [ax for ax in fig.axes if ax.get_visible()]
        assert len(visible_axes) == 4
        # A 2x2 grid means axes pair up on 2 distinct x-positions and 2
        # distinct y-positions in figure coordinates, not 4 of each (a 1x4 row).
        positions = [ax.get_position() for ax in visible_axes]
        assert len({round(p.x0, 3) for p in positions}) == 2
        assert len({round(p.y0, 3) for p in positions}) == 2
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_unused_grid_cells_are_hidden_not_left_as_blank_empty_axes(self):
        alpha_results = {
            metric: {'site': {'groups': {'siteA': [1, 2], 'siteB': [3, 4]}}}
            for metric in ('faith_pd', 'observed_features', 'shannon')
        }
        fig = report._alpha_boxplot_figure(alpha_results, 'site')
        visible_axes = [ax for ax in fig.axes if ax.get_visible()]
        assert len(visible_axes) == 3
        assert len(fig.axes) == 4  # 2x2 grid, one cell hidden
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_many_groups_rotate_labels_fully_vertical(self):
        """Regression test: 45-degree labels' horizontal footprint grows
        with label length, so a column with many (real example: 14)
        categories crowded adjacent rotated labels into each other and made
        them unreadable in a real generated report. Fully vertical labels
        take the same small footprint regardless of label length."""
        import matplotlib.pyplot as plt
        groups = {f'group{i}': [1, 2] for i in range(report._MANY_SAMPLES_THRESHOLD + 1)}
        alpha_results = {'shannon': {'site': {'groups': groups}}}
        fig = report._alpha_boxplot_figure(alpha_results, 'site')
        ax = [a for a in fig.axes if a.get_visible()][0]
        assert ax.xaxis.get_ticklabels()[0].get_rotation() == 90
        plt.close(fig)

    def test_few_groups_keep_the_45_degree_rotation(self):
        import matplotlib.pyplot as plt
        alpha_results = {'shannon': {'site': {'groups': {'siteA': [1, 2], 'siteB': [3, 4]}}}}
        fig = report._alpha_boxplot_figure(alpha_results, 'site')
        ax = [a for a in fig.axes if a.get_visible()][0]
        assert ax.xaxis.get_ticklabels()[0].get_rotation() == 45
        plt.close(fig)


_SAMPLE_RUN_METADATA = {
    'pipeline': {
        'qiime2_its_version': '0.2.0',
        'command_line': 'qiime2-its \\\n    -q rachis-qiime2-2026.7 \\\n    -i in \\\n    -o out \\\n'
                        '    -m meta.tsv \\\n    -c clf.qza \\\n    -pe',
        'start_time': '2026-09-16T12:00:00+00:00',
        'end_time': '2026-09-16T12:02:30+00:00',
        'duration_seconds': 150.0,
    },
    'environment': {
        'username': 'bioinfo', 'hostname': 'workstation', 'platform': 'Linux-x86_64',
        'conda_env': 'rachis-qiime2-2026.7', 'python_version': '3.12.13',
        'qiime2_framework_version': '2026.7.0', 'bbmap_version': '39.80',
    },
    'qiime2_plugins': {'dada2': '2026.7.0', 'itsxpress': '2.2.0'},
    'inputs': {
        'input_folder': '/in', 'metadata_file': 'meta.tsv', 'classifier_file': 'clf.qza',
        'output_folder': '/out',
        'samples': [{'sample_id': 'sampleA', 'files': ['sampleA_S1_L001_R1_001.fastq.gz']}],
    },
    'parameters': {'max_ee': 4.0, 'allow_one_off': True, 'input': '/in', 'qiime2': 'rachis-qiime2-2026.7'},
}


class TestRunMetadataRows:
    def test_run_info_rows_includes_key_qa_fields(self):
        rows = dict(report._run_info_rows(_SAMPLE_RUN_METADATA))
        assert rows['Run by'] == 'bioinfo@workstation'
        assert rows['QIIME2 framework version'] == '2026.7.0'
        assert rows['Run duration'] == '2m30s'
        assert 'qiime2-its' in rows['Command invoked']
        assert rows['BBMap (bbduk.sh) version'] == '39.80'

    def test_run_duration_omits_zero_leading_units(self):
        """A sub-hour run shouldn't show '0d0h' -- and a multi-day one should
        show real days, not an absurd number of minutes."""
        under_an_hour = {**_SAMPLE_RUN_METADATA,
                          'pipeline': {**_SAMPLE_RUN_METADATA['pipeline'], 'duration_seconds': 45.0}}
        rows = dict(report._run_info_rows(under_an_hour))
        assert rows['Run duration'] == '45s'

        over_a_day = {**_SAMPLE_RUN_METADATA,
                       'pipeline': {**_SAMPLE_RUN_METADATA['pipeline'], 'duration_seconds': 90000.0}}
        rows = dict(report._run_info_rows(over_a_day))
        assert rows['Run duration'] == '1d1h'

    def test_run_info_rows_reports_bbmap_not_installed(self):
        metadata = {**_SAMPLE_RUN_METADATA,
                     'environment': {**_SAMPLE_RUN_METADATA['environment'], 'bbmap_version': None}}
        rows = dict(report._run_info_rows(metadata))
        assert rows['BBMap (bbduk.sh) version'] == 'not installed'

    def test_parameter_rows_excludes_fields_shown_on_run_info_page(self):
        rows = report._parameter_rows(_SAMPLE_RUN_METADATA)
        keys = [row[0] for row in rows]
        assert 'input' not in keys
        assert 'qiime2' not in keys
        assert ['max-ee', 4.0] in rows
        assert ['allow-one-off', True] in rows

    def test_sample_file_rows_one_row_per_file(self):
        assert report._sample_file_rows(_SAMPLE_RUN_METADATA) == [['sampleA', 'sampleA_S1_L001_R1_001.fastq.gz']]

    def test_plugin_rows_sorted(self):
        assert report._plugin_rows(_SAMPLE_RUN_METADATA) == [['dada2', '2026.7.0'], ['itsxpress', '2.2.0']]


class TestBuildReport:
    def test_produces_a_valid_multi_page_pdf_with_every_expected_section(self, tmp_path, mocker):
        output_folder, metadata_path = _build_synthetic_output_folder(tmp_path)
        titles = _page_titles(mocker)

        report_path = report.build_report(output_folder, metadata_path)

        assert report_path == output_folder / 'report.pdf'
        assert report_path.read_bytes()[:4] == b'%PDF'
        # One title per section the synthetic fixture provides data for --
        # a byte-size/page-count proxy can't tell "every section rendered"
        # apart from "one section silently never ran but the PDF is still
        # plausibly sized" (e.g. a typo'd qzv glob, or a dropped `if
        # classifier_rows:` block).
        for expected in ('DADA2 read retention',
                          'Alpha diversity group significance (Kruskal-Wallis)',
                          'Alpha diversity by site',
                          'Beta diversity group significance (PERMANOVA)',
                          'site: Bray-Curtis PCoA',
                          'site: Bray-Curtis sample clustering',
                          'Genus-level relative abundance',
                          'Rarefaction curve',
                          'Sample classifier results'):
            assert expected in titles, f'missing page: {expected!r}'

    def test_column_names_sanitized_for_file_names_are_mapped_back(self, tmp_path, mocker):
        """Regression test: the pipeline names its per-column outputs after
        safe_filename_component(column), and q2-diversity URL-quotes the
        column in its alpha jsonp file names, so "Host Plant" came back as
        "Host_Plant" / "Host%20Plant" -- matching nothing in the metadata
        (every sample lost its group) and printed mangled in the tables."""
        output_folder, metadata_path = _build_synthetic_output_folder(tmp_path)
        metadata_path.write_text(
            'sample-id\tsite\tHost Plant\n#q2:types\tcategorical\tcategorical\n'
            'sampleA\tsiteA\toak\nsampleB\tsiteA\tpine\nsampleC\tsiteB\toak\nsampleD\tsiteB\tpine\n')
        _write_beta_qzv(output_folder / 'beta-group-significance-Host_Plant-bray_curtis.qzv', {
            'method name': 'PERMANOVA', 'test statistic name': 'pseudo-F',
            'sample size': '4', 'number of groups': '2', 'test statistic': '9.5', 'p-value': '0.01',
        })
        _write_alpha_qzv(output_folder / 'alpha-group-significance-shannon.qzv', {
            'Host%20Plant': (5.0, 0.01, {'oak (n=2)': [0.9, 1.5], 'pine (n=2)': [2.2, 2.3]}),
        })
        (output_folder / 'sample-classifier-Host_Plant').mkdir()
        _write_classifier_qzv(output_folder / 'sample-classifier-Host_Plant' / 'accuracy_results.qzv',
                               'Overall Accuracy\t\t\t0.75\n')
        titles = _page_titles(mocker)
        tables = mocker.spy(report._ReportPDF, 'add_table_page')
        pcoa = mocker.spy(report, '_pcoa_figure')

        report.build_report(output_folder, metadata_path)

        assert 'Host Plant: Bray-Curtis PCoA' in titles
        assert 'Alpha diversity by Host Plant' in titles
        assert not any('Host_Plant' in title or 'Host%20Plant' in title for title in titles)
        host_plant_calls = [call for call in pcoa.call_args_list if call.args[3] == 'Host Plant']
        assert host_plant_calls and host_plant_calls[0].kwargs['group_color'].keys() == {'oak', 'pine'}
        table_cells = {str(cell) for call in tables.call_args_list for row in call.args[3] for cell in row}
        assert 'Host Plant' in table_cells
        assert 'Host_Plant' not in table_cells and 'Host%20Plant' not in table_cells

    def test_auto_picks_first_eligible_column_when_not_specified(self, tmp_path, mocker):
        output_folder, metadata_path = _build_synthetic_output_folder(tmp_path)
        titles = _page_titles(mocker)

        report_path = report.build_report(output_folder, metadata_path, report_column=None)

        assert report_path.exists()
        # Not just "didn't raise" -- 'site' (the fixture's only eligible
        # column) must actually have been picked and driven real pages,
        # not silently landed on report_column=None (which would still
        # produce a same-sized-ish PDF, just missing every site-grouped page).
        assert 'Alpha diversity by site' in titles
        assert 'site: Bray-Curtis PCoA' in titles

    def test_invalid_report_column_falls_back_with_a_warning_instead_of_degrading_silently(self, tmp_path, capsys):
        """Regression test: a typo'd/nonexistent --report-metadata-column
        used to fail silently -- every group-by-report_column lookup just
        returns nothing for a column that doesn't exist, producing a
        degraded report (blank groups) with no indication anything was
        wrong. It should fall back to auto-selection instead, visibly."""
        output_folder, metadata_path = _build_synthetic_output_folder(tmp_path)

        report_path = report.build_report(output_folder, metadata_path, report_column='not-a-real-column')

        assert report_path.exists()
        assert 'not-a-real-column' in capsys.readouterr().out

    def test_ineligible_report_column_falls_back_with_a_warning(self, tmp_path, capsys):
        """Regression test: a real column that exists in the metadata file
        but isn't eligible for grouping (here, numeric rather than
        categorical) used to pass the old "does this column exist at all"
        check and reach the PCoA/dendrogram loop -- which only checks that
        the per-metric ordination/distance artifacts exist, not that any
        significance test was ever run for this particular column -- so the
        figures got built and captioned as "did not reach statistical
        significance" even though no such test exists for a numeric column.
        Must fall back to auto-selection instead, same as a nonexistent
        column."""
        output_folder, metadata_path = _build_synthetic_output_folder(tmp_path)
        metadata_path.write_text(
            'sample-id\tsite\televation\n#q2:types\tcategorical\tnumeric\n'
            'sampleA\tsiteA\t100\nsampleB\tsiteA\t150\nsampleC\tsiteB\t200\nsampleD\tsiteB\t250\n'
        )

        report_path = report.build_report(output_folder, metadata_path, report_column='elevation')

        assert report_path.exists()
        out = capsys.readouterr().out
        assert 'elevation' in out
        assert 'eligible' in out

    def test_significant_non_default_column_adds_its_own_pages(self, tmp_path, mocker):
        """A metadata column other than the default/auto-picked one, but
        that comes back statistically significant, should get its own
        PCoA/dendrogram pages too -- not just the default column's."""
        output_folder, metadata_path = _build_synthetic_output_folder(tmp_path)

        # Add a second metadata column with a genuinely significant result,
        # distinct from 'site' (the auto-picked default for this fixture).
        metadata_path.write_text(
            'sample-id\tsite\ttreatment\n#q2:types\tcategorical\tcategorical\n'
            'sampleA\tsiteA\ttreatA\nsampleB\tsiteA\ttreatA\n'
            'sampleC\tsiteB\ttreatB\nsampleD\tsiteB\ttreatB\n'
        )
        _write_beta_qzv(output_folder / 'beta-group-significance-treatment-bray_curtis.qzv', {
            'method name': 'PERMANOVA', 'test statistic name': 'pseudo-F',
            'sample size': '4', 'number of groups': '2', 'test statistic': '9.0', 'p-value': '0.01',
        })
        titles = _page_titles(mocker)

        report_path = report.build_report(output_folder, metadata_path)

        assert report_path.exists()
        # The actual claim the docstring makes -- not a PDF-size proxy,
        # which a broken `beta_columns = sorted({report_column})` (i.e.
        # significant columns never added) doesn't move enough to notice:
        # the extra PERMANOVA table row alone grows the PDF regardless of
        # whether treatment's own PCoA/dendrogram pages were ever built.
        assert 'treatment: Bray-Curtis PCoA' in titles
        assert 'treatment: Bray-Curtis sample clustering' in titles
        assert 'site: Bray-Curtis PCoA' in titles  # the default column's pages still there too

    def test_works_with_missing_optional_artifacts(self, tmp_path):
        """A --skip-advanced-stats run won't have any of the group-
        significance/classifier/pcoa-export artifacts -- the report should
        degrade gracefully to just the pages it has data for, not crash."""
        out = tmp_path / 'output'
        (out / 'dada2_stats').mkdir(parents=True)
        (out / 'dada2_stats' / 'stats.tsv').write_text(
            'sample-id\tinput\n#q2:types\tnumeric\nsampleA\t100\n'
        )
        metadata_path = tmp_path / 'metadata.tsv'
        metadata_path.write_text('sample-id\tsite\nsampleA\tsiteA\n')

        report_path = report.build_report(out, metadata_path)

        assert report_path.exists()
        assert report_path.read_bytes()[:4] == b'%PDF'
