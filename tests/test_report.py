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


def _write_alpha_qzv(path, column_stats):
    with zipfile.ZipFile(path, 'w') as zf:
        for column, (h, p, groups) in column_stats.items():
            group_data = json.dumps({'name': None, 'index': list(groups), 'data': list(groups.values())})
            content = (f"load_data('{column}',{group_data},{{}},"
                       f'{{"H": {h}, "p": {p}}},\'<table></table>\',\'x.csv\', \'metric\');')
            zf.writestr(f'uuid1/data/column-{column}.jsonp', content)


def _write_beta_qzv(path, overview_rows):
    html = '<html><body><table>' + ''.join(
        f'<th>{k}</th><td>{v}</td>' for k, v in overview_rows.items()) + '</table></body></html>'
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
        assert rows['Run duration'] == '2.5 min'
        assert 'qiime2-its' in rows['Command invoked']
        assert rows['BBMap (bbduk.sh) version'] == '39.80'

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
    def test_produces_a_valid_multi_page_pdf(self, tmp_path):
        output_folder, metadata_path = _build_synthetic_output_folder(tmp_path)

        report_path = report.build_report(output_folder, metadata_path)

        assert report_path == output_folder / 'report.pdf'
        assert report_path.stat().st_size > 1000
        assert report_path.read_bytes()[:4] == b'%PDF'

    def test_auto_picks_first_eligible_column_when_not_specified(self, tmp_path):
        output_folder, metadata_path = _build_synthetic_output_folder(tmp_path)
        # Should not raise, and should pick up 'site' automatically.
        report_path = report.build_report(output_folder, metadata_path, report_column=None)
        assert report_path.exists()

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
