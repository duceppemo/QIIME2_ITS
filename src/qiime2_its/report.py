"""Assembles a lightweight PDF summary of a qiime2-its run from its
already-produced output: run/sample summary, DADA2 retention, alpha
diversity (boxplots + group-significance p-values), beta diversity (PCoA +
PERMANOVA), a genus-level composition chart, the rarefaction curve, and a
sample-classifier accuracy summary if any column was successfully modeled.

Pure PDF assembly -- all the parsing this depends on lives in report_data.py
so it can be unit-tested without matplotlib/fpdf2 involved.
"""
import importlib.resources
import io
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402 (backend must be set before this import)
import pandas as pd  # noqa: E402
from fpdf import FPDF  # noqa: E402

from qiime2_its import metadata_utils, report_data  # noqa: E402


def _fig_to_png_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return buf


def _logo_bytes():
    try:
        data = importlib.resources.files('qiime2_its').joinpath('assets/logo.png').read_bytes()
    except (FileNotFoundError, ModuleNotFoundError):
        return None
    return io.BytesIO(data)


class _ReportPDF(FPDF):
    def section_title(self, title):
        self.set_font('Helvetica', 'B', 14)
        self.cell(0, 10, title, new_x='LMARGIN', new_y='NEXT')
        self.ln(1)

    def add_title_page(self, title, lines, logo=None):
        self.add_page()
        if logo is not None:
            logo_width = 56
            self.image(logo, x=(self.w - logo_width) / 2, w=logo_width)
            self.ln(4)
        self.set_font('Helvetica', 'B', 20)
        self.cell(0, 14, title, new_x='LMARGIN', new_y='NEXT')
        self.ln(6)
        self.set_font('Helvetica', '', 11)
        for line in lines:
            # multi_cell() doesn't reset x to the left margin before
            # rendering the way cell(new_x='LMARGIN') does, only after -- so
            # a run of them drifts x rightward until there's no room left.
            self.set_x(self.l_margin)
            self.multi_cell(0, 7, line)

    def _fit_cell_text(self, text, width):
        # cell() doesn't wrap or clip -- text wider than the cell just draws
        # past its border and overlaps whatever's in the next cell. Truncate
        # with an ellipsis instead, so a longer-than-expected value (a long
        # metadata column name, an unexpectedly long sample id, ...) stays
        # inside its own cell rather than corrupting its neighbor.
        # '...' rather than the Unicode ellipsis: the base Helvetica font
        # (WinAnsi-encoded) can't represent U+2026.
        ellipsis = '...'
        pad = 2 * self.c_margin
        if self.get_string_width(text) <= width - pad:
            return text
        while text and self.get_string_width(text + ellipsis) > width - pad:
            text = text[:-1]
        return (text + ellipsis) if text else ellipsis

    def add_table_page(self, title, header, rows):
        # Wide tables (many columns, e.g. DADA2 retention) get cramped in
        # portrait; landscape gives them roughly 40% more width.
        orientation = 'L' if len(header) > 5 else 'P'
        self.add_page(orientation=orientation)
        self.section_title(title)
        if not rows:
            self.set_font('Helvetica', '', 10)
            self.cell(0, 8, '(no data)', new_x='LMARGIN', new_y='NEXT')
            return
        n_cols = len(header)
        col_width = (self.w - 2 * self.l_margin) / n_cols
        self.set_font('Helvetica', 'B', 8)
        for h in header:
            self.cell(col_width, 7, self._fit_cell_text(str(h), col_width), border=1)
        self.ln()
        self.set_font('Helvetica', '', 8)
        for row in rows:
            for value in row:
                self.cell(col_width, 6, self._fit_cell_text(str(value), col_width), border=1)
            self.ln()

    def add_keyvalue_page(self, title, rows, label_width=45):
        """A label/value list, one entry per line -- for QA/provenance
        fields whose values vary too much in length for add_table_page's
        fixed-width, truncate-if-too-long cells (a full file path or command
        line must stay intact, not get an ellipsis)."""
        self.add_page()
        self.section_title(title)
        for label, value in rows:
            self.set_font('Helvetica', 'B', 9)
            self.set_x(self.l_margin)
            self.cell(label_width, 6, f'{label}:', new_x='RIGHT', new_y='TOP')
            self.set_font('Helvetica', '', 9)
            self.multi_cell(0, 6, str(value))

    def add_figure_page(self, title, fig, width=180):
        self.add_page()
        self.section_title(title)
        self.image(_fig_to_png_bytes(fig), w=width)


_DADA2_COLUMN_LABELS = {
    'percentage of input passed filter': '% passed filter',
    'percentage of input merged': '% merged',
    'percentage of input non-chimeric': '% non-chimeric',
}

# Matches the beta-group-significance-{column}-{metric}.qzv naming pipeline.py builds
# (cli/pipeline.py's `_run_advanced_stats`). Listed explicitly rather than rsplit() on
# '-' because metadata column names can themselves contain hyphens (e.g. "host-plant").
_BETA_DISTANCE_METRICS = ('bray_curtis', 'unweighted_unifrac')


def _split_beta_stem(stem):
    stem = stem.removeprefix('beta-group-significance-')
    for metric in _BETA_DISTANCE_METRICS:
        suffix = f'-{metric}'
        if stem.endswith(suffix):
            return stem[:-len(suffix)], metric
    return stem, ''


# Already shown as their own rows on the "Run information" page -- excluded
# from "Pipeline parameters" so the two pages don't just repeat each other.
_PARAMETERS_SHOWN_ELSEWHERE = {'input', 'output', 'metadata', 'classifier', 'qiime2'}


def _run_info_rows(run_metadata):
    pipeline = run_metadata.get('pipeline', {})
    env = run_metadata.get('environment', {})
    inputs = run_metadata.get('inputs', {})
    duration = pipeline.get('duration_seconds')
    return [
        ('Run started', pipeline.get('start_time', '')),
        ('Run finished', pipeline.get('end_time', '')),
        ('Run duration', f'{duration / 60:.1f} min' if duration is not None else ''),
        ('Run by', f"{env.get('username', '')}@{env.get('hostname', '')}"),
        ('Platform', env.get('platform', '')),
        ('Conda environment', env.get('conda_env', '')),
        ('qiime2-its version', pipeline.get('qiime2_its_version', '')),
        ('QIIME2 framework version', env.get('qiime2_framework_version', '')),
        ('BBMap (bbduk.sh) version', env.get('bbmap_version') or 'not installed'),
        ('Python version', env.get('python_version', '')),
        ('Input folder', inputs.get('input_folder', '')),
        ('Metadata file', inputs.get('metadata_file', '')),
        ('Classifier file', inputs.get('classifier_file', '')),
        ('Output folder', inputs.get('output_folder', '')),
        ('Command invoked', pipeline.get('command_line', '')),
    ]


def _parameter_rows(run_metadata):
    parameters = run_metadata.get('parameters', {})
    return [[key.replace('_', '-'), value] for key, value in sorted(parameters.items())
            if key not in _PARAMETERS_SHOWN_ELSEWHERE]


def _sample_file_rows(run_metadata):
    rows = []
    for entry in run_metadata.get('inputs', {}).get('samples', []):
        for filename in entry.get('files', []):
            rows.append([entry.get('sample_id', ''), filename])
    return rows


def _plugin_rows(run_metadata):
    return [[name, version] for name, version in sorted(run_metadata.get('qiime2_plugins', {}).items())]


def _dada2_summary_table(output_folder):
    stats_path = output_folder / 'dada2_stats' / 'stats.tsv'
    if not stats_path.exists():
        return None
    df = report_data.parse_dada2_stats(stats_path)
    header = ['sample'] + [_DADA2_COLUMN_LABELS.get(c, c) for c in df.columns]
    rows = []
    for sample_id, row in df.iterrows():
        rows.append([sample_id] + [f'{v:.1f}' if pd.notna(v) else '' for v in row])
    return df, header, rows


def _alpha_boxplot_figure(alpha_results, report_column):
    metrics = sorted(alpha_results)
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4), squeeze=False)
    for ax, metric in zip(axes[0], metrics):
        groups = alpha_results[metric].get(report_column, {}).get('groups', {})
        labels = sorted(groups)
        data = [groups[label] for label in labels]
        if data:
            ax.boxplot(data, tick_labels=labels)
        ax.set_title(metric)
        ax.tick_params(axis='x', rotation=45)
    fig.suptitle(f'Alpha diversity by {report_column}')
    fig.tight_layout()
    return fig


def _pcoa_figure(sample_coords, proportion_explained, metadata_table, report_column, title):
    fig, ax = plt.subplots(figsize=(6, 5))
    groups = sorted({metadata_table.get(sid, {}).get(report_column, '?') for sid in sample_coords})
    for group in groups:
        xs, ys = [], []
        for sid, (x, y) in sample_coords.items():
            if metadata_table.get(sid, {}).get(report_column, '?') == group:
                xs.append(x)
                ys.append(y)
        ax.scatter(xs, ys, label=group)
    ax.set_xlabel(f'PC1 ({proportion_explained[0] * 100:.1f}%)')
    ax.set_ylabel(f'PC2 ({proportion_explained[1] * 100:.1f}%)')
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def _genus_barplot_figure(genus_table):
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(genus_table.columns)), 5))
    genus_table.T.plot(kind='bar', stacked=True, ax=ax, legend=True)
    ax.set_ylabel('Relative abundance')
    ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()
    return fig


def _rarefaction_figure(curves, metric):
    fig, ax = plt.subplots(figsize=(6, 5))
    for sample_id, points in curves.items():
        if points:
            ax.plot([p[0] for p in points], [p[1] for p in points], marker='o', markersize=3,
                    label=sample_id)
    ax.set_xlabel('Sequencing depth')
    ax.set_ylabel(metric)
    ax.set_title(f'Rarefaction curve ({metric})')
    ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()
    return fig


def build_report(output_folder, metadata_file, report_column=None):
    """Build `<output_folder>/report.pdf`. Returns its path.

    Discovers what to include by scanning `output_folder` for the artifacts
    the pipeline may have produced (alpha-group-significance-*.qzv,
    beta-group-significance-*.qzv, etc.) rather than requiring an explicit
    list, so it degrades gracefully if `--skip-advanced-stats` was used.
    """
    output_folder = Path(output_folder)
    metadata_file = Path(metadata_file)
    metadata_table = metadata_utils.read_metadata_table(metadata_file)

    sample_freq_path = output_folder / 'sample_frequencies' / 'metadata.tsv'
    sample_frequencies = (report_data.parse_sample_frequencies(sample_freq_path)
                           if sample_freq_path.exists() else {})
    final_sample_ids = [sid for sid, freq in sample_frequencies.items() if freq > 0] \
        or list(metadata_table)

    if report_column is None:
        eligible = metadata_utils.eligible_categorical_columns(metadata_file, final_sample_ids)
        report_column = eligible[0] if eligible else None

    pdf = _ReportPDF()
    pdf.set_auto_page_break(auto=True, margin=15)

    # 1. Title / summary
    summary_lines = [f'Samples: {len(final_sample_ids)}', f'Report grouping column: {report_column or "(none eligible)"}']
    dada2 = _dada2_summary_table(output_folder)
    if dada2 is not None:
        df, _header, _rows = dada2
        if 'non-chimeric' in df.columns:
            summary_lines.append(f'Median non-chimeric reads/sample: {df["non-chimeric"].median():.0f}')
    pdf.add_title_page('QIIME2-ITS run summary', summary_lines, logo=_logo_bytes())

    # 1b. Run provenance / QA -- absent for a report built from an output
    # folder that predates this (or a run that crashed before it was written).
    run_metadata = report_data.parse_run_metadata(output_folder / 'run_metadata.json')
    if run_metadata is not None:
        pdf.add_keyvalue_page('Run information', _run_info_rows(run_metadata))
        parameter_rows = _parameter_rows(run_metadata)
        if parameter_rows:
            pdf.add_table_page('Pipeline parameters', ['parameter', 'value'], parameter_rows)
        sample_rows = _sample_file_rows(run_metadata)
        if sample_rows:
            pdf.add_table_page('Input sample files', ['sample', 'file'], sample_rows)
        plugin_rows = _plugin_rows(run_metadata)
        if plugin_rows:
            pdf.add_table_page('Installed QIIME2 plugins', ['plugin', 'version'], plugin_rows)

    # 2. DADA2 retention table
    if dada2 is not None:
        _df, header, rows = dada2
        pdf.add_table_page('DADA2 read retention', header, rows)

    # 3. Alpha diversity: group-significance p-values + boxplots
    alpha_results = {}
    for qzv_path in sorted(output_folder.glob('alpha-group-significance-*.qzv')):
        metric = qzv_path.stem.removeprefix('alpha-group-significance-')
        alpha_results[metric] = report_data.parse_alpha_group_significance(qzv_path)

    if alpha_results:
        rows = []
        for metric, by_column in alpha_results.items():
            for column, stats in by_column.items():
                rows.append([metric, column, f'{stats["h_statistic"]:.3f}', f'{stats["p_value"]:.3f}'])
        pdf.add_table_page('Alpha diversity group significance (Kruskal-Wallis)',
                            ['metric', 'metadata column', 'H', 'p-value'], rows)

        if report_column and any(report_column in by_column for by_column in alpha_results.values()):
            pdf.add_figure_page('Alpha diversity by group', _alpha_boxplot_figure(alpha_results, report_column))

    # 4. Beta diversity: PCoA + PERMANOVA
    beta_rows = []
    for qzv_path in sorted(output_folder.glob('beta-group-significance-*.qzv')):
        stats = report_data.parse_beta_group_significance(qzv_path)
        column, metric = _split_beta_stem(qzv_path.stem)
        beta_rows.append([column, metric, stats.get('test_statistic_name'),
                           f'{stats["test_statistic"]:.3f}' if stats.get('test_statistic') is not None else '',
                           f'{stats["p_value"]:.3f}' if stats.get('p_value') is not None else ''])
    if beta_rows:
        pdf.add_table_page('Beta diversity group significance (PERMANOVA)',
                            ['metadata column', 'distance metric', 'statistic', 'value', 'p-value'], beta_rows)

    if report_column:
        for metric in ('bray_curtis', 'unweighted_unifrac'):
            ordination_path = output_folder / 'core-metrics-results' / f'{metric}_pcoa_export' / 'ordination.txt'
            if ordination_path.exists():
                sample_coords, proportion_explained = report_data.parse_ordination(ordination_path)
                pdf.add_figure_page(f'{metric} PCoA', _pcoa_figure(
                    sample_coords, proportion_explained, metadata_table, report_column,
                    f'{metric} (colored by {report_column})'))

    # 5. Genus-level composition
    biom_taxo_path = output_folder / 'biom_table' / 'table-with-taxonomy.biom.tsv'
    if biom_taxo_path.exists():
        genus_table = report_data.build_genus_abundance_table(biom_taxo_path)
        pdf.add_figure_page('Genus-level relative abundance', _genus_barplot_figure(genus_table))

    # 6. Rarefaction curve
    rarefaction_qzv = output_folder / 'alpha-rarefaction.qzv'
    if rarefaction_qzv.exists():
        curves = report_data.parse_rarefaction_curve(rarefaction_qzv, 'observed_features')
        pdf.add_figure_page('Rarefaction curve', _rarefaction_figure(curves, 'observed_features'))

    # 7. Sample classifier accuracy, if any
    classifier_rows = []
    for accuracy_qzv in sorted(output_folder.glob('sample-classifier-*/accuracy_results.qzv')):
        column = accuracy_qzv.parent.name.removeprefix('sample-classifier-')
        accuracy = report_data.parse_classifier_accuracy(accuracy_qzv)
        for key, label in (('overall_accuracy', 'Overall accuracy'),
                            ('baseline_accuracy', 'Baseline accuracy'),
                            ('accuracy_ratio', 'Accuracy ratio')):
            if key in accuracy:
                classifier_rows.append([column, label, f'{accuracy[key]:.3f}'])
    if classifier_rows:
        pdf.add_table_page('Sample classifier results', ['metadata column', 'metric', 'value'],
                            classifier_rows)

    report_path = output_folder / 'report.pdf'
    pdf.output(str(report_path))
    return report_path
