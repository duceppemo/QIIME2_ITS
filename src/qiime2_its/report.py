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
            logo_width = 28
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
            self.cell(col_width, 7, str(h), border=1)
        self.ln()
        self.set_font('Helvetica', '', 8)
        for row in rows:
            for value in row:
                self.cell(col_width, 6, str(value), border=1)
            self.ln()

    def add_figure_page(self, title, fig, width=180):
        self.add_page()
        self.section_title(title)
        self.image(_fig_to_png_bytes(fig), w=width)


_DADA2_COLUMN_LABELS = {
    'percentage of input passed filter': '% passed filter',
    'percentage of input merged': '% merged',
    'percentage of input non-chimeric': '% non-chimeric',
}


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
        beta_rows.append([qzv_path.stem, stats.get('test_statistic_name'),
                           f'{stats["test_statistic"]:.3f}' if stats.get('test_statistic') is not None else '',
                           f'{stats["p_value"]:.3f}' if stats.get('p_value') is not None else ''])
    if beta_rows:
        pdf.add_table_page('Beta diversity group significance (PERMANOVA)',
                            ['comparison', 'statistic', 'value', 'p-value'], beta_rows)

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
