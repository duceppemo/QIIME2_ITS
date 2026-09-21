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
import math
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402 (backend must be set before this import)
from matplotlib.lines import Line2D  # noqa: E402
import pandas as pd  # noqa: E402
from fpdf import FPDF  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402

from qiime2_its import metadata_utils, report_data, timing  # noqa: E402

# CARTOColors "Safe" qualitative scheme (https://github.com/CartoDB/CartoColor,
# tagged "colorblind" in its own source), derived from Paul Tol's colorblind-
# safe schemes. Confirmed against CartoColor's own source (src/carto.ts) --
# its "Safe" palette's largest defined set is 11 informative hues (its 12th
# color, #888888, is a fixed grey always appended as a separate "other"
# indicator regardless of how many hues are requested, not a 12th rotating
# category color, so it's excluded here in favor of this report's own
# _NEUTRAL_GREY_* colors for that role). Replaces the 8-color Okabe-Ito
# palette used earlier in this report's design: real report columns have
# run past 8 groups, and more colorblind-safe hues before any repeat back
# to a color already used is a direct improvement for that case. Used for
# every categorical color in this report instead of matplotlib's default
# cycle.
_COLORBLIND_PALETTE = [
    '#88CCEE',  # cyan
    '#CC6677',  # rose
    '#DDCC77',  # sand
    '#117733',  # green
    '#332288',  # indigo
    '#AA4499',  # purple
    '#44AA99',  # teal
    '#999933',  # olive
    '#882255',  # wine
    '#661100',  # dark red
    '#6699CC',  # blue
]
# Fixed, non-cycled colors for a chart's "Other"/"Unclassified" catch-all
# categories, so neither ever coincides with a real category's color and both
# read consistently as "not a specific answer" the way grey conventionally
# does. Two different shades (rather than one grey for both) keep the two
# catch-alls -- which mean different things and can each be a large segment
# -- visually distinguishable from each other. Deliberately darker than a
# typical "muted" grey (#999999): checked against a real rendered report,
# that lighter grey read as visually indistinguishable from the page's white
# background, making a real (and often large, for "Unclassified") segment
# look like a gap in the bar instead.
_NEUTRAL_GREY_UNCLASSIFIED = '#4D4D4D'
_NEUTRAL_GREY_OTHER = '#A6A6A6'

plt.rcParams.update({
    'axes.prop_cycle': plt.cycler(color=_COLORBLIND_PALETTE),
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.5,
    'figure.dpi': 150,
})


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

    def add_intro_text(self, text):
        """A short explanatory paragraph under a section title -- what the
        analysis on this page actually shows, set apart from the data itself
        (italic, grey) the way a figure caption or methods note would be."""
        self.set_font('Helvetica', 'I', 9)
        self.set_text_color(90, 90, 90)
        self.set_x(self.l_margin)
        self.multi_cell(0, 5, text)
        self.set_text_color(0, 0, 0)
        self.ln(2)

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

    def _column_widths(self, header, rows):
        """Size each column to its own content (header + every cell in that
        column) instead of always splitting the page width evenly -- a
        short two-column table (e.g. "parameter"/"value") otherwise stretches
        across the full page with a wide gap between short values. Capped at
        what even division would give the column, so no single long column
        can crowd the others out or push the table past the page width."""
        equal_width = (self.w - 2 * self.l_margin) / len(header)
        pad = 2 * self.c_margin + 2
        col_widths = []
        for col_idx, h in enumerate(header):
            self.set_font('Helvetica', 'B', 9)
            natural = self.get_string_width(str(h))
            self.set_font('Helvetica', '', 9)
            for row in rows:
                natural = max(natural, self.get_string_width(str(row[col_idx])))
            col_widths.append(min(equal_width, natural + pad))
        return col_widths

    def add_table_page(self, title, header, rows, intro=None):
        # Wide tables (many columns, e.g. DADA2 retention) get cramped in
        # portrait; landscape gives them roughly 40% more width.
        orientation = 'L' if len(header) > 5 else 'P'
        self.add_page(orientation=orientation)
        self.section_title(title)
        if intro:
            self.add_intro_text(intro)
        if not rows:
            self.set_font('Helvetica', '', 10)
            self.cell(0, 8, '(no data)', new_x='LMARGIN', new_y='NEXT')
            return

        # A classic academic "three-line" table (rule / header / rule ...
        # rule) instead of a full grid: no vertical borders at all, subtle
        # zebra striping on the data rows instead (this report's tables can
        # run to 50+ rows, where striping matters more for readability than
        # it does for a handful).
        col_widths = self._column_widths(header, rows)
        table_width = sum(col_widths)
        x_start = self.l_margin

        self.set_draw_color(50, 50, 50)
        self.set_line_width(0.4)
        self.line(x_start, self.get_y(), x_start + table_width, self.get_y())

        self.set_font('Helvetica', 'B', 9)
        self.set_x(x_start)
        for h, w in zip(header, col_widths):
            self.cell(w, 7.5, self._fit_cell_text(str(h), w), border=0)
        self.ln()

        self.set_line_width(0.25)
        self.line(x_start, self.get_y(), x_start + table_width, self.get_y())

        self.set_font('Helvetica', '', 9)
        self.set_fill_color(242, 242, 242)
        for i, row in enumerate(rows):
            self.set_x(x_start)
            fill = (i % 2 == 1)
            for value, w in zip(row, col_widths):
                self.cell(w, 6.5, self._fit_cell_text(str(value), w), border=0, fill=fill)
            self.ln()

        self.set_line_width(0.4)
        self.line(x_start, self.get_y(), x_start + table_width, self.get_y())
        self.ln(3)

    def add_keyvalue_page(self, title, rows, label_width=None, intro=None):
        """A label/value list, one entry per line -- for QA/provenance
        fields whose values vary too much in length for add_table_page's
        fixed-width, truncate-if-too-long cells (a full file path or command
        line must stay intact, not get an ellipsis)."""
        self.add_page()
        self.section_title(title)
        if intro:
            self.add_intro_text(intro)
        self.set_font('Helvetica', 'B', 10)
        if label_width is None:
            # cell() doesn't wrap or clip -- a label wider than a fixed
            # label_width would bleed into the value column instead of
            # being cut off, so size the column to the widest label
            # actually present rather than guessing a fixed value.
            label_width = max((self.get_string_width(f'{label}:') for label, _ in rows), default=0) + 3
        for label, value in rows:
            value = str(value)
            self.set_font('Helvetica', 'B', 10)
            self.set_x(self.l_margin)
            self.cell(label_width, 6.5, f'{label}:', new_x='RIGHT', new_y='TOP')
            # A value with embedded newlines (e.g. the multi-line "Command
            # invoked") reads better as a monospaced block than justified
            # proportional text.
            if '\n' in value:
                self.set_font('Courier', '', 8.5)
            else:
                self.set_font('Helvetica', '', 10)
            self.multi_cell(0, 6.5, value)

    def add_figure_page(self, title, fig, width=180, intro=None):
        self.add_page()
        self.section_title(title)
        if intro:
            self.add_intro_text(intro)
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

_METRIC_DISPLAY_NAMES = {
    'bray_curtis': 'Bray-Curtis',
    'unweighted_unifrac': 'Unweighted UniFrac',
}


def _real_column_names(metadata_file):
    """{sanitized file-name component: real column name}. The pipeline
    names beta-group-significance-<column>-<metric>.qzv and
    sample-classifier-<column>/ after metadata_utils.
    safe_filename_component(column), so a column like "Host Plant" comes
    back from those paths as "Host_Plant" -- which matches nothing in the
    metadata: every sample then lost its group on that column's PCoA/
    dendrogram pages, and the tables showed the mangled name."""
    return {metadata_utils.safe_filename_component(column): column
            for column in metadata_utils.parse_metadata_columns(metadata_file)}


def _split_beta_stem(stem):
    stem = stem.removeprefix('beta-group-significance-')
    for metric in _BETA_DISTANCE_METRICS:
        suffix = f'-{metric}'
        if stem.endswith(suffix):
            return stem[:-len(suffix)], metric
    return stem, ''


_SIGNIFICANCE_THRESHOLD = 0.05


def _significant_columns(pvalues_by_column):
    """Columns with at least one p-value below _SIGNIFICANCE_THRESHOLD.
    `pvalues_by_column` is {column: [p_value, ...]}; entries may be None
    (a test that couldn't be computed for that column/metric)."""
    return sorted(column for column, pvalues in pvalues_by_column.items()
                  if any(p is not None and p < _SIGNIFICANCE_THRESHOLD for p in pvalues))


def _significance_note(column, significant_columns):
    """One extra sentence on why this particular column's figure appears --
    once more than one column can drive a figure, an undifferentiated run of
    same-shaped pages reads as a dump rather than a report; this makes each
    page's reason for being there explicit."""
    if column in significant_columns:
        return (f' Shown here because {column} came back statistically significant '
                f'(p < {_SIGNIFICANCE_THRESHOLD}) above.')
    return (f' Shown here as the report\'s default grouping column, though it did not reach '
            f'statistical significance for this comparison.')


# Already shown as their own rows on the "Run information" page -- excluded
# from "Pipeline parameters" so the two pages don't just repeat each other.
_PARAMETERS_SHOWN_ELSEWHERE = {'input', 'output', 'metadata', 'classifier', 'qiime2'}

# One short explanatory paragraph per analysis section, printed under its
# title -- what the analysis actually shows and how to read it, for a reader
# who isn't already familiar with these specific QIIME2 outputs.
_INTRO_DADA2 = (
    'DADA2 denoises raw reads into exact amplicon sequence variants (ASVs) and removes chimeric '
    'artifacts. This table shows what fraction of each sample\'s reads survived quality filtering, '
    'denoising, merging (paired-end only), and chimera removal -- a sample retaining very few reads '
    'at any step may need different DADA2 parameters, or exclusion from downstream analysis.'
)
_INTRO_ALPHA_TABLE = (
    'Alpha diversity measures within-sample richness and evenness of taxa (Faith\'s phylogenetic '
    'diversity, observed features, Shannon, and evenness). A Kruskal-Wallis test checks whether a '
    'metric differs significantly between the groups of a metadata column; a low p-value '
    '(conventionally < 0.05) indicates the groups differ.'
)
_INTRO_ALPHA_BOXPLOT = (
    'Boxplots of each alpha diversity metric, split by this metadata column -- the same comparison '
    'as the group-significance table above, shown visually.'
)
_INTRO_BETA_TABLE = (
    'Beta diversity compares community composition (not just richness) between samples. PERMANOVA '
    'tests whether samples within the same metadata group are more similar to each other than to '
    'samples in other groups, for both Bray-Curtis (abundance-weighted) and unweighted UniFrac '
    '(phylogenetic presence/absence) distances.'
)
_INTRO_PCOA = (
    'Principal coordinates analysis (PCoA) projects the pairwise distance matrix into two dimensions '
    'for visualization; samples that cluster together have more similar community composition. The '
    'percentage on each axis is the proportion of total variance it explains.'
)
_INTRO_DENDROGRAM = (
    'UPGMA hierarchical clustering of samples from the same pairwise distance matrix used for the '
    'PCoA above; samples that join at a shorter distance are more similar in community composition. '
    'Unlike PCoA, this preserves the exact pairwise distances rather than a two-dimensional '
    'approximation of them.'
)
_INTRO_GENUS = (
    'Relative abundance of the most abundant genera in each sample, collapsed from the classifier\'s '
    'taxonomic assignments. "Unclassified" covers features the classifier could not assign to a '
    'genus; "Other" covers genera outside the most abundant ones shown individually.'
)
_INTRO_RAREFACTION = (
    'Number of observed features (ASVs) as a function of sequencing depth, per sample. A curve that '
    'plateaus indicates sequencing depth was sufficient to capture most of that sample\'s diversity; '
    'one still rising steeply at the sampling depth used suggests deeper sequencing would likely '
    'reveal more.'
)
_INTRO_CLASSIFIER = (
    'Tests whether community composition alone can predict a sample\'s metadata category, using a '
    'random-forest classifier. "Baseline accuracy" is what a naive classifier predicting only the '
    'most common class would achieve; an "accuracy ratio" above 1 means community composition '
    'carries real predictive signal for that column.'
)
_INTRO_SEQ_LENGTH = (
    'Length distribution of the representative sequences (one per ASV) produced by DADA2. ITS regions '
    'are naturally variable in length, but a distribution that is unexpectedly narrow, wide, or '
    'off-target suggests the primer/extraction step included flanking regions it should not have, or '
    'trimmed too aggressively.'
)
_INTRO_CONFIDENCE = (
    'Classifier confidence score for each ASV\'s taxonomic assignment (1.0 is fully confident). A '
    'distribution skewed toward low values suggests many ASVs are only weakly resolved by the '
    'reference database used to train the classifier.'
)


def _run_info_rows(run_metadata):
    pipeline = run_metadata.get('pipeline', {})
    env = run_metadata.get('environment', {})
    inputs = run_metadata.get('inputs', {})
    duration = pipeline.get('duration_seconds')
    return [
        ('Run started', pipeline.get('start_time', '')),
        ('Run finished', pipeline.get('end_time', '')),
        ('Run duration', timing.format_elapsed(duration) if duration is not None else ''),
        # Both can be blanked in run_metadata.json before sharing a report
        # (rebuild it with report.build_report()); a bare "@" reads as a bug.
        ('Run by', f"{env.get('username', '')}@{env.get('hostname', '')}"
         if env.get('username') or env.get('hostname') else '(not recorded)'),
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


def _alpha_boxplot_figure(alpha_results, report_column, group_color=None):
    """`group_color` ({group: color}, see _group_color_map) keeps each
    group's box the same color as its PCoA dots/dendrogram labels; boxes
    used to be colored by position among the groups present, which shifts
    every color as soon as one group is absent from the alpha results."""
    group_color = group_color or {}
    metrics = sorted(alpha_results)
    n = len(metrics)
    # A near-square grid (2x2 for the usual 4 metrics) instead of a single
    # wide row -- each subplot gets roughly 4x the area, since it now scales
    # with both page dimensions instead of just getting squeezed narrower.
    n_cols = math.ceil(math.sqrt(n)) if n else 1
    n_rows = math.ceil(n / n_cols) if n else 1
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 4.5 * n_rows), squeeze=False)
    flat_axes = axes.flatten()
    n_groups = max((len(alpha_results[m].get(report_column, {}).get('groups', {})) for m in metrics), default=0)
    # 45-degree labels' horizontal footprint grows with label length -- fine
    # for a handful of short categories, but a column with many (and often
    # long, e.g. env-local-scale's "self-heating coal mine waste material")
    # category names crowds adjacent labels into each other. Fully vertical
    # labels take the same small, length-independent horizontal footprint
    # regardless of how long the label text is.
    rotation = 90 if n_groups > _MANY_SAMPLES_THRESHOLD else 45
    for ax, metric in zip(flat_axes, metrics):
        groups = alpha_results[metric].get(report_column, {}).get('groups', {})
        labels = sorted(groups)
        data = [groups[label] for label in labels]
        if data:
            bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, medianprops={'color': 'black'})
            for i, patch in enumerate(bp['boxes']):
                patch.set_facecolor(group_color.get(labels[i], _COLORBLIND_PALETTE[i % len(_COLORBLIND_PALETTE)]))
                patch.set_alpha(0.75)
        ax.set_title(metric)
        ax.tick_params(axis='x', rotation=rotation)
    for ax in flat_axes[n:]:
        ax.set_visible(False)
    fig.suptitle(f'Alpha diversity by {report_column}', fontsize=13)
    fig.tight_layout()
    return fig


# Label for a sample with no value in the grouping column (or no metadata
# row at all). Never in _group_color_map/_group_marker_map, so it always
# takes their callers' neutral fallback (grey, plain circle) instead of
# using up a palette color. An empty-string label would also be silently
# dropped from the legend by matplotlib.
_MISSING_GROUP = '(missing)'


def _sample_group(metadata_table, sample_id, column):
    return metadata_table.get(sample_id, {}).get(column) or _MISSING_GROUP


def _column_groups(metadata_table, column):
    return sorted({row.get(column) for row in metadata_table.values()} - {None, ''})


def _group_color_map(metadata_table, report_column):
    """{group: color}, built once from every value of `report_column`
    across the *whole* metadata table -- not just one figure's own sample
    subset -- so the same group gets the same color on every figure that
    colors by this column (PCoA, dendrogram), whatever subset of samples
    that particular figure happens to plot."""
    if not report_column:
        return {}
    groups = _column_groups(metadata_table, report_column)
    return {group: _COLORBLIND_PALETTE[i % len(_COLORBLIND_PALETTE)] for i, group in enumerate(groups)}


# Cycled once per full pass through _COLORBLIND_PALETTE (see _group_marker_map)
# -- 11 colors x 8 markers covers up to 88 groups before a group repeats both
# its color and its marker, comfortably past any realistic metadata column's
# cardinality.
_MARKER_CYCLE = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']


def _group_marker_map(metadata_table, report_column):
    """{group: marker}, using the same group order as _group_color_map so
    the two stay aligned. A column with more groups than _COLORBLIND_PALETTE
    has colors makes two groups share a color; cycling the marker shape
    once per full pass through the palette means any two groups sharing a
    color always differ in shape instead. Only meaningful for a scatter
    plot (PCoA) -- the dendrogram encodes group via colored text, which has
    no marker-shape equivalent."""
    if not report_column:
        return {}
    groups = _column_groups(metadata_table, report_column)
    return {group: _MARKER_CYCLE[(i // len(_COLORBLIND_PALETTE)) % len(_MARKER_CYCLE)]
            for i, group in enumerate(groups)}


def _wrap_legend_label(label, width=18):
    """Wrap a long legend label onto multiple lines. A single very long
    metadata value (e.g. "self-heating coal mine waste material") used as
    one legend entry otherwise forces the whole legend column wide enough
    to squeeze the actual plot into a fraction of the page -- every other
    entry pays for the longest one's width whether it needs to or not."""
    return '\n'.join(textwrap.wrap(str(label), width=width)) or str(label)


def _readable_text_color(hex_color):
    """Darken `hex_color` if it's too pale to read as small text on a white
    page (the palette's sand, most visibly) -- used only where a palette
    color labels text directly (dendrogram leaf labels/legend), never for
    the dots/bars/lines elsewhere that same color fills, which read fine at
    full brightness."""
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (1, 3, 5))
    brightness = (299 * r + 587 * g + 114 * b) / 1000
    if brightness <= 190:
        return hex_color
    scale = 190 / brightness
    return '#{:02x}{:02x}{:02x}'.format(round(r * scale), round(g * scale), round(b * scale))


def _pcoa_figure(sample_coords, proportion_explained, metadata_table, report_column, title,
                  group_color=None, group_marker=None):
    fig, ax = plt.subplots(figsize=(6, 5))
    groups = sorted({_sample_group(metadata_table, sid, report_column) for sid in sample_coords})
    group_color = group_color if group_color is not None else _group_color_map(metadata_table, report_column)
    group_marker = group_marker if group_marker is not None else _group_marker_map(metadata_table, report_column)
    for group in groups:
        xs, ys = [], []
        for sid, (x, y) in sample_coords.items():
            if _sample_group(metadata_table, sid, report_column) == group:
                xs.append(x)
                ys.append(y)
        color = group_color.get(group, '#666666')
        # Semi-transparent fill + a thin outline: overlapping points (common
        # with pilot-scale sample counts) stay distinguishable from each
        # other and from a solid single-point marker instead of merging
        # into one opaque blob. Marker shape (not just color) also encodes
        # group past _COLORBLIND_PALETTE's length, where two groups
        # otherwise share a color.
        ax.scatter(xs, ys, label=_wrap_legend_label(group), color=color, marker=group_marker.get(group, 'o'),
                   alpha=0.7, s=55, edgecolors='black', linewidths=0.6)
    ax.set_xlabel(f'PC1 ({proportion_explained[0] * 100:.1f}%)')
    ax.set_ylabel(f'PC2 ({proportion_explained[1] * 100:.1f}%)')
    ax.set_title(title)
    # Placed outside the axes (rather than loc='best' inside it) so the
    # legend never sits on top of a data point, whatever the point cloud's
    # shape happens to be for a given run's data.
    ax.legend(fontsize=9, frameon=True, bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()
    return fig


def _dendrogram_figure(distance_df, metadata_table, report_column, title, group_color=None):
    """UPGMA (average-linkage) hierarchical clustering of samples from a
    square distance matrix -- the standard complement to a PCoA scatter for
    the same distance metric. Leaf labels (not link colors) are colored by
    `report_column`: scipy's own automatic link coloring is a cluster-shape
    heuristic unrelated to metadata groups, and would visually compete with
    a metadata-driven color coding here."""
    sample_ids = list(distance_df.index)
    condensed = squareform(distance_df.values, checks=False)
    link = linkage(condensed, method='average')

    many_samples = len(sample_ids) > _MANY_SAMPLES_THRESHOLD
    if many_samples:
        # Leaves along the y-axis (one row per sample), like the genus
        # barplot's many-sample layout, so labels stay readable instead of
        # crowding along a fixed-width x-axis.
        fig, ax = plt.subplots(figsize=(7, max(5, 0.25 * len(sample_ids))))
        orientation = 'left'
    else:
        fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(sample_ids)), 5))
        orientation = 'top'

    dendrogram(link, labels=sample_ids, ax=ax, orientation=orientation,
               leaf_rotation=0 if many_samples else 90, leaf_font_size=9,
               color_threshold=0, above_threshold_color='#444444')

    if report_column:
        groups = sorted({_sample_group(metadata_table, sid, report_column) for sid in sample_ids})
        group_color = group_color if group_color is not None else _group_color_map(metadata_table, report_column)
        # Darkened for text/legend readability (some palette colors, e.g.
        # yellow, are too pale to read as small text even though they're
        # fine as a dot/bar/line fill) -- the underlying group_color mapping
        # itself stays untouched so it matches the PCoA page's dot colors.
        # Same '#666666' fallback _pcoa_figure uses for a group missing from
        # group_color (a sample in the distance matrix but not metadata_table)
        # -- group_color_map's own docstring promises the same group gets the
        # same color on every figure; two different hardcoded fallback colors
        # here and there would quietly break that for this edge case.
        text_color = {group: _readable_text_color(group_color.get(group, '#666666')) for group in groups}
        tick_labels = ax.get_yticklabels() if orientation == 'left' else ax.get_xticklabels()
        for tick_label in tick_labels:
            group = _sample_group(metadata_table, tick_label.get_text(), report_column)
            tick_label.set_color(text_color.get(group, 'black'))
        handles = [Line2D([0], [0], color=text_color[group], lw=4, label=_wrap_legend_label(group))
                   for group in groups]
        if orientation == 'left':
            # The leaf labels themselves (not just the plotted lines) sit at
            # the right edge of the axes here, so a legend placed just
            # outside that edge (as used for every other figure in this
            # report) would overlap the top labels' text instead of clearing
            # it. Placed above the axes instead, where nothing else is
            # drawn. Left-aligned rather than mode='expand': stretching a
            # short handful of entries across the full axes width read as
            # an odd gap-filled bar rather than a normal legend.
            ax.legend(handles=handles, fontsize=9, frameon=True, loc='lower left',
                      bbox_to_anchor=(0, 1.01), ncol=min(len(handles), 4))
        else:
            ax.legend(handles=handles, fontsize=9, frameon=True, bbox_to_anchor=(1.02, 1), loc='upper left')

    # A figure-level title (rather than ax.set_title): tight_layout reserves
    # room for it above everything else, including the legend placed just
    # above the axes in the 'left'-orientation case -- an axes-level title
    # sits right at the axes edge instead, which the legend there would
    # overlap.
    fig.suptitle(title, fontsize=13)
    if orientation == 'left':
        ax.set_xlabel('Distance')
    else:
        ax.set_ylabel('Distance')
    fig.tight_layout()
    return fig


_MANY_SAMPLES_THRESHOLD = 10


def _fit_width_mm(fig, max_height_mm=230, max_width_mm=180):
    """The fpdf placement width (mm) that keeps `fig` at most `max_height_mm`
    tall once placed on the page. add_figure_page's own default (width=180)
    assumes a roughly landscape-shaped figure; a much taller one needs a
    correspondingly narrower placement width or it runs off the page --
    fpdf places images by width alone and scales height to match the
    figure's own aspect ratio."""
    fig_w_in, fig_h_in = fig.get_size_inches()
    return min(max_width_mm, max_height_mm * fig_w_in / fig_h_in)


def _genus_colors(index):
    """One color per genus, cycling the colorblind palette -- except 'Other'
    and 'Unclassified', which always get their own fixed neutral greys
    (different from each other) rather than whatever color they'd land on
    next in the cycle."""
    colors = []
    next_color = 0
    for label in index:
        if label == 'Unclassified':
            colors.append(_NEUTRAL_GREY_UNCLASSIFIED)
        elif label == 'Other':
            colors.append(_NEUTRAL_GREY_OTHER)
        else:
            colors.append(_COLORBLIND_PALETTE[next_color % len(_COLORBLIND_PALETTE)])
            next_color += 1
    return colors


def _genus_barplot_figure(genus_table):
    n_samples = len(genus_table.columns)
    colors = _genus_colors(genus_table.index)
    if n_samples > _MANY_SAMPLES_THRESHOLD:
        # Horizontal bars (one row per sample) scale by adding figure
        # height, not width -- vertical bars for this many samples squeeze
        # every bar, and its sample label, down to an unreadable sliver. A
        # wider figure (11in, up from 8in) also matters here: add_figure_page
        # sizes its placement width off this figure's own aspect ratio
        # (_fit_width_mm), and the previous 8:height ratio left the bars
        # using well under the available page width once a real (>10-sample)
        # run made the figure this tall.
        fig, ax = plt.subplots(figsize=(11, min(24, max(6, 0.3 * n_samples))))
        genus_table.T.plot(kind='barh', stacked=True, ax=ax, legend=True, width=0.85, color=colors,
                           edgecolor='white', linewidth=0.4)
        ax.set_xlabel('Relative abundance', fontsize=12)
        ax.invert_yaxis()  # first sample at the top, matching reading order
    else:
        fig, ax = plt.subplots(figsize=(max(6, 0.6 * n_samples), 5))
        genus_table.T.plot(kind='bar', stacked=True, ax=ax, legend=True, width=0.85, color=colors,
                           edgecolor='white', linewidth=0.4)
        ax.set_ylabel('Relative abundance', fontsize=12)
    ax.tick_params(axis='both', labelsize=11)
    ax.legend(fontsize=12, bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True)
    fig.tight_layout()
    return fig


def _rarefaction_figure(curves, metric):
    # A one-line-per-sample legend easily runs to dozens of entries, which
    # dwarfed the actual plot (squeezed into a small corner) when this used
    # a single wide-and-short figure. A taller figure, plus a 2-column
    # legend past the sample-count threshold above, fixes both at once --
    # but not too tall: (7, 11) read as excessively elongated in practice,
    # so this is scaled back to roughly two-thirds of that height.
    many_samples = len(curves) > _MANY_SAMPLES_THRESHOLD
    fig, ax = plt.subplots(figsize=(7, 7) if many_samples else (6, 5))
    for i, (sample_id, points) in enumerate(curves.items()):
        if points:
            color = _COLORBLIND_PALETTE[i % len(_COLORBLIND_PALETTE)]
            ax.plot([p[0] for p in points], [p[1] for p in points], marker='o', markersize=3,
                    color=color, label=sample_id)
    ax.set_xlabel('Sequencing depth')
    ax.set_ylabel(metric)
    ax.set_title(f'Rarefaction curve ({metric})')
    ax.legend(fontsize=8.5 if many_samples else 9, ncol=2 if many_samples else 1,
              bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True)
    fig.tight_layout()
    return fig


def _seq_length_histogram_figure(lengths):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(lengths, bins='auto', color=_COLORBLIND_PALETTE[0], edgecolor='white', linewidth=0.4)
    ax.set_xlabel('Sequence length (bp)')
    ax.set_ylabel('Number of ASVs')
    ax.set_title('Representative sequence length distribution')
    fig.tight_layout()
    return fig


def _confidence_histogram_figure(confidences):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(confidences, bins='auto', range=(0, 1), color=_COLORBLIND_PALETTE[1], edgecolor='white',
             linewidth=0.4)
    ax.set_xlabel('Classification confidence')
    ax.set_ylabel('Number of ASVs')
    ax.set_title('Taxonomic classification confidence distribution')
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
    final_sample_ids = metadata_utils.final_sample_ids(sample_frequencies, metadata_table)

    # Same eligibility rule cli/pipeline.py uses to decide which columns
    # beta-group-significance/alpha-group-significance actually get run
    # against -- computed once, used both to validate an explicit
    # report_column and to auto-pick one below.
    eligible_columns = metadata_utils.eligible_categorical_columns(metadata_file, final_sample_ids)

    if report_column is not None and report_column not in eligible_columns:
        # A typo'd/nonexistent --report-metadata-column used to fail
        # silently (every group-by-report_column lookup downstream just
        # returns nothing for a column that doesn't exist). A real but
        # ineligible column (numeric, or too few/too-small groups) was
        # worse: PCoA/dendrogram figures were still built for it -- the
        # metric-level artifacts they read exist independent of any
        # particular column -- captioned as "did not reach statistical
        # significance", which is wrong: no significance test was ever run
        # for a column eligibility itself excludes from testing.
        print(f'Warning: --report-metadata-column "{report_column}" is not an eligible categorical column '
              f'(>=2 distinct values, >=2 samples each) in {metadata_file} -- falling back to auto-selecting '
              f'an eligible column instead.')
        report_column = None

    if report_column is None:
        report_column = eligible_columns[0] if eligible_columns else None

    real_column = _real_column_names(metadata_file)

    pdf = _ReportPDF()
    pdf.set_auto_page_break(auto=True, margin=15)

    # 1. Title / summary
    summary_lines = [f'Samples: {len(final_sample_ids)}',
                      f'Default report grouping column: {report_column or "(none eligible)"}',
                      '(alpha/beta diversity figures below also cover any other column with a '
                      'statistically significant result)']
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
        pdf.add_table_page('DADA2 read retention', header, rows, intro=_INTRO_DADA2)

    # 2b. Representative sequence length distribution
    rep_seqs_fasta = output_folder / 'rep_seqs_export' / 'dna-sequences.fasta'
    if rep_seqs_fasta.exists():
        lengths = report_data.parse_fasta_sequence_lengths(rep_seqs_fasta)
        if lengths:
            seq_length_fig = _seq_length_histogram_figure(lengths)
            pdf.add_figure_page('Sequence length distribution', seq_length_fig, intro=_INTRO_SEQ_LENGTH)

    # 3. Alpha diversity: group-significance p-values + boxplots
    alpha_results = {}
    for qzv_path in sorted(output_folder.glob('alpha-group-significance-*.qzv')):
        metric = qzv_path.stem.removeprefix('alpha-group-significance-')
        alpha_results[metric] = report_data.parse_alpha_group_significance(qzv_path)

    if alpha_results:
        rows = []
        alpha_pvalues_by_column = {}
        for metric, by_column in alpha_results.items():
            for column, stats in by_column.items():
                rows.append([metric, column, f'{stats["h_statistic"]:.3f}', f'{stats["p_value"]:.3f}'])
                alpha_pvalues_by_column.setdefault(column, []).append(stats.get('p_value'))
        pdf.add_table_page('Alpha diversity group significance (Kruskal-Wallis)',
                            ['metric', 'metadata column', 'H', 'p-value'], rows, intro=_INTRO_ALPHA_TABLE)

        # One boxplot page per metadata column that's either the report's
        # default grouping column or came back significant above -- not
        # every eligible column, which would bury the columns that actually
        # show something behind however many don't.
        alpha_significant = _significant_columns(alpha_pvalues_by_column)
        alpha_columns = sorted(set(alpha_significant) | ({report_column} if report_column else set()))
        for column in alpha_columns:
            if any(column in by_column for by_column in alpha_results.values()):
                alpha_fig = _alpha_boxplot_figure(alpha_results, column,
                                                   group_color=_group_color_map(metadata_table, column))
                intro = _INTRO_ALPHA_BOXPLOT + _significance_note(column, alpha_significant)
                pdf.add_figure_page(f'Alpha diversity by {column}', alpha_fig, width=_fit_width_mm(alpha_fig),
                                     intro=intro)

    # 4. Beta diversity: PCoA + PERMANOVA
    beta_rows = []
    beta_pvalues_by_column = {}
    for qzv_path in sorted(output_folder.glob('beta-group-significance-*.qzv')):
        stats = report_data.parse_beta_group_significance(qzv_path)
        column, metric = _split_beta_stem(qzv_path.stem)
        column = real_column.get(column, column)
        beta_rows.append([column, metric, stats.get('test_statistic_name'),
                           f'{stats["test_statistic"]:.3f}' if stats.get('test_statistic') is not None else '',
                           f'{stats["p_value"]:.3f}' if stats.get('p_value') is not None else ''])
        beta_pvalues_by_column.setdefault(column, []).append(stats.get('p_value'))
    if beta_rows:
        pdf.add_table_page('Beta diversity group significance (PERMANOVA)',
                            ['metadata column', 'distance metric', 'statistic', 'value', 'p-value'], beta_rows,
                            intro=_INTRO_BETA_TABLE)

    # One sub-section per metadata column that's either the report's default
    # grouping column or came back significant above (same selection as the
    # alpha boxplots) -- grouped by column rather than interleaving columns
    # and metrics, so everything about one metadata factor (both metrics,
    # both PCoA and dendrogram) reads together instead of scattering related
    # pages apart.
    beta_significant = _significant_columns(beta_pvalues_by_column)
    beta_columns = sorted(set(beta_significant) | ({report_column} if report_column else set()))
    for column in beta_columns:
        # Computed once per column so the same group gets the same color
        # (and, for PCoA, the same marker shape) on both the PCoA and
        # dendrogram pages, for both distance metrics, rather than each
        # figure picking its own colors from whatever subset of
        # samples/groups it happens to see.
        group_color = _group_color_map(metadata_table, column)
        group_marker = _group_marker_map(metadata_table, column)
        note = _significance_note(column, beta_significant)
        # Only worth telling the reader about when it actually kicks in --
        # with few enough groups that every one gets its own palette color,
        # every marker is a plain circle and the sentence would just be
        # confusing noise.
        marker_note = (' Marker shape also distinguishes groups that share a color (this column has '
                        'more groups than the palette has colors).' if len(set(group_marker.values())) > 1 else '')
        for metric in ('bray_curtis', 'unweighted_unifrac'):
            metric_display = _METRIC_DISPLAY_NAMES.get(metric, metric)
            ordination_path = output_folder / 'core-metrics-results' / f'{metric}_pcoa_export' / 'ordination.txt'
            if ordination_path.exists():
                sample_coords, proportion_explained = report_data.parse_ordination(ordination_path)
                pdf.add_figure_page(f'{column}: {metric_display} PCoA', _pcoa_figure(
                    sample_coords, proportion_explained, metadata_table, column,
                    f'{metric_display} (colored by {column})', group_color=group_color,
                    group_marker=group_marker), intro=_INTRO_PCOA + marker_note + note)

            distance_path = output_folder / 'core-metrics-results' / f'{metric}_distance_export' / \
                'distance-matrix.tsv'
            if distance_path.exists():
                distance_df = report_data.parse_distance_matrix(distance_path)
                dendrogram_fig = _dendrogram_figure(distance_df, metadata_table, column,
                                                     f'{metric_display} (UPGMA, colored by {column})',
                                                     group_color=group_color)
                pdf.add_figure_page(f'{column}: {metric_display} sample clustering', dendrogram_fig,
                                     width=_fit_width_mm(dendrogram_fig), intro=_INTRO_DENDROGRAM + note)

    # 5. Genus-level composition
    biom_taxo_path = output_folder / 'biom_table' / 'table-with-taxonomy.biom.tsv'
    if biom_taxo_path.exists():
        genus_table = report_data.build_genus_abundance_table(biom_taxo_path)
        genus_fig = _genus_barplot_figure(genus_table)
        pdf.add_figure_page('Genus-level relative abundance', genus_fig, width=_fit_width_mm(genus_fig),
                             intro=_INTRO_GENUS)

    # 5b. Taxonomic classification confidence distribution
    taxonomy_tsv_path = output_folder / 'biom_table' / 'taxonomy.tsv'
    if taxonomy_tsv_path.exists():
        confidences = report_data.parse_taxonomy_confidence(taxonomy_tsv_path)
        if confidences:
            confidence_fig = _confidence_histogram_figure(confidences)
            pdf.add_figure_page('Classification confidence', confidence_fig, intro=_INTRO_CONFIDENCE)

    # 6. Rarefaction curve
    rarefaction_qzv = output_folder / 'alpha-rarefaction.qzv'
    if rarefaction_qzv.exists():
        curves = report_data.parse_rarefaction_curve(rarefaction_qzv, 'observed_features')
        rarefaction_fig = _rarefaction_figure(curves, 'observed_features')
        pdf.add_figure_page('Rarefaction curve', rarefaction_fig, width=_fit_width_mm(rarefaction_fig),
                             intro=_INTRO_RAREFACTION)

    # 7. Sample classifier accuracy, if any
    classifier_rows = []
    for accuracy_qzv in sorted(output_folder.glob('sample-classifier-*/accuracy_results.qzv')):
        column = accuracy_qzv.parent.name.removeprefix('sample-classifier-')
        column = real_column.get(column, column)
        accuracy = report_data.parse_classifier_accuracy(accuracy_qzv)
        for key, label in (('overall_accuracy', 'Overall accuracy'),
                            ('baseline_accuracy', 'Baseline accuracy'),
                            ('accuracy_ratio', 'Accuracy ratio')):
            if key in accuracy:
                classifier_rows.append([column, label, f'{accuracy[key]:.3f}'])
    if classifier_rows:
        pdf.add_table_page('Sample classifier results', ['metadata column', 'metric', 'value'],
                            classifier_rows, intro=_INTRO_CLASSIFIER)

    report_path = output_folder / 'report.pdf'
    pdf.output(str(report_path))
    return report_path
