"""Pure data-extraction helpers for the PDF summary report.

Each function parses one already-on-disk pipeline artifact (a plain-text
export, or a .qzv read directly as a zip) into plain Python/pandas
structures. None of them touch matplotlib, fpdf2, or subprocess, so they're
fully unit-testable with small fixture files -- report.py is the thin layer
that turns this data into pages.

Internal .qzv layouts (jsonp variable names, the beta-group-significance
"Overview" HTML table, the per-metric rarefaction CSVs) were confirmed by
generating real files against real data and inspecting them, not guessed
from documentation.
"""
import csv
import json
import re
import zipfile

import pandas as pd


def parse_sample_frequencies(sample_frequencies_tsv_path):
    """Parse the plain-text export of `sample-frequencies.qza` (produced by
    `qiime feature-table summarize`) into {sample_id: total_reads}.

    Needed because a near-empty input sample can survive DADA2 as a
    zero-read row in the final feature table (retain-all-samples defaults to
    True) -- QIIME2's own beta-group-significance/sample-classifier actions
    then treat that as a group of size zero for that sample's class, which
    can silently turn a nominally-eligible metadata column into an
    ineligible one. Callers should exclude zero-frequency samples before
    computing group eligibility.

    QIIME2 writes this file with thousands separators once a sample's
    frequency reaches four digits (e.g. "110,406.0", confirmed against a
    real 53-sample production run whose samples had ~100k+ reads each --
    every fixture used before that had frequencies too small to trigger it).
    """
    frequencies = {}
    with open(sample_frequencies_tsv_path) as f:
        f.readline()  # header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if not fields or fields[0] in ('', '#q2:types'):
                continue
            frequencies[fields[0]] = float(fields[1].replace(',', ''))
    return frequencies


def parse_dada2_stats(stats_tsv_path):
    """DADA2 denoising-stats.tsv export as a DataFrame indexed by sample-id,
    numeric columns coerced to numbers."""
    df = pd.read_csv(stats_tsv_path, sep='\t')
    df = df[df['sample-id'] != '#q2:types']
    df = df.set_index('sample-id')
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


def parse_ordination(ordination_txt_path):
    """Parse a skbio OrdinationResults text file (as exported from a QIIME2
    PCoAResults artifact). Returns (sample_coords, proportion_explained):
    sample_coords is {sample_id: (pc1, pc2)}, proportion_explained is
    (pc1_fraction, pc2_fraction).
    """
    with open(ordination_txt_path) as f:
        lines = [line.rstrip('\n') for line in f]

    proportion_explained = (0.0, 0.0)
    sample_coords = {}

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('Proportion explained'):
            values = [float(v) for v in lines[i + 1].split('\t') if v]
            proportion_explained = (values[0], values[1] if len(values) > 1 else 0.0)
            i += 2
        elif line.startswith('Site\t'):
            n_samples = int(line.split('\t')[1])
            for row in lines[i + 1:i + 1 + n_samples]:
                fields = row.split('\t')
                coords = [float(v) for v in fields[1:]]
                sample_coords[fields[0]] = (coords[0], coords[1] if len(coords) > 1 else 0.0)
            i += 1 + n_samples
        else:
            i += 1
    return sample_coords, proportion_explained


_ALPHA_STAT_RE = re.compile(r'\{"H":\s*([0-9.eE+-]+),\s*"p":\s*([0-9.eE+-]+)\}')


def _extract_json_object(text, start):
    """Extract one balanced `{...}` JSON object from `text`, `start` pointing
    at the opening brace. A plain regex can't do this since the object
    contains nested braces/brackets."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return json.loads(text[start:i + 1])
    raise ValueError(f'Unbalanced JSON object starting at index {start}')


def parse_alpha_group_significance(qzv_path):
    """Parse an alpha-group-significance .qzv. Returns
    {column_name: {'h_statistic': float, 'p_value': float, 'groups': {group_name: [values]}}}
    for every categorical column QIIME2 tested. `groups` is the same raw
    per-group alpha-diversity values QIIME2 itself plots, reused here instead
    of needing a separate alpha-vector export just to draw a boxplot.
    """
    results = {}
    with zipfile.ZipFile(qzv_path) as zf:
        for name in zf.namelist():
            match = re.search(r'/data/column-(.+)\.jsonp$', name)
            if not match:
                continue
            column = match.group(1)
            content = zf.read(name).decode('utf-8')

            stat_match = _ALPHA_STAT_RE.search(content)
            if not stat_match:
                continue

            groups = {}
            args_start = content.index('(') + 1
            first_comma = content.index(',', args_start)
            group_data_start = content.find('{', first_comma)
            if group_data_start != -1:
                group_data = _extract_json_object(content, group_data_start)
                for label, values in zip(group_data.get('index') or [], group_data.get('data') or []):
                    group_name = label.rsplit(' (n=', 1)[0] if ' (n=' in label else label
                    groups[group_name] = values

            results[column] = {
                'h_statistic': float(stat_match.group(1)),
                'p_value': float(stat_match.group(2)),
                'groups': groups,
            }
    return results


_BETA_ROW_RE = re.compile(r'<th>([^<]+)</th>\s*<td>([^<]+)</td>')


def parse_beta_group_significance(qzv_path):
    """Parse a beta-group-significance .qzv's "Overview" table (method,
    test statistic, sample/group counts, p-value)."""
    with zipfile.ZipFile(qzv_path) as zf:
        index_name = next(n for n in zf.namelist() if n.endswith('/data/index.html'))
        html = zf.read(index_name).decode('utf-8')

    raw = dict(_BETA_ROW_RE.findall(html))

    def _num(key, cast):
        value = raw.get(key)
        return cast(value) if value not in (None, '') else None

    return {
        'method_name': raw.get('method name'),
        'test_statistic_name': raw.get('test statistic name'),
        'sample_size': _num('sample size', int),
        'number_of_groups': _num('number of groups', int),
        'test_statistic': _num('test statistic', float),
        'p_value': _num('p-value', float),
    }


def parse_classifier_accuracy(qzv_path):
    """Parse a sample-classifier `classify-samples` accuracy_results.qzv's
    plain-text `predictive_accuracy.tsv` (a confusion matrix followed by
    summary rows). Returns {'overall_accuracy', 'baseline_accuracy',
    'accuracy_ratio'} -> float, omitting any not present.

    Deliberately not the HTML confusion-matrix table on the same qzv's
    index.html: that table has a variable number of columns (one per class)
    with the accuracy value in whichever is last, which a fixed
    `<th>...</th><td>...</td>` pair extraction (as used for
    beta-group-significance's simpler two-column "Overview" table) would
    silently misread as an empty leading cell.
    """
    with zipfile.ZipFile(qzv_path) as zf:
        tsv_name = next(n for n in zf.namelist() if n.endswith('/data/predictive_accuracy.tsv'))
        content = zf.read(tsv_name).decode('utf-8')

    wanted = {
        'Overall Accuracy': 'overall_accuracy',
        'Baseline Accuracy': 'baseline_accuracy',
        'Accuracy Ratio': 'accuracy_ratio',
    }
    result = {}
    for line in content.splitlines():
        fields = line.split('\t')
        key = wanted.get(fields[0])
        if key is None:
            continue
        values = [f for f in fields[1:] if f.strip()]
        if values:
            result[key] = float(values[-1])
    return result


def parse_rarefaction_curve(qzv_path, metric):
    """Parse one metric's raw per-sample rarefaction CSV from an
    alpha-rarefaction .qzv. Returns {sample_id: [(depth, mean_value), ...]},
    points sorted by depth, averaged across iterations at each depth.
    """
    with zipfile.ZipFile(qzv_path) as zf:
        csv_name = next(n for n in zf.namelist() if n.endswith(f'/data/{metric}.csv'))
        content = zf.read(csv_name).decode('utf-8').splitlines()

    reader = csv.DictReader(content)
    depth_cols = [c for c in reader.fieldnames if c.startswith('depth-')]
    depths = sorted({int(c.split('_')[0].split('-')[1]) for c in depth_cols})

    curves = {}
    for row in reader:
        points = []
        for depth in depths:
            values = [float(row[c]) for c in depth_cols
                      if c.startswith(f'depth-{depth}_') and row[c]]
            if values:
                points.append((depth, sum(values) / len(values)))
        curves[row['sample-id']] = points
    return curves


def parse_biom_taxonomy_tsv(biom_tsv_path):
    """table-with-taxonomy.biom.tsv export as a DataFrame indexed by feature
    id, with one column per sample plus a 'taxonomy' column."""
    df = pd.read_csv(biom_tsv_path, sep='\t', skiprows=1)
    return df.rename(columns={df.columns[0]: 'feature_id'}).set_index('feature_id')


def _genus_from_taxonomy(taxonomy):
    for field in taxonomy.split(';'):
        field = field.strip()
        if field.startswith('g__') and field != 'g__' and 'unidentified' not in field:
            return field[3:]
    return 'Unclassified'


def build_genus_abundance_table(biom_tsv_path, top_n=10):
    """Collapse table-with-taxonomy.biom.tsv to genus-level relative
    abundance. Returns a DataFrame: rows = genus (top `top_n`, plus 'Other'
    if anything was dropped), columns = samples, values = that sample's
    fraction of reads belonging to that genus.
    """
    df = parse_biom_taxonomy_tsv(biom_tsv_path)
    sample_cols = [c for c in df.columns if c != 'taxonomy']

    df = df.copy()
    df['genus'] = df['taxonomy'].map(_genus_from_taxonomy)
    genus_table = df.groupby('genus')[sample_cols].sum()

    totals = genus_table.sum(axis=0)
    relative = genus_table.div(totals.replace(0, 1), axis=1)

    top_genera = relative.sum(axis=1).sort_values(ascending=False).head(top_n).index
    other = relative.drop(index=top_genera).sum(axis=0)
    result = relative.loc[top_genera]
    if other.sum() > 0:
        result.loc['Other'] = other
    return result


def parse_fasta_sequence_lengths(fasta_path):
    """Per-record sequence lengths from a plain-text FASTA export (e.g.
    `qiime tools export`'s dna-sequences.fasta for rep-seqs.qza). Returns a
    list of ints, one per record; each record's sequence may itself be
    wrapped across several lines, so lines are accumulated between '>'
    headers rather than assumed to be one sequence per line."""
    lengths = []
    length = 0
    seen_header = False
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seen_header:
                    lengths.append(length)
                seen_header = True
                length = 0
            else:
                length += len(line)
    if seen_header:
        lengths.append(length)
    return lengths


def parse_taxonomy_confidence(taxonomy_tsv_path):
    """Per-feature classifier confidence scores from the exported, biom-
    header-rewritten taxonomy.tsv (biom_utils.rewrite_taxonomy_header's
    '#OTUID\\ttaxonomy\\tconfidence'). Returns a list of floats."""
    df = pd.read_csv(taxonomy_tsv_path, sep='\t')
    confidence_col = next(c for c in df.columns if c.strip().lower() == 'confidence')
    return pd.to_numeric(df[confidence_col], errors='coerce').dropna().tolist()


def parse_distance_matrix(distance_matrix_tsv_path):
    """Square distance matrix export (e.g. bray_curtis_distance_matrix.qza,
    via `qiime tools export`) as a DataFrame indexed and columned by sample
    id, in the same order as the file's own row/column order."""
    return pd.read_csv(distance_matrix_tsv_path, sep='\t', index_col=0)


def parse_run_metadata(run_metadata_json_path):
    """Parse run_metadata.json (written by cli/pipeline.py via
    provenance.write_run_metadata()) into a plain dict, or None if it
    doesn't exist -- e.g. an output folder from before this existed -- or
    can't be parsed as JSON -- e.g. a leftover truncated file from an older
    run that predates write_run_metadata()'s atomic write, or a process
    killed while the file itself was being read. Either way this is one of
    several optional QA/provenance pages; build_report() degrades to
    skipping them rather than the whole report crashing over one missing or
    unreadable file."""
    try:
        with open(run_metadata_json_path) as fh:
            return json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return None
