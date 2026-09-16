"""Parsing of QIIME2 metadata TSV files: column type detection and grouping
eligibility.

Pure Python, no QIIME2 dependency -- lets the pipeline decide which metadata
columns to run group-significance tests and sample-classifier against without
needing the qiime2 Python API (which would tie this module to running inside
an activated QIIME2 environment just to be imported).
"""
from collections import defaultdict


def _is_numeric(value):
    if value == '':
        return True  # a missing value doesn't disqualify numeric inference
    try:
        float(value)
        return True
    except ValueError:
        return False


def read_metadata_rows(path):
    """Return (header, types_row, data_rows) for a QIIME2 metadata TSV.

    `header` is the full header row (including the sample-id column).
    `types_row` is the `#q2:types` directive's per-column values, or None if
    the file doesn't have one. `data_rows` is every remaining row, each a
    list of fields (sample-id first).
    """
    with open(path) as f:
        lines = [line.rstrip('\n') for line in f if line.strip()]
    header = lines[0].split('\t')

    types_row = None
    data_rows = []
    for line in lines[1:]:
        fields = line.split('\t')
        if fields[0] == '#q2:types':
            types_row = fields[1:]
        else:
            data_rows.append(fields)
    return header, types_row, data_rows


def parse_metadata_columns(path):
    """Return {column_name: 'categorical'|'numeric'} for every non-id column.

    Honors an explicit `#q2:types` directive if present; otherwise infers a
    column as numeric only if every non-missing value in it parses as a
    float, matching QIIME2's own type-inference rule.
    """
    header, types_row, data_rows = read_metadata_rows(path)
    columns = header[1:]

    if types_row is not None:
        return {col: ('numeric' if t.strip().lower() == 'numeric' else 'categorical')
                for col, t in zip(columns, types_row)}

    result = {}
    for i, col in enumerate(columns):
        values = [row[i + 1] for row in data_rows if len(row) > i + 1]
        result[col] = 'numeric' if values and all(_is_numeric(v) for v in values) else 'categorical'
    return result


def read_metadata_table(path):
    """Return {sample_id: {column: value}}."""
    header, _types_row, data_rows = read_metadata_rows(path)
    columns = header[1:]
    return {row[0]: dict(zip(columns, row[1:])) for row in data_rows}


def class_sizes(path, column, sample_ids):
    """{value: count} for `column`, restricted to `sample_ids`."""
    table = read_metadata_table(path)
    counts = defaultdict(int)
    for sample_id in sample_ids:
        if sample_id in table:
            counts[table[sample_id].get(column, '')] += 1
    return dict(counts)


def eligible_categorical_columns(path, sample_ids, min_per_group=2):
    """Categorical columns with >=2 distinct values, each held by at least
    `min_per_group` of the given `sample_ids`. Preserves metadata-file column
    order.
    """
    column_types = parse_metadata_columns(path)
    eligible = []
    for col, col_type in column_types.items():
        if col_type != 'categorical':
            continue
        counts = class_sizes(path, col, sample_ids)
        groups = [n for n in counts.values() if n >= min_per_group]
        if len(groups) >= 2:
            eligible.append(col)
    return eligible


def has_alpha_group_significance_column(path, sample_ids):
    """True if at least one categorical column is usable for
    `qiime diversity alpha-group-significance`.

    That action tests every categorical column in the metadata file in one
    call and fails outright (for the whole file, not just skipping a column)
    if none qualify -- its own requirement is looser than
    `eligible_categorical_columns`'s >=2-per-group rule: a column just needs
    2+ distinct values with at least one repeated (not every value unique,
    which reads as an ID column) and not every sample sharing one value.
    """
    column_types = parse_metadata_columns(path)
    for col, col_type in column_types.items():
        if col_type != 'categorical':
            continue
        counts = class_sizes(path, col, sample_ids)
        if len(counts) >= 2 and any(n >= 2 for n in counts.values()):
            return True
    return False
