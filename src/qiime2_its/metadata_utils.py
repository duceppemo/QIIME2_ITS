"""Parsing of QIIME2 metadata TSV files: column type detection and grouping
eligibility.

Pure Python, no QIIME2 dependency -- lets the pipeline decide which metadata
columns to run group-significance tests and sample-classifier against without
needing the qiime2 Python API (which would tie this module to running inside
an activated QIIME2 environment just to be imported).
"""
import csv
import re
from collections import defaultdict

# Legacy ID-column headers QIIME2 accepts that start with "#" -- every other
# "#"-prefixed line is a comment (or a "#q2:" directive).
_HASH_ID_HEADERS = {'#SampleID', '#Sample ID', '#OTUID', '#OTU ID'}


def safe_filename_component(text):
    """Replace characters outside [A-Za-z0-9._-] with '_' for safe use as a
    filesystem path component. QIIME2 doesn't restrict metadata column names
    from containing '/' (or other characters with filesystem meaning), and
    the pipeline builds output paths directly from the column name -- a
    column like "site/plot" would otherwise be interpreted as a
    subdirectory, and QIIME2's own writers don't create missing parent
    directories for a single-artifact output, so it fails outright (a
    crash, not silently writing outside the output folder, but a fragile
    invariant to depend on regardless)."""
    return re.sub(r'[^A-Za-z0-9._-]', '_', text)


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
    # Same rules as QIIME2's own metadata reader: excel-tab CSV dialect
    # (quoted cells), leading/trailing whitespace stripped from every cell,
    # empty rows ignored, and "#"-prefixed rows are comments wherever they
    # appear -- except the legacy "#SampleID"-style ID headers and "#q2:"
    # directives. Comment rows used to be read as samples here.
    header = None
    types_row = None
    data_rows = []
    with open(path, newline='') as f:
        for fields in csv.reader(f, dialect='excel-tab'):
            fields = [field.strip() for field in fields]
            if not any(fields):
                continue
            first = fields[0]
            if header is None:
                if first.startswith('#') and first not in _HASH_ID_HEADERS:
                    continue
                header = fields
            elif first == '#q2:types':
                types_row = fields[1:]
            elif first.startswith('#'):
                continue
            else:
                data_rows.append(fields)
    if header is None:
        raise ValueError(f'{path} has no header row.')
    return header, types_row, data_rows


def _column_types(header, types_row, data_rows):
    columns = header[1:]

    if types_row is not None:
        return {col: ('numeric' if t.strip().lower() == 'numeric' else 'categorical')
                for col, t in zip(columns, types_row)}

    result = {}
    for i, col in enumerate(columns):
        values = [row[i + 1] for row in data_rows if len(row) > i + 1]
        result[col] = 'numeric' if values and all(_is_numeric(v) for v in values) else 'categorical'
    return result


def parse_metadata_columns(path):
    """Return {column_name: 'categorical'|'numeric'} for every non-id column.

    Honors an explicit `#q2:types` directive if present; otherwise infers a
    column as numeric only if every non-missing value in it parses as a
    float, matching QIIME2's own type-inference rule.
    """
    return _column_types(*read_metadata_rows(path))


def _table(header, data_rows):
    columns = header[1:]
    return {row[0]: dict(zip(columns, row[1:])) for row in data_rows}


def read_metadata_table(path):
    """Return {sample_id: {column: value}}."""
    header, _types_row, data_rows = read_metadata_rows(path)
    return _table(header, data_rows)


def _class_sizes(table, column, sample_ids):
    """Missing (empty) values are not a group: QIIME2 drops those samples
    from group-significance tests and sample classification. Counting ''
    as a class made a column with one real group plus two or more blanks
    look "eligible", and `qiime diversity beta-group-significance` then
    failed the whole run on it."""
    counts = defaultdict(int)
    for sample_id in sample_ids:
        value = table.get(sample_id, {}).get(column, '')
        if value != '':
            counts[value] += 1
    return dict(counts)


def class_sizes(path, column, sample_ids):
    """{value: count} for `column`, restricted to `sample_ids`."""
    return _class_sizes(read_metadata_table(path), column, sample_ids)


def _categorical_class_sizes(path, sample_ids):
    """Yield (column, {value: count}) for every categorical column, in file
    order, from a single read of `path` -- eligible_categorical_columns()
    and has_alpha_group_significance_column() previously re-read and
    re-parsed the whole file once per categorical column (via class_sizes()
    -> read_metadata_table()) on top of parse_metadata_columns()'s own read.
    """
    header, types_row, data_rows = read_metadata_rows(path)
    table = _table(header, data_rows)
    for col, col_type in _column_types(header, types_row, data_rows).items():
        if col_type == 'categorical':
            yield col, _class_sizes(table, col, sample_ids)


def final_sample_ids(sample_frequencies, fallback_sample_ids):
    """Sample ids to use for group-eligibility calculations (class_sizes,
    eligible_categorical_columns, etc.): every sample_frequencies key with a
    nonzero frequency, or every id in `fallback_sample_ids` if
    sample_frequencies is empty (its export is missing -- an older output
    folder, or a --skip-advanced-stats run).

    Deliberately keyed on sample_frequencies being empty, not on the
    filtered result being empty: a near-empty input sample can survive
    DADA2 as a zero-read row (retain-all-samples defaults to True), and
    those must be excluded from group-eligibility counts. If every sample
    happened to end up zero-read, falling back to "use every sample
    anyway" would silently un-exclude them all instead of correctly
    reporting nothing eligible.
    """
    if not sample_frequencies:
        return list(fallback_sample_ids)
    return [sid for sid, freq in sample_frequencies.items() if freq > 0]


def eligible_categorical_columns(path, sample_ids, min_per_group=2):
    """Categorical columns with >=2 distinct values, each held by at least
    `min_per_group` of the given `sample_ids`. Preserves metadata-file column
    order.
    """
    eligible = []
    for col, counts in _categorical_class_sizes(path, sample_ids):
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
    return any(len(counts) >= 2 and any(n >= 2 for n in counts.values())
               for _col, counts in _categorical_class_sizes(path, sample_ids))
