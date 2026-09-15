import gzip

import pytest


def _write_fastq(path, records, gz):
    lines = []
    for header, seq, qual in records:
        lines.extend([header, seq, '+', qual])
    text = '\n'.join(lines) + '\n'
    if gz:
        with gzip.open(path, 'wt') as f:
            f.write(text)
    else:
        with open(path, 'w') as f:
            f.write(text)


@pytest.fixture
def write_fastq():
    """Write a small fastq file. Usage: write_fastq(path, [(header, seq, qual), ...], gz=True)."""
    return _write_fastq
