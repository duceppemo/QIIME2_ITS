"""Pure-Python helpers for listing, parsing, and cleaning fastq files.

None of these functions shell out to external tools, so they are fully
unit-testable without QIIME2/ITSxpress/BBTools installed.
"""
import gzip
import os
from collections import defaultdict, namedtuple
from concurrent import futures
from pathlib import Path

FASTQ_EXTENSIONS = ('.fastq', '.fastq.gz', '.fq', '.fq.gz')

CleanStats = namedtuple('CleanStats', ['total', 'kept', 'empty'])

CASAVA_NAMING_ERROR = (
    'File name must be as is:\n'
    '\t1. the sample identifier,\n'
    '\t2. the barcode sequence or a barcode identifier,\n'
    '\t3. the lane number starting with "L" followed by 3 digits,\n'
    '\t4. the direction of the read ("R1" or "R2"), and\n'
    '\t5. the set number (always "001").'
)

# Handles upper and lower case, and ambiguous "N"/"n".
_COMPLEMENT_TABLE = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def list_fastq(directory):
    """Recursively list fastq files (.fastq[.gz]/.fq[.gz]) under `directory`."""
    fastq_list = []
    for root, _dirs, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(FASTQ_EXTENSIONS):
                fastq_list.append(Path(root) / filename)
    return sorted(fastq_list)


def parse_fastq_list(fastq_list):
    """Bucket fastq files by sample name using the QIIME2 Casava naming scheme.

    ``<sample>_<barcode>_L<lane>_R[12]_001.fastq[.gz]`` -> {sample: [R1, (R2)]}
    """
    sample_dict = defaultdict(list)
    for fq in fastq_list:
        fields = Path(fq).name.split('_')
        sample = fields[0]
        direction = fields[3]
        if 'R1' in direction:
            sample_dict[sample].insert(0, fq)
        elif 'R2' in direction:
            sample_dict[sample].insert(1, fq)
    return dict(sample_dict)


def validate_casava_filenames(fastq_list):
    """Raise ValueError if any file doesn't follow the QIIME2 Casava naming scheme:
    ``<sample>_<barcode>_L<lane>_R[12]_001.fastq[.gz]``.
    """
    for fq in fastq_list:
        name = Path(fq).name
        parts = name.split('.')
        ext = parts[-2] if name.endswith('.gz') else parts[-1]
        if ext != 'fastq':
            raise ValueError(CASAVA_NAMING_ERROR)

        fields = parts[0].split('_')
        if len(fields) != 5 or not fields[2].startswith('L') or fields[3] not in ('R1', 'R2') \
                or fields[4] != '001':
            raise ValueError(CASAVA_NAMING_ERROR)


def _open_fastq(path, mode='rt'):
    """Open a fastq file for text read/write, transparently gzip-aware."""
    opener = gzip.open if str(path).endswith('.gz') else open
    return opener(path, mode)


def iter_fastq_records(file_handle):
    """Yield (header, sequence, plus, quality) 4-tuples from an open fastq handle."""
    while True:
        header = file_handle.readline().rstrip('\n')
        if not header:
            return
        sequence = file_handle.readline().rstrip('\n')
        plus = file_handle.readline().rstrip('\n')
        quality = file_handle.readline().rstrip('\n')
        yield header, sequence, plus, quality


def remove_empties_se(fastq_path):
    """Remove fastq records with an empty sequence (an ITSxpress artifact).

    Rewrites the file in place (same compression as the input) and returns a
    CleanStats(total, kept, empty) summary.
    """
    fastq_path = Path(fastq_path)
    tmp_path = fastq_path.with_name(fastq_path.name + '.clean')
    # tmp_path's own name doesn't end in .gz, so pick the (de)compressor from
    # the real fastq_path for both sides instead of letting _open_fastq guess.
    opener = gzip.open if fastq_path.name.endswith('.gz') else open

    total = kept = empty = 0
    with opener(fastq_path, 'rt') as in_f, opener(tmp_path, 'wt') as out_f:
        for header, sequence, plus, quality in iter_fastq_records(in_f):
            total += 1
            if sequence == '':
                empty += 1
                continue
            out_f.write(f'{header}\n{sequence}\n{plus}\n{quality}\n')
            kept += 1

    tmp_path.replace(fastq_path)
    print(f'{fastq_path.name}: total sequences: {total}, good sequences: {kept}, empty sequences: {empty}')
    return CleanStats(total, kept, empty)


def remove_empties_pe(r1_path, r2_path):
    """Remove read pairs where either mate has an empty sequence, keeping R1/R2 in sync."""
    r1_path, r2_path = Path(r1_path), Path(r2_path)
    tmp_r1 = r1_path.with_name(r1_path.name + '.clean')
    tmp_r2 = r2_path.with_name(r2_path.name + '.clean')
    opener_r1 = gzip.open if r1_path.name.endswith('.gz') else open
    opener_r2 = gzip.open if r2_path.name.endswith('.gz') else open

    total = kept = empty = 0
    with opener_r1(r1_path, 'rt') as in_r1, opener_r2(r2_path, 'rt') as in_r2, \
            opener_r1(tmp_r1, 'wt') as out_r1, opener_r2(tmp_r2, 'wt') as out_r2:
        for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(iter_fastq_records(in_r1), iter_fastq_records(in_r2)):
            total += 1
            if s1 == '' or s2 == '':
                empty += 1
                continue
            out_r1.write(f'{h1}\n{s1}\n{p1}\n{q1}\n')
            out_r2.write(f'{h2}\n{s2}\n{p2}\n{q2}\n')
            kept += 1

    tmp_r1.replace(r1_path)
    tmp_r2.replace(r2_path)
    print(f'{r1_path.name.split("_")[0]}: input sequences: {total}, good: {kept}, empty: {empty}')
    return CleanStats(total, kept, empty)


def remove_empties_se_parallel(fastq_list, parallel):
    with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
        list(executor.map(remove_empties_se, fastq_list))


def remove_empties_pe_parallel(sample_dict, parallel):
    with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
        pairs = ((reads[0], reads[1]) for reads in sample_dict.values())
        list(executor.map(lambda p: remove_empties_pe(*p), pairs))


def reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT_TABLE)[::-1]


def rc_fastq(fastq_path, output_dir):
    """Write a reverse-complemented, gzip-compressed copy of `fastq_path` into `output_dir`."""
    fastq_path = Path(fastq_path)
    output_dir = Path(output_dir)
    out_name = fastq_path.name if fastq_path.name.endswith('.gz') else fastq_path.name + '.gz'
    out_path = output_dir / out_name

    print(f'\t{fastq_path.name}')
    with _open_fastq(fastq_path, 'rt') as in_f, gzip.open(out_path, 'wt') as out_f:
        for header, sequence, plus, quality in iter_fastq_records(in_f):
            if not header.startswith('@'):
                raise ValueError(f'Invalid fastq file: {fastq_path}')
            out_f.write(f'{header}\n{reverse_complement(sequence)}\n{plus}\n{quality[::-1]}\n')
    return out_path


def rc_fastq_parallel(fastq_list, output_dir, parallel):
    with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
        list(executor.map(lambda fq: rc_fastq(fq, output_dir), fastq_list))
