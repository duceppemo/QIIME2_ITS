"""Pure-Python helpers for listing, parsing, and cleaning fastq files.

None of these functions shell out to external tools, so they are fully
unit-testable without QIIME2/ITSxpress/BBTools installed.
"""
import gzip
import os
from collections import Counter, defaultdict, namedtuple
from concurrent import futures
from pathlib import Path

FASTQ_EXTENSIONS = ('.fastq', '.fastq.gz', '.fq', '.fq.gz')

# gzip.open()'s default compresslevel is 9. Measured on a real 150k-read
# Illumina fastq (77 MB uncompressed): level 9 took 35s to write vs 9.5s at
# level 6 (gzip(1)'s default) and 1.9s at level 4, for 23.5 / 24.3 / 25.9 MB
# outputs -- and reading+parsing the same file takes 0.5s, so at level 9
# remove_empties_se()/rc_fastq() spent ~98% of their time in zlib. Every
# file written here is either a pipeline intermediate `qiime tools import`
# reads exactly once (exported_reads/, rc_reads/) or the standalone fastq-rc
# tool's output, so a ~10% smaller file isn't worth ~18x the write time. 4 is
# also bcl2fastq's own default for the fastq files it writes.
GZIP_COMPRESS_LEVEL = 4

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


def strip_non_fastq_files(directory, keep):
    """Delete every file in `directory` that isn't in `keep`.

    `qiime tools export` of a demultiplexed-reads artifact writes a MANIFEST
    and metadata.yml alongside the fastq files; CasavaOneEightSingleLanePerSampleDirFmt
    rejects re-importing a directory containing anything but Casava-named fastq
    files, so these must be stripped before re-importing exported/post-processed reads.
    """
    keep = {Path(p) for p in keep}
    for entry in Path(directory).iterdir():
        if entry.is_file() and entry not in keep:
            entry.unlink()


def validate_casava_filenames(fastq_list):
    """Raise ValueError if any file doesn't follow the QIIME2 Casava naming scheme:
    ``<sample>_<barcode>_L<lane>_R[12]_001.fastq[.gz]``.
    """
    for fq in fastq_list:
        name = Path(fq).name
        if name.endswith('.fastq.gz'):
            stem = name[:-len('.fastq.gz')]
        elif name.endswith('.fastq'):
            stem = name[:-len('.fastq')]
        else:
            raise ValueError(CASAVA_NAMING_ERROR)

        # Split only the trailing .fastq[.gz] off first -- real-world sample
        # identifiers (e.g. SRA-derived ones like "K.BeL.1.1") can themselves
        # contain dots, which a naive name.split('.')[0] would mistake for
        # the start of the extension and truncate the sample identifier.
        fields = stem.split('_')
        if len(fields) != 5 or not fields[2].startswith('L') or fields[3] not in ('R1', 'R2') \
                or fields[4] != '001':
            raise ValueError(CASAVA_NAMING_ERROR)


def _open_fastq(path, mode='rt', gz=None):
    """Open a fastq file for text read/write, transparently gzip-aware.

    `gz` overrides the by-extension detection, for a temp file whose own name
    doesn't end in .gz but must use its final destination's compression
    (remove_empties_se/pe's `.clean` files). Writes use GZIP_COMPRESS_LEVEL.
    """
    if gz is None:
        gz = str(path).endswith('.gz')
    if not gz:
        return open(path, mode)
    if 'w' in mode:
        return gzip.open(path, mode, compresslevel=GZIP_COMPRESS_LEVEL)
    return gzip.open(path, mode)


def is_empty_fastq(path):
    """True if `path` contains zero fastq records.

    Only reads the first line, so this is cheap even for large real files --
    a zero-read sample otherwise fails several steps deep (ITSxpress's HMM
    search errors out on an empty/misformatted input file) with a cryptic
    external-tool stack trace instead of a clear, immediate message.
    """
    with _open_fastq(path, 'rt') as f:
        return f.readline() == ''


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
    gz = fastq_path.name.endswith('.gz')

    total = kept = empty = 0
    with _open_fastq(fastq_path, 'rt', gz=gz) as in_f, _open_fastq(tmp_path, 'wt', gz=gz) as out_f:
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
    gz_r1 = r1_path.name.endswith('.gz')
    gz_r2 = r2_path.name.endswith('.gz')

    total = kept = empty = 0
    with _open_fastq(r1_path, 'rt', gz=gz_r1) as in_r1, _open_fastq(r2_path, 'rt', gz=gz_r2) as in_r2, \
            _open_fastq(tmp_r1, 'wt', gz=gz_r1) as out_r1, _open_fastq(tmp_r2, 'wt', gz=gz_r2) as out_r2:
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


def _rc_output_path(fastq_path, output_dir):
    """Where rc_fastq() writes `fastq_path`'s copy: flat in `output_dir`, gzipped."""
    fastq_path = Path(fastq_path)
    out_name = fastq_path.name if fastq_path.name.endswith('.gz') else fastq_path.name + '.gz'
    return Path(output_dir) / out_name


def rc_fastq(fastq_path, output_dir):
    """Write a reverse-complemented, gzip-compressed copy of `fastq_path` into `output_dir`."""
    fastq_path = Path(fastq_path)
    out_path = _rc_output_path(fastq_path, output_dir)

    print(f'\t{fastq_path.name}')
    with _open_fastq(fastq_path, 'rt') as in_f, _open_fastq(out_path, 'wt', gz=True) as out_f:
        for header, sequence, plus, quality in iter_fastq_records(in_f):
            if not header.startswith('@'):
                raise ValueError(f'Invalid fastq file: {fastq_path}')
            out_f.write(f'{header}\n{reverse_complement(sequence)}\n{plus}\n{quality[::-1]}\n')
    return out_path


def rc_fastq_parallel(fastq_list, output_dir, parallel):
    # list_fastq() searches recursively, but every output lands flat in
    # `output_dir` under its input's basename -- so two same-named inputs
    # from different subfolders (run1/S1_..._R1_001.fastq.gz and
    # run2/S1_..._R1_001.fastq.gz, or a.fastq next to a.fastq.gz) would be
    # written to the *same* output path by two workers at once, each
    # truncating and overwriting the other's bytes: a corrupt gzip stream
    # containing neither file's reads. Refused before any worker starts.
    out_names = Counter(_rc_output_path(fq, output_dir).name for fq in fastq_list)
    duplicates = sorted(name for name, n in out_names.items() if n > 1)
    if duplicates:
        raise ValueError(
            'Two or more input fastq files (in different subfolders, or a .fastq next to its .fastq.gz) '
            'would be written to the same output file name -- rename them or process the subfolders '
            'separately: {}'.format(', '.join(duplicates)))
    with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
        list(executor.map(lambda fq: rc_fastq(fq, output_dir), fastq_list))
