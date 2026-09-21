import gzip

import pytest

from qiime2_its import fastq_utils
from qiime2_its.cli import fastq_rc


def test_raises_if_output_equals_input(tmp_path):
    with pytest.raises(ValueError):
        fastq_rc.run(tmp_path, tmp_path, threads=1)


def test_raises_if_output_is_inside_input(tmp_path):
    """The fastq search is recursive: a nested output folder's files would
    be found (and re-processed) as input on the next run."""
    with pytest.raises(ValueError, match='not inside'):
        fastq_rc.run(tmp_path, tmp_path / 'rc', threads=1)


def test_raises_if_input_missing(tmp_path):
    with pytest.raises(ValueError):
        fastq_rc.run(tmp_path / 'does-not-exist', tmp_path / 'out', threads=1)


def test_raises_if_no_fastq_files(tmp_path):
    in_dir = tmp_path / 'in'
    in_dir.mkdir()
    out_dir = tmp_path / 'out'
    with pytest.raises(ValueError):
        fastq_rc.run(in_dir, out_dir, threads=1)


def test_reverse_complements_all_files(tmp_path, write_fastq):
    in_dir = tmp_path / 'in'
    in_dir.mkdir()
    out_dir = tmp_path / 'out'
    write_fastq(in_dir / 'sample_bc_L001_R1_001.fastq', [('@r1', 'AACG', 'IIJJ')], gz=False)

    fastq_rc.run(in_dir, out_dir, threads=2)

    out_files = fastq_utils.list_fastq(out_dir)
    assert len(out_files) == 1
    with gzip.open(out_files[0], 'rt') as f:
        _header, seq, _plus, _qual = next(fastq_utils.iter_fastq_records(f))
    assert seq == 'CGTT'
