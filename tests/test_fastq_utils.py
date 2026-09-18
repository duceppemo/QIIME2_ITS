import gzip

import pytest

from qiime2_its import fastq_utils


def _gzip_xfl(path):
    """The gzip header's XFL byte (RFC 1952): zlib writes 2 for its slowest,
    maximum-compression level 9, 4 for fastest (level 1), 0 otherwise -- so a
    written file records whether it came out of gzip.open()'s default level 9."""
    return path.read_bytes()[8]


_GZIP_MAX_COMPRESSION_XFL = 2


class TestIsEmptyFastq:
    def test_true_for_zero_read_gzipped_file(self, tmp_path):
        fq = tmp_path / 'empty.fastq.gz'
        with gzip.open(fq, 'wt'):
            pass  # valid gzip stream, zero records -- a real "sample had 0 reads" case
        assert fastq_utils.is_empty_fastq(fq) is True

    def test_true_for_zero_read_plain_file(self, tmp_path):
        fq = tmp_path / 'empty.fastq'
        fq.write_text('')
        assert fastq_utils.is_empty_fastq(fq) is True

    def test_false_for_non_empty_file(self, tmp_path, write_fastq):
        fq = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
        write_fastq(fq, [('@r1', 'ACGT', 'IIII')], gz=True)
        assert fastq_utils.is_empty_fastq(fq) is False


class TestListFastq:
    def test_finds_all_accepted_extensions(self, tmp_path):
        names = ['a.fastq', 'b.fastq.gz', 'c.fq', 'd.fq.gz', 'e.txt']
        for name in names:
            (tmp_path / name).write_text('x')
        found = {p.name for p in fastq_utils.list_fastq(tmp_path)}
        assert found == {'a.fastq', 'b.fastq.gz', 'c.fq', 'd.fq.gz'}

    def test_recurses_into_subdirectories(self, tmp_path):
        sub = tmp_path / 'sub'
        sub.mkdir()
        (sub / 'nested.fastq.gz').write_text('x')
        found = fastq_utils.list_fastq(tmp_path)
        assert len(found) == 1
        assert found[0].name == 'nested.fastq.gz'

    def test_empty_directory_returns_empty_list(self, tmp_path):
        assert fastq_utils.list_fastq(tmp_path) == []


class TestParseFastqList:
    def test_pairs_r1_and_r2_by_sample(self):
        fastq_list = [
            'sampleA_bc1_L001_R2_001.fastq.gz',
            'sampleA_bc1_L001_R1_001.fastq.gz',
            'sampleB_bc2_L001_R1_001.fastq.gz',
        ]
        sample_dict = fastq_utils.parse_fastq_list(fastq_list)
        assert sample_dict['sampleA'] == ['sampleA_bc1_L001_R1_001.fastq.gz',
                                           'sampleA_bc1_L001_R2_001.fastq.gz']
        assert sample_dict['sampleB'] == ['sampleB_bc2_L001_R1_001.fastq.gz']


class TestValidateCasavaFilenames:
    def test_accepts_valid_names(self):
        fastq_utils.validate_casava_filenames([
            'L2S357_15_L001_R1_001.fastq.gz',
            'L2S357_15_L001_R2_001.fastq',
        ])  # no raise

    @pytest.mark.parametrize('bad_name', [
        'sample_barcode_L001_R1_001.txt',          # wrong extension
        'sample_barcode_L001_R1.fastq.gz',          # missing set number field
        'sample_barcode_001_R1_001.fastq.gz',       # lane doesn't start with L
        'sample_barcode_L001_R3_001.fastq.gz',      # bad direction
        'sample_barcode_L001_R1_002.fastq.gz',      # set number not 001
        'siteA_rep1_S1_L001_R1_001.fastq.gz',       # underscore inside the sample identifier itself
    ])
    def test_rejects_invalid_names(self, bad_name):
        with pytest.raises(ValueError):
            fastq_utils.validate_casava_filenames([bad_name])

    def test_hyphen_in_sample_identifier_is_accepted(self):
        """Real-world workaround for the underscore restriction above: hyphens
        inside the sample identifier are fine, since they don't affect the
        underscore-delimited field count."""
        fastq_utils.validate_casava_filenames(['siteA-rep1_S1_L001_R1_001.fastq.gz'])  # no raise

    def test_dot_in_sample_identifier_is_accepted(self):
        """Regression test: real SRA-derived sample identifiers (e.g. this
        one, taken from PRJNA767765) commonly contain dots. The old
        implementation split the whole filename on '.' to find the
        extension, which mistook the first dot in the sample identifier for
        the start of the extension and truncated it -- 'K.BeL.1.1_S1_L001_
        R1_001.fastq.gz' become stem 'K' (1 field, not 5) and was rejected
        as invalid even though it's a perfectly valid Casava name."""
        fastq_utils.validate_casava_filenames(['K.BeL.1.1_S1_L001_R1_001.fastq.gz'])  # no raise


class TestStripNonFastqFiles:
    def test_removes_manifest_and_metadata_keeps_fastq(self, tmp_path):
        """Regression test: `qiime tools export` writes MANIFEST/metadata.yml
        alongside the fastq files, which broke re-importing exported reads
        with CasavaOneEightSingleLanePerSampleDirFmt."""
        fq = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
        fq.write_text('x')
        manifest = tmp_path / 'MANIFEST'
        manifest.write_text('sample-id,filename,direction\n')
        metadata_yml = tmp_path / 'metadata.yml'
        metadata_yml.write_text('{}')

        fastq_utils.strip_non_fastq_files(tmp_path, keep=[fq])

        assert fq.exists()
        assert not manifest.exists()
        assert not metadata_yml.exists()


class TestRemoveEmptiesSe:
    def test_drops_empty_records_keeps_good_ones(self, tmp_path, write_fastq):
        fq = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
        write_fastq(fq, [
            ('@read1', 'ACGT', 'IIII'),
            ('@read2', '', ''),
            ('@read3', 'TTTT', 'JJJJ'),
        ], gz=True)

        stats = fastq_utils.remove_empties_se(fq)

        assert stats == (3, 2, 1)
        with gzip.open(fq, 'rt') as f:
            records = list(fastq_utils.iter_fastq_records(f))
        assert [r[0] for r in records] == ['@read1', '@read3']

    def test_handles_plain_text_input(self, tmp_path, write_fastq):
        fq = tmp_path / 'sample_bc_L001_R1_001.fastq'
        write_fastq(fq, [('@read1', 'ACGT', 'IIII')], gz=False)
        stats = fastq_utils.remove_empties_se(fq)
        assert stats == (1, 1, 0)
        assert fq.read_text().startswith('@read1')

    def test_does_not_rewrite_at_gzip_maximum_compression(self, tmp_path, write_fastq):
        """Regression test: gzip.open()'s default compresslevel=9 was ~98% of
        this function's runtime on a real 150k-read Illumina file (35s
        end-to-end, vs 1.9s at fastq_utils.GZIP_COMPRESS_LEVEL=4) for a ~10%
        smaller intermediate that `qiime tools import` reads exactly once.
        The fixture writes the input at the default level 9 (XFL=2), so the
        rewrite must come out at something else."""
        fq = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
        write_fastq(fq, [('@read1', 'ACGT', 'IIII')], gz=True)
        assert _gzip_xfl(fq) == _GZIP_MAX_COMPRESSION_XFL  # precondition: input is level 9

        fastq_utils.remove_empties_se(fq)

        assert _gzip_xfl(fq) != _GZIP_MAX_COMPRESSION_XFL


class TestRemoveEmptiesPe:
    def test_drops_pair_if_either_mate_empty(self, tmp_path, write_fastq):
        r1 = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
        r2 = tmp_path / 'sample_bc_L001_R2_001.fastq.gz'
        write_fastq(r1, [('@a', 'ACGT', 'IIII'), ('@b', '', ''), ('@c', 'GGGG', 'IIII')], gz=True)
        write_fastq(r2, [('@a', 'TTTT', 'IIII'), ('@b', 'CCCC', 'IIII'), ('@c', '', '')], gz=True)

        stats = fastq_utils.remove_empties_pe(r1, r2)

        assert stats == (3, 1, 2)
        with gzip.open(r1, 'rt') as f:
            r1_headers = [r[0] for r in fastq_utils.iter_fastq_records(f)]
        with gzip.open(r2, 'rt') as f:
            r2_headers = [r[0] for r in fastq_utils.iter_fastq_records(f)]
        assert r1_headers == r2_headers == ['@a']

    def test_does_not_rewrite_either_mate_at_gzip_maximum_compression(self, tmp_path, write_fastq):
        """Same as the single-end regression test, for both mates."""
        r1 = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
        r2 = tmp_path / 'sample_bc_L001_R2_001.fastq.gz'
        write_fastq(r1, [('@a', 'ACGT', 'IIII')], gz=True)
        write_fastq(r2, [('@a', 'TTTT', 'IIII')], gz=True)
        assert _gzip_xfl(r1) == _gzip_xfl(r2) == _GZIP_MAX_COMPRESSION_XFL

        fastq_utils.remove_empties_pe(r1, r2)

        assert _gzip_xfl(r1) != _GZIP_MAX_COMPRESSION_XFL
        assert _gzip_xfl(r2) != _GZIP_MAX_COMPRESSION_XFL


class TestReverseComplement:
    @pytest.mark.parametrize('seq,expected', [
        ('ACGT', 'ACGT'),
        ('AACCGGTT', 'AACCGGTT'),
        ('AAAA', 'TTTT'),
        ('ACGTN', 'NACGT'),
        ('acgtn', 'nacgt'),
    ])
    def test_reverse_complement(self, seq, expected):
        assert fastq_utils.reverse_complement(seq) == expected


class TestRcFastq:
    def test_writes_reverse_complemented_gzipped_copy(self, tmp_path, write_fastq):
        src = tmp_path / 'sample_bc_L001_R1_001.fastq'
        write_fastq(src, [('@r1', 'AACG', 'IIJJ')], gz=False)

        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        out_path = fastq_utils.rc_fastq(src, out_dir)

        assert out_path.name.endswith('.gz')
        with gzip.open(out_path, 'rt') as f:
            header, seq, plus, qual = next(fastq_utils.iter_fastq_records(f))
        assert header == '@r1'
        assert seq == 'CGTT'
        assert qual == 'JJII'

    def test_raises_on_malformed_header(self, tmp_path):
        src = tmp_path / 'bad_bc_L001_R1_001.fastq'
        src.write_text('not-a-header\nACGT\n+\nIIII\n')
        with pytest.raises(ValueError):
            fastq_utils.rc_fastq(src, tmp_path)

    def test_does_not_write_at_gzip_maximum_compression(self, tmp_path, write_fastq):
        """Same regression as TestRemoveEmptiesSe's: rc_fastq() on a real
        150k-read Illumina file took 36.5s end-to-end at gzip.open()'s default
        level 9, of which ~0.5s was reading/reverse-complementing."""
        src = tmp_path / 'sample_bc_L001_R1_001.fastq'
        write_fastq(src, [('@r1', 'AACG', 'IIJJ')], gz=False)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()

        out_path = fastq_utils.rc_fastq(src, out_dir)

        assert _gzip_xfl(out_path) != _GZIP_MAX_COMPRESSION_XFL


class TestRcFastqParallel:
    def test_reverse_complements_every_file(self, tmp_path, write_fastq):
        a = tmp_path / 'a_bc_L001_R1_001.fastq'
        b = tmp_path / 'b_bc_L001_R1_001.fastq.gz'
        write_fastq(a, [('@a', 'AACG', 'IIJJ')], gz=False)
        write_fastq(b, [('@b', 'TTTA', 'IIJJ')], gz=True)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()

        fastq_utils.rc_fastq_parallel([a, b], out_dir, parallel=2)

        assert sorted(p.name for p in out_dir.iterdir()) == ['a_bc_L001_R1_001.fastq.gz',
                                                             'b_bc_L001_R1_001.fastq.gz']
        with gzip.open(out_dir / 'b_bc_L001_R1_001.fastq.gz', 'rt') as f:
            assert next(fastq_utils.iter_fastq_records(f))[1] == 'TAAA'

    def test_rejects_same_named_inputs_from_different_subfolders(self, tmp_path, write_fastq):
        """list_fastq() searches recursively, but every output lands flat in
        output_dir under its input's basename -- so two same-named inputs
        from different subfolders would be reverse-complemented into the
        *same* output file by two worker threads at once, each truncating
        and overwriting the other's bytes into a corrupt gzip stream that
        holds neither file's reads. Must be refused before any worker
        starts, not discovered as a garbled file at import time."""
        name = 'sample_bc_L001_R1_001.fastq'
        a = tmp_path / 'run1' / name
        b = tmp_path / 'run2' / name
        a.parent.mkdir()
        b.parent.mkdir()
        write_fastq(a, [('@a', 'AACG', 'IIJJ')], gz=False)
        write_fastq(b, [('@b', 'TTTA', 'IIJJ')], gz=False)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()

        with pytest.raises(ValueError, match='sample_bc_L001_R1_001.fastq.gz'):
            fastq_utils.rc_fastq_parallel([a, b], out_dir, parallel=2)

        assert list(out_dir.iterdir()) == []  # refused up front, nothing written

    def test_rejects_plain_and_gzipped_copies_of_the_same_name(self, tmp_path, write_fastq):
        """a.fastq and a.fastq.gz both map to the output name a.fastq.gz."""
        a = tmp_path / 'sample_bc_L001_R1_001.fastq'
        a_gz = tmp_path / 'sample_bc_L001_R1_001.fastq.gz'
        write_fastq(a, [('@a', 'AACG', 'IIJJ')], gz=False)
        write_fastq(a_gz, [('@a', 'AACG', 'IIJJ')], gz=True)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()

        with pytest.raises(ValueError):
            fastq_utils.rc_fastq_parallel([a, a_gz], out_dir, parallel=2)
