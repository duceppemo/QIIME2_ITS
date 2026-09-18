import gzip

import pytest

from qiime2_its import fastq_utils


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
