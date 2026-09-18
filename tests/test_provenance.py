"""Unit tests for provenance.py: parsing `qiime info`'s real output shape,
building the sample-file manifest, and assembling run_metadata.json's
contents. No subprocess/QIIME2 needed -- qiime_wrapper.qiime_info() is what
actually runs `qiime info`; this only covers what happens to its text.
"""
import json
from datetime import datetime, timezone

from qiime2_its import provenance

# Trimmed but verbatim-shaped copy of real `qiime info` output from
# rachis-qiime2-2026.7 (see report_data.py's module docstring convention:
# internal formats are confirmed against real output, not guessed).
_REAL_QIIME_INFO_TEXT = """System versions
Python version: 3.12.13
Parsl version: 2026.2.23
rachis release: 2026.7
rachis version: 2026.7.0
q2cli version: 2026.7.0

Installed plugins
alignment: 2026.7.0
dada2: 2026.7.0
itsxpress: 2.2.0
taxa: 2026.7.0

Application config directory
/home/bioinfo/miniconda3/envs/rachis-qiime2-2026.7/var/q2cli

Config
Config Source: /home/bioinfo/miniconda3/envs/rachis-qiime2-2026.7/etc/qiime2_config.toml

Getting help
To find help and learning resources, visit https://qiime2.org.
"""

# Verbatim shape of real `bbduk.sh --version` output (written to stderr, not stdout).
_REAL_BBDUK_VERSION_TEXT = (
    'java -ea  --add-modules jdk.incubator.vector -Xmx95681m -Xms95681m '
    '-cp /opt/bbmap-39.80-0/bbtools.jar bbduk.BBDukS --version\n'
    'WARNING: Using incubator modules: jdk.incubator.vector\n'
    'BBTools version 39.80\n'
    'For help, please run the shellscript with no parameters, or look in /docs/.\n'
)


class TestFormatCommandLine:
    def test_one_flag_and_value_per_line(self):
        argv = ['qiime2-its', '-q', 'rachis-qiime2-2026.7', '-i', 'in', '-o', 'out',
                '-m', 'meta.tsv', '-c', 'clf.qza', '-pe', '--extract-its2', '--taxa', 'Fungi']
        assert provenance.format_command_line(argv) == (
            'qiime2-its \\\n'
            '    -q rachis-qiime2-2026.7 \\\n'
            '    -i in \\\n'
            '    -o out \\\n'
            '    -m meta.tsv \\\n'
            '    -c clf.qza \\\n'
            '    -pe \\\n'
            '    --extract-its2 \\\n'
            '    --taxa Fungi'
        )

    def test_program_only(self):
        assert provenance.format_command_line(['qiime2-its']) == 'qiime2-its'

    def test_empty_argv(self):
        assert provenance.format_command_line([]) == ''

    def test_quotes_values_needing_it(self):
        result = provenance.format_command_line(['qiime2-its', '-m', 'a file with spaces.tsv'])
        assert "'a file with spaces.tsv'" in result


class TestParseBbdukVersion:
    def test_extracts_version_number(self):
        assert provenance.parse_bbduk_version(_REAL_BBDUK_VERSION_TEXT) == '39.80'

    def test_none_when_not_installed(self):
        assert provenance.parse_bbduk_version(None) is None

    def test_none_when_text_does_not_match(self):
        assert provenance.parse_bbduk_version('command not found') is None


class TestParseQiimeInfo:
    def test_parses_system_versions(self):
        result = provenance.parse_qiime_info(_REAL_QIIME_INFO_TEXT)
        assert result['system']['Python version'] == '3.12.13'
        assert result['system']['rachis version'] == '2026.7.0'

    def test_parses_installed_plugins(self):
        result = provenance.parse_qiime_info(_REAL_QIIME_INFO_TEXT)
        assert result['plugins'] == {
            'alignment': '2026.7.0', 'dada2': '2026.7.0', 'itsxpress': '2.2.0', 'taxa': '2026.7.0',
        }

    def test_ignores_trailing_non_key_value_sections(self):
        # "Config Source: ..." looks like a key:value line but belongs to
        # the "Config" section, which isn't "System versions" or "Installed
        # plugins" -- it must not leak into either parsed dict.
        result = provenance.parse_qiime_info(_REAL_QIIME_INFO_TEXT)
        assert 'Config Source' not in result['system']
        assert 'Config Source' not in result['plugins']

    def test_empty_text_returns_empty_dicts(self):
        assert provenance.parse_qiime_info('') == {'system': {}, 'plugins': {}}


class TestSampleFileManifest:
    def test_extracts_filenames_sorted_by_sample(self):
        manifest = provenance.sample_file_manifest({
            'sampleB': ['/data/sampleB_S2_L001_R1_001.fastq.gz'],
            'sampleA': ['/data/sampleA_S1_L001_R1_001.fastq.gz',
                        '/data/sampleA_S1_L001_R2_001.fastq.gz'],
        })
        assert manifest == [
            {'sample_id': 'sampleA',
             'files': ['sampleA_S1_L001_R1_001.fastq.gz', 'sampleA_S1_L001_R2_001.fastq.gz']},
            {'sample_id': 'sampleB', 'files': ['sampleB_S2_L001_R1_001.fastq.gz']},
        ]


class TestBuildRunMetadata:
    def _build(self, **overrides):
        kwargs = dict(
            qiime2_its_version='0.2.0',
            command_line='qiime2-its -q env -i in -o out -m meta.tsv -c clf.qza -pe',
            start_time=datetime(2026, 9, 16, 12, 0, 0, tzinfo=timezone.utc),
            end_time=datetime(2026, 9, 16, 12, 2, 30, tzinfo=timezone.utc),
            username='bioinfo', hostname='workstation',
            platform_string='Linux-x86_64', conda_env='rachis-qiime2-2026.7',
            qiime_info_text=_REAL_QIIME_INFO_TEXT,
            bbduk_version_text=_REAL_BBDUK_VERSION_TEXT,
            input_folder='/in', metadata_file='meta.tsv', classifier_file='clf.qza',
            output_folder='/out',
            sample_dict={'sampleA': ['/in/sampleA_S1_L001_R1_001.fastq.gz']},
            parameters={'max_ee': 4.0, 'allow_one_off': True},
        )
        kwargs.update(overrides)
        return provenance.build_run_metadata(**kwargs)

    def test_computes_duration_from_start_and_end(self):
        assert self._build()['pipeline']['duration_seconds'] == 150.0

    def test_pulls_framework_version_from_qiime_info(self):
        assert self._build()['environment']['qiime2_framework_version'] == '2026.7.0'

    def test_pulls_bbmap_version(self):
        assert self._build()['environment']['bbmap_version'] == '39.80'

    def test_bbmap_version_none_when_not_installed(self):
        assert self._build(bbduk_version_text=None)['environment']['bbmap_version'] is None

    def test_carries_plugin_versions_through(self):
        assert self._build()['qiime2_plugins']['dada2'] == '2026.7.0'

    def test_builds_sample_manifest(self):
        assert self._build()['inputs']['samples'] == [
            {'sample_id': 'sampleA', 'files': ['sampleA_S1_L001_R1_001.fastq.gz']}
        ]

    def test_stores_parameters_as_given(self):
        assert self._build()['parameters'] == {'max_ee': 4.0, 'allow_one_off': True}

    def test_result_is_json_round_trippable(self, tmp_path):
        result = self._build()
        path = tmp_path / 'run_metadata.json'
        provenance.write_run_metadata(path, result)
        assert json.loads(path.read_text()) == result

    def test_writes_atomically_leaving_no_leftover_temp_file(self, tmp_path):
        """Regression test: a plain write_text() used to leave a truncated,
        unparseable run_metadata.json on disk if the process was killed or
        the disk filled up mid-write -- report_data.parse_run_metadata()
        would then hit a corrupt file it couldn't distinguish from "never
        finished writing". Written via a temp file + rename instead, so the
        target path is always either the previous complete file or the new
        one, never a partial write."""
        path = tmp_path / 'run_metadata.json'
        provenance.write_run_metadata(path, self._build())
        assert list(tmp_path.iterdir()) == [path]  # no leftover .tmp sibling
