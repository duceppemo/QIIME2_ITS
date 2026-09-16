"""QA/audit provenance for a pipeline run: who ran it, when, with what exact
command and parameters, against which inputs, and against which QIIME2/
plugin versions.

cli/pipeline.py collects the raw inputs (subprocess output, environment,
timing) and calls build_run_metadata() here to assemble them into one
JSON-serializable dict, written to run_metadata.json in the output folder.
report.py reads it back (via report_data.parse_run_metadata()) to build the
report's QA pages. Kept pure/testable the same way report_data.py is: the
only subprocess call this depends on (`qiime info`) happens in
qiime_wrapper.py, which just hands this module the raw text to parse.
"""
import json
import re
import shlex
from pathlib import Path

_QIIME_INFO_SECTION_HEADERS = {'System versions': 'system', 'Installed plugins': 'plugins'}
_BBDUK_VERSION_RE = re.compile(r'BBTools version (\S+)')


def format_command_line(argv):
    """Format `sys.argv` as a readable, multi-line shell invocation, one
    flag (and its value, if any) per line -- e.g.:

        qiime2-its \\
            -q rachis-qiime2-2026.7 \\
            -i /path/to/input \\
            -pe

    A single `shlex.join(argv)`-style line reads fine on a terminal but
    wraps illegibly once it's justified inside a report table cell.
    """
    if not argv:
        return ''
    program, *args = argv
    lines = [shlex.quote(program)]
    current = []
    for arg in args:
        if arg.startswith('-') and current:
            lines.append(' '.join(current))
            current = []
        current.append(shlex.quote(arg))
    if current:
        lines.append(' '.join(current))
    return ' \\\n    '.join(lines)


def parse_qiime_info(text):
    """Parse `qiime info`'s stdout into {'system': {...}, 'plugins': {...}}.

    Format confirmed against real `qiime info` output from
    rachis-qiime2-2026.7, not assumed: a "System versions" section and an
    "Installed plugins" section of "key: value" lines, followed by
    "Application config directory"/"Config"/"Getting help" sections this
    doesn't need, whose content lines (bare paths/URLs) don't contain ":".
    """
    system = {}
    plugins = {}
    section = None
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line in _QIIME_INFO_SECTION_HEADERS:
            section = _QIIME_INFO_SECTION_HEADERS[line]
            continue
        if ':' not in line:
            # A bare header ("Config", "Getting help", ...) or one of the
            # plain paths/URLs that follow them -- either way, not part of
            # the section we're collecting.
            section = None
            continue
        if section:
            key, _, value = line.partition(':')
            (system if section == 'system' else plugins)[key.strip()] = value.strip()
    return {'system': system, 'plugins': plugins}


def parse_bbduk_version(text):
    """Extract "39.80" from BBTools' "BBTools version 39.80" banner (its
    `--version` output, which is on stderr -- see size_filter.bbduk_version()).
    None if `text` is None (bbduk.sh isn't installed) or doesn't match.
    """
    if not text:
        return None
    match = _BBDUK_VERSION_RE.search(text)
    return match.group(1) if match else None


def sample_file_manifest(sample_dict):
    """{sample: [Path, ...]} -> [{'sample_id': ..., 'files': [name, ...]}, ...], sorted by sample."""
    return [
        {'sample_id': sample, 'files': [Path(f).name for f in files]}
        for sample, files in sorted(sample_dict.items())
    ]


def build_run_metadata(*, qiime2_its_version, command_line, start_time, end_time,
                        username, hostname, platform_string, conda_env, qiime_info_text,
                        bbduk_version_text, input_folder, metadata_file, classifier_file,
                        output_folder, sample_dict, parameters):
    """Assemble one JSON-serializable provenance record for this run.

    `start_time`/`end_time` are timezone-aware datetimes; `qiime_info_text`
    is `qiime info`'s raw stdout (parsed here via parse_qiime_info());
    `bbduk_version_text` is `bbduk.sh --version`'s raw stderr, or None if it
    isn't installed (parsed via parse_bbduk_version() -- BBTools/BBMap is an
    external, non-QIIME2 dependency `qiime info` doesn't know about, only
    needed for --min-len/--max-len); `parameters` is the run's CLI arguments
    as a plain dict (e.g. vars(argparse.Namespace)), stored as-is for exact
    reproducibility.
    """
    qiime_info = parse_qiime_info(qiime_info_text)
    return {
        'pipeline': {
            'qiime2_its_version': qiime2_its_version,
            'command_line': command_line,
            'start_time': start_time.isoformat(),
            'end_time': end_time.isoformat(),
            'duration_seconds': round((end_time - start_time).total_seconds(), 1),
        },
        'environment': {
            'username': username,
            'hostname': hostname,
            'platform': platform_string,
            'conda_env': conda_env,
            'python_version': qiime_info['system'].get('Python version'),
            'qiime2_framework_version': qiime_info['system'].get('rachis version')
            or qiime_info['system'].get('q2cli version'),
            'bbmap_version': parse_bbduk_version(bbduk_version_text),
        },
        'qiime2_plugins': qiime_info['plugins'],
        'inputs': {
            'input_folder': str(input_folder),
            'metadata_file': str(metadata_file),
            'classifier_file': str(classifier_file),
            'output_folder': str(output_folder),
            'samples': sample_file_manifest(sample_dict),
        },
        'parameters': parameters,
    }


def write_run_metadata(path, run_metadata):
    Path(path).write_text(json.dumps(run_metadata, indent=2, default=str))
