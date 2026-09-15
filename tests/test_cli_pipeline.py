import pytest

from qiime2_its.cli.pipeline import build_parser

REQUIRED = ['-q', 'rachis-qiime2-2026.7', '-i', 'in/', '-o', 'out/', '-m', 'meta.tsv', '-c', 'clf.qza']


def test_requires_se_or_pe():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED)


def test_se_and_pe_are_mutually_exclusive():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED + ['-se', '-pe'])


def test_se_alone_parses(capsys):
    parser = build_parser()
    args = parser.parse_args(REQUIRED + ['-se'])
    assert args.se is True
    assert args.pe is False


def test_extract_its1_and_its2_are_mutually_exclusive():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED + ['-se', '--extract-its1', '--extract-its2'])


def test_defaults():
    parser = build_parser()
    args = parser.parse_args(REQUIRED + ['-se'])
    assert args.threads == 4
    assert args.parallel_processes == 1
    assert args.min_len == 0
    assert args.max_len == 0
    assert args.taxa == 'Fungi'
    assert args.reverse_complement is False


def test_missing_required_argument_exits():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(['-se'])


def test_invalid_taxa_choice_exits():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(REQUIRED + ['-se', '--taxa', 'NotARealTaxon'])
