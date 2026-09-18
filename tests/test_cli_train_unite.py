from pathlib import Path

import pytest

from qiime2_its.cli import train_unite


class TestParseUniteFilename:
    def test_parses_version_clustering_release_date(self):
        version, clustering, release_date = train_unite.parse_unite_filename(
            'sh_refs_qiime_ver10_99_19.02.2025.fasta')
        assert version == 'ver10'
        assert clustering == '99'
        assert release_date == '19.02.2025'

    def test_extension_does_not_leak_into_release_date(self):
        """Regression test: parsing used to split the filename (with its
        .fasta extension still attached) before extracting fields, so the
        release date field ended up as "19.02.2025.fasta"."""
        _version, _clustering, release_date = train_unite.parse_unite_filename(
            'sh_refs_qiime_ver10_99_19.02.2025.fasta')
        assert not release_date.endswith('.fasta')


class TestFindUniteFiles:
    def test_finds_developer_refs_and_taxonomy_files(self, tmp_path):
        developer = tmp_path / 'sh_qiime_release_19.02.2025' / 'developer'
        developer.mkdir(parents=True)
        seq = developer / 'sh_refs_qiime_ver10_99_19.02.2025.fasta'
        seq.write_text('>a\nACGT\n')
        taxo = developer / 'sh_taxonomy_qiime_ver10_99_19.02.2025.txt'
        taxo.write_text('a\tk__Fungi\n')

        found_seq, found_taxo = train_unite._find_unite_files(tmp_path)

        assert found_seq == seq
        assert found_taxo == taxo

    def test_does_not_confuse_the_singleton_inclusive_all_variant(self, tmp_path):
        """Regression test: real UNITE archives ship both
        sh_refs_qiime_ver10_99_<date>.fasta and the singleton-inclusive
        sh_refs_qiime_ver10_99_all_<date>.fasta side by side; matching by
        substring alone (both contain "sh_refs_qiime" and "_99_") picked
        whichever os.walk() happened to see last."""
        developer = tmp_path / 'sh_qiime_release_19.02.2025' / 'developer'
        developer.mkdir(parents=True)
        seq = developer / 'sh_refs_qiime_ver10_99_19.02.2025.fasta'
        seq.write_text('>a\nACGT\n')
        (developer / 'sh_refs_qiime_ver10_99_all_19.02.2025.fasta').write_text('>b\nTTTT\n')
        taxo = developer / 'sh_taxonomy_qiime_ver10_99_19.02.2025.txt'
        taxo.write_text('a\tk__Fungi\n')
        (developer / 'sh_taxonomy_qiime_ver10_99_all_19.02.2025.txt').write_text('a\tk__Fungi\nb\tk__Fungi\n')

        found_seq, found_taxo = train_unite._find_unite_files(tmp_path)

        assert found_seq == seq
        assert found_taxo == taxo

    def test_does_not_confuse_other_clustering_thresholds(self, tmp_path):
        developer = tmp_path / 'sh_qiime_release_19.02.2025' / 'developer'
        developer.mkdir(parents=True)
        seq = developer / 'sh_refs_qiime_ver10_99_19.02.2025.fasta'
        seq.write_text('>a\nACGT\n')
        (developer / 'sh_refs_qiime_ver10_97_19.02.2025.fasta').write_text('>b\nTTTT\n')
        (developer / 'sh_refs_qiime_ver10_dynamic_19.02.2025.fasta').write_text('>c\nGGGG\n')
        taxo = developer / 'sh_taxonomy_qiime_ver10_99_19.02.2025.txt'
        taxo.write_text('a\tk__Fungi\n')

        found_seq, _found_taxo = train_unite._find_unite_files(tmp_path)

        assert found_seq == seq

    def test_raises_when_files_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            train_unite._find_unite_files(tmp_path)


class TestFixFasta:
    def test_uppercases_sequence_lines_and_strips_spaces(self, tmp_path):
        input_fasta = tmp_path / 'in.fasta'
        input_fasta.write_text('>header with spaces\nacgt n n\n>h2\nTTTT\n')
        output_fasta = tmp_path / 'out.fasta'

        train_unite.fix_fasta(input_fasta, output_fasta)

        lines = output_fasta.read_text().splitlines()
        assert lines[0] == '>headerwithspaces'
        assert lines[1] == 'ACGTNN'
        assert lines[2] == '>h2'
        assert lines[3] == 'TTTT'


class TestRun:
    def test_full_run_orchestrates_expected_calls(self, tmp_path, mocker):
        output_folder = tmp_path / 'out'
        archive = tmp_path / 'unite.tgz'
        archive.write_bytes(b'fake')

        mocker.patch('qiime2_its.cli.train_unite.env_checks.check_qiime2_env_active',
                     return_value='rachis-qiime2-2026.7')

        def fake_extract(_archive, output_path):
            developer = Path(output_path) / 'developer'
            developer.mkdir(parents=True, exist_ok=True)
            (developer / 'sh_refs_qiime_ver10_99_19.02.2025.fasta').write_text('>a\nacgt\n')
            (developer / 'sh_taxonomy_qiime_ver10_99_19.02.2025.txt').write_text('a\tk__Fungi\n')

        mocker.patch('qiime2_its.cli.train_unite.downloader.extract_targz', side_effect=fake_extract)
        mock_import_seq = mocker.patch('qiime2_its.cli.train_unite.qiime_wrapper.import_sequences')
        mock_import_taxo = mocker.patch('qiime2_its.cli.train_unite.qiime_wrapper.import_taxonomy')
        mock_train = mocker.patch('qiime2_its.cli.train_unite.qiime_wrapper.train_naive_bayes_classifier')

        train_unite.run(str(archive), output_folder, 'rachis-qiime2-2026.7')

        mock_import_seq.assert_called_once()
        mock_import_taxo.assert_called_once()
        mock_train.assert_called_once()
        classifier_path = mock_train.call_args.args[2]
        assert str(classifier_path).endswith('unite-ver10-99-classifier-19.02.2025.qza')

    def test_warns_when_declared_env_does_not_match_the_active_one(self, tmp_path, mocker, capsys):
        """Regression test: -q/--qiime2 used to be parsed but never passed to
        run() at all, so it silently did nothing no matter what -- including
        actually being a different environment than the one really active."""
        output_folder = tmp_path / 'out'
        archive = tmp_path / 'unite.tgz'
        archive.write_bytes(b'fake')

        mocker.patch('qiime2_its.cli.train_unite.env_checks.check_qiime2_env_active',
                     return_value='rachis-qiime2-2026.7')

        def fake_extract(_archive, output_path):
            developer = Path(output_path) / 'developer'
            developer.mkdir(parents=True, exist_ok=True)
            (developer / 'sh_refs_qiime_ver10_99_19.02.2025.fasta').write_text('>a\nacgt\n')
            (developer / 'sh_taxonomy_qiime_ver10_99_19.02.2025.txt').write_text('a\tk__Fungi\n')

        mocker.patch('qiime2_its.cli.train_unite.downloader.extract_targz', side_effect=fake_extract)
        mocker.patch('qiime2_its.cli.train_unite.qiime_wrapper.import_sequences')
        mocker.patch('qiime2_its.cli.train_unite.qiime_wrapper.import_taxonomy')
        mocker.patch('qiime2_its.cli.train_unite.qiime_wrapper.train_naive_bayes_classifier')

        train_unite.run(str(archive), output_folder, 'some-other-env')

        assert 'some-other-env' in capsys.readouterr().out


class TestMain:
    def test_every_parsed_arg_reaches_run(self, mocker):
        """Regression test: -q/--qiime2 was parsed (and required) but main()
        called run(args.url, args.output_folder) -- silently dropping it, so
        it had zero effect no matter what the user passed. Nothing caught
        this because nothing tested the parse_args() -> main() -> run()
        wiring itself, only run() called directly."""
        mock_run = mocker.patch('qiime2_its.cli.train_unite.run')
        mocker.patch('sys.argv', ['qiime2-its-train-unite', '-u', 'unite.tgz',
                                   '-o', '/out', '-q', 'rachis-qiime2-2026.7'])

        train_unite.main()

        mock_run.assert_called_once_with('unite.tgz', '/out', 'rachis-qiime2-2026.7')
