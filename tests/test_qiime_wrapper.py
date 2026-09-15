import pytest

from qiime2_its import qiime_wrapper


@pytest.fixture
def mock_run(mocker):
    return mocker.patch('qiime2_its.qiime_wrapper.subprocess.run')


def test_import_sequences(mock_run):
    qiime_wrapper.import_sequences('seq.fasta', 'seq.qza')
    mock_run.assert_called_once_with(
        ['qiime', 'tools', 'import', '--type', 'FeatureData[Sequence]',
         '--input-path', 'seq.fasta', '--output-path', 'seq.qza'],
        check=True)


def test_import_taxonomy(mock_run):
    qiime_wrapper.import_taxonomy('taxo.txt', 'taxo.qza')
    mock_run.assert_called_once_with(
        ['qiime', 'tools', 'import', '--type', 'FeatureData[Taxonomy]',
         '--input-format', 'HeaderlessTSVTaxonomyFormat',
         '--input-path', 'taxo.txt', '--output-path', 'taxo.qza'],
        check=True)


def test_train_naive_bayes_classifier_default(mock_run):
    qiime_wrapper.train_naive_bayes_classifier('seq.qza', 'taxo.qza', 'clf.qza')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'feature-classifier', 'fit-classifier-naive-bayes',
                    '--i-reference-reads', 'seq.qza', '--i-reference-taxonomy', 'taxo.qza',
                    '--o-classifier', 'clf.qza']


def test_train_naive_bayes_classifier_verbose(mock_run):
    qiime_wrapper.train_naive_bayes_classifier('seq.qza', 'taxo.qza', 'clf.qza', verbose=True)
    cmd = mock_run.call_args.args[0]
    assert cmd[-1] == '--verbose'


def test_import_fastq_se(mock_run):
    qiime_wrapper.import_fastq_se('reads/', 'reads.qza')
    cmd = mock_run.call_args.args[0]
    assert 'SampleData[SequencesWithQuality]' in cmd
    assert 'CasavaOneEightSingleLanePerSampleDirFmt' in cmd


def test_import_fastq_pe(mock_run):
    qiime_wrapper.import_fastq_pe('reads/', 'reads.qza')
    cmd = mock_run.call_args.args[0]
    assert 'SampleData[PairedEndSequencesWithQuality]' in cmd


def test_dada2_denoise_single_no_trimming(mock_run):
    qiime_wrapper.dada2_denoise_single('r.qza', 'rep.qza', 't.qza', 's.qza')
    cmd = mock_run.call_args.args[0]
    assert '--p-trim-left' in cmd and cmd[cmd.index('--p-trim-left') + 1] == '0'
    assert '--p-trunc-len' in cmd and cmd[cmd.index('--p-trunc-len') + 1] == '0'


def test_dada2_denoise_paired_no_trimming(mock_run):
    qiime_wrapper.dada2_denoise_paired('r.qza', 'rep.qza', 't.qza', 's.qza')
    cmd = mock_run.call_args.args[0]
    for flag in ('--p-trim-left-f', '--p-trim-left-r', '--p-trunc-len-f', '--p-trunc-len-r'):
        assert flag in cmd and cmd[cmd.index(flag) + 1] == '0'


def test_core_diversity_actually_runs(mock_run):
    """Regression test: the original code built this command but never called
    subprocess.run(), so the core-diversity step silently no-op'd."""
    qiime_wrapper.core_diversity(4, 'metadata.tsv', 'rooted.qza', 'table.qza', 'out')
    mock_run.assert_called_once()
    cmd = mock_run.call_args.args[0]
    assert cmd[:3] == ['qiime', 'diversity', 'core-metrics-phylogenetic']
    assert '--output-dir' in cmd and cmd[cmd.index('--output-dir') + 1] == 'out/core-metrics-results'


def test_all_wrapper_calls_use_check_true(mock_run):
    qiime_wrapper.export('x.qza', 'out/')
    assert mock_run.call_args.kwargs.get('check') is True


def test_classify(mock_run):
    qiime_wrapper.classify('clf.qza', 'rep.qza', 'taxo.qza')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'feature-classifier', 'classify-sklearn',
                    '--p-n-jobs', '-1', '--i-classifier', 'clf.qza',
                    '--i-reads', 'rep.qza', '--o-classification', 'taxo.qza']
