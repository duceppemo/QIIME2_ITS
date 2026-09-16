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
    qiime_wrapper.dada2_denoise_single('r.qza', 'rep.qza', 't.qza', 's.qza', 'bts.qza')
    cmd = mock_run.call_args.args[0]
    assert '--p-trim-left' in cmd and cmd[cmd.index('--p-trim-left') + 1] == '0'
    assert '--p-trunc-len' in cmd and cmd[cmd.index('--p-trunc-len') + 1] == '0'
    assert '--o-base-transition-stats' in cmd and cmd[cmd.index('--o-base-transition-stats') + 1] == 'bts.qza'


def test_dada2_denoise_paired_no_trimming(mock_run):
    qiime_wrapper.dada2_denoise_paired('r.qza', 'rep.qza', 't.qza', 's.qza', 'bts.qza')
    cmd = mock_run.call_args.args[0]
    for flag in ('--p-trim-left-f', '--p-trim-left-r', '--p-trunc-len-f', '--p-trunc-len-r'):
        assert flag in cmd and cmd[cmd.index(flag) + 1] == '0'
    assert '--o-base-transition-stats' in cmd and cmd[cmd.index('--o-base-transition-stats') + 1] == 'bts.qza'


def test_dada2_denoise_single_defaults_match_qiime2(mock_run):
    qiime_wrapper.dada2_denoise_single('r.qza', 'rep.qza', 't.qza', 's.qza', 'bts.qza')
    cmd = mock_run.call_args.args[0]
    assert cmd[cmd.index('--p-max-ee') + 1] == '2.0'
    assert cmd[cmd.index('--p-trunc-q') + 1] == '2'
    assert cmd[cmd.index('--p-pooling-method') + 1] == 'independent'
    assert cmd[cmd.index('--p-chimera-method') + 1] == 'consensus'
    assert '--p-no-allow-one-off' in cmd


def test_dada2_denoise_single_ion_torrent_style_overrides(mock_run):
    qiime_wrapper.dada2_denoise_single('r.qza', 'rep.qza', 't.qza', 's.qza', 'bts.qza',
                                        max_ee=5.0, allow_one_off=True)
    cmd = mock_run.call_args.args[0]
    assert cmd[cmd.index('--p-max-ee') + 1] == '5.0'
    assert '--p-allow-one-off' in cmd
    assert '--p-no-allow-one-off' not in cmd


def test_dada2_denoise_paired_uses_separate_max_ee_per_direction(mock_run):
    qiime_wrapper.dada2_denoise_paired('r.qza', 'rep.qza', 't.qza', 's.qza', 'bts.qza', max_ee_f=1.0, max_ee_r=3.0)
    cmd = mock_run.call_args.args[0]
    assert cmd[cmd.index('--p-max-ee-f') + 1] == '1.0'
    assert cmd[cmd.index('--p-max-ee-r') + 1] == '3.0'


def test_core_diversity_sampling_depth_override(mock_run):
    qiime_wrapper.core_diversity(4, 'metadata.tsv', 'rooted.qza', 'table.qza', 'out', sampling_depth=50)
    cmd = mock_run.call_args.args[0]
    assert cmd[cmd.index('--p-sampling-depth') + 1] == '50'


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


def test_sample_summarize(mock_run):
    qiime_wrapper.sample_summarize('meta.tsv', 'table.qza', 'table.qzv', 'ff.qza', 'sf.qza')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'feature-table', 'summarize',
                    '--m-metadata-file', 'meta.tsv', '--i-table', 'table.qza',
                    '--o-summary', 'table.qzv', '--o-feature-frequencies', 'ff.qza',
                    '--o-sample-frequencies', 'sf.qza']


def test_classify(mock_run):
    qiime_wrapper.classify('clf.qza', 'rep.qza', 'taxo.qza')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'feature-classifier', 'classify-sklearn',
                    '--p-n-jobs', '0', '--i-classifier', 'clf.qza',
                    '--i-reads', 'rep.qza', '--o-classification', 'taxo.qza']


def test_alpha_group_significance(mock_run):
    qiime_wrapper.alpha_group_significance('shannon.qza', 'meta.tsv', 'out.qzv')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'diversity', 'alpha-group-significance',
                    '--i-alpha-diversity', 'shannon.qza',
                    '--m-metadata-file', 'meta.tsv',
                    '--o-visualization', 'out.qzv']


def test_beta_group_significance_default_method(mock_run):
    qiime_wrapper.beta_group_significance('bray.qza', 'meta.tsv', 'site', 'out.qzv')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'diversity', 'beta-group-significance',
                    '--i-distance-matrix', 'bray.qza',
                    '--m-metadata-file', 'meta.tsv',
                    '--m-metadata-column', 'site',
                    '--p-method', 'permanova',
                    '--o-visualization', 'out.qzv']


def test_beta_group_significance_custom_method(mock_run):
    qiime_wrapper.beta_group_significance('bray.qza', 'meta.tsv', 'site', 'out.qzv', method='anosim')
    cmd = mock_run.call_args.args[0]
    assert cmd[cmd.index('--p-method') + 1] == 'anosim'


def test_taxa_collapse(mock_run):
    qiime_wrapper.taxa_collapse('table.qza', 'taxo.qza', 6, 'collapsed.qza')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'taxa', 'collapse',
                    '--i-table', 'table.qza', '--i-taxonomy', 'taxo.qza',
                    '--p-level', '6', '--o-collapsed-table', 'collapsed.qza']


def test_relative_frequency(mock_run):
    qiime_wrapper.relative_frequency('table.qza', 'rel.qza')
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'feature-table', 'relative-frequency',
                    '--i-table', 'table.qza', '--o-relative-frequency-table', 'rel.qza']


def test_classify_samples(mock_run):
    qiime_wrapper.classify_samples('table.qza', 'meta.tsv', 'site', 'out_dir', cv=3)
    cmd = mock_run.call_args.args[0]
    assert cmd == ['qiime', 'sample-classifier', 'classify-samples',
                    '--i-table', 'table.qza', '--m-metadata-file', 'meta.tsv',
                    '--m-metadata-column', 'site', '--p-cv', '3',
                    '--p-n-estimators', '100', '--output-dir', 'out_dir']


def test_classify_samples_custom_n_estimators(mock_run):
    qiime_wrapper.classify_samples('table.qza', 'meta.tsv', 'site', 'out_dir', cv=2, n_estimators=50)
    cmd = mock_run.call_args.args[0]
    assert cmd[cmd.index('--p-n-estimators') + 1] == '50'
