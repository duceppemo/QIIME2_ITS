import pytest

from qiime2_its import env_checks


class TestClampCpu:
    @pytest.mark.parametrize('requested,available,expected', [
        (4, 8, 4),      # within range, unchanged
        (8, 8, 8),      # exactly at the limit
        (0, 8, 8),      # below 1 -> clamp to available
        (-1, 8, 8),     # negative -> clamp to available
        (100, 8, 8),    # above available -> clamp to available
        (1, 8, 1),      # lower boundary
    ])
    def test_clamp(self, requested, available, expected):
        assert env_checks.clamp_cpu(requested, available_cpu=available) == expected


class TestClampParallel:
    def test_parallel_capped_to_cpu(self):
        assert env_checks.clamp_parallel(10, cpu=4) == 4

    def test_parallel_within_cpu_unchanged(self):
        assert env_checks.clamp_parallel(2, cpu=4) == 2

    def test_zero_clamped_to_one(self):
        """Regression test: 0 used to pass straight through (min(0, cpu) ==
        0), which downstream code divides by (ZeroDivisionError in
        size_filter.py) or passes to ThreadPoolExecutor(max_workers=0)
        (ValueError) -- both an opaque crash instead of a clear message."""
        assert env_checks.clamp_parallel(0, cpu=4) == 1

    def test_negative_clamped_to_one(self):
        assert env_checks.clamp_parallel(-3, cpu=4) == 1


class TestCheckQiime2EnvActive:
    def test_raises_when_no_env_active(self, monkeypatch):
        monkeypatch.delenv('CONDA_DEFAULT_ENV', raising=False)
        with pytest.raises(EnvironmentError):
            env_checks.check_qiime2_env_active()

    def test_raises_when_unrelated_env_active(self, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'base')
        with pytest.raises(EnvironmentError):
            env_checks.check_qiime2_env_active()

    def test_accepts_legacy_naming(self, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'qiime2-2022.8')
        assert env_checks.check_qiime2_env_active() == 'qiime2-2022.8'

    def test_accepts_rachis_naming(self, monkeypatch):
        monkeypatch.setenv('CONDA_DEFAULT_ENV', 'rachis-qiime2-2026.7')
        assert env_checks.check_qiime2_env_active() == 'rachis-qiime2-2026.7'


class TestCheckExecutable:
    def test_raises_with_hint_when_missing(self, mocker):
        mocker.patch('qiime2_its.env_checks.shutil.which', return_value=None)
        with pytest.raises(EnvironmentError, match='bbduk.sh.*Install BBTools'):
            env_checks.check_executable('bbduk.sh', 'Install BBTools/BBMap.')

    def test_raises_without_hint_when_missing(self, mocker):
        mocker.patch('qiime2_its.env_checks.shutil.which', return_value=None)
        with pytest.raises(EnvironmentError):
            env_checks.check_executable('bbduk.sh')

    def test_does_not_raise_when_present(self, mocker):
        mocker.patch('qiime2_its.env_checks.shutil.which', return_value='/usr/bin/bbduk.sh')
        env_checks.check_executable('bbduk.sh')  # no raise
