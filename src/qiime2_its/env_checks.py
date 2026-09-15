"""Sanity checks shared by all CLIs: CPU/parallelism clamping and conda env validation."""
import os
from multiprocessing import cpu_count


def clamp_cpu(requested_cpu, available_cpu=None):
    """Return a valid CPU count: `available_cpu` if `requested_cpu` is out of [1, available_cpu]."""
    available_cpu = available_cpu if available_cpu is not None else cpu_count()
    if requested_cpu < 1 or requested_cpu > available_cpu:
        return available_cpu
    return requested_cpu


def clamp_parallel(requested_parallel, cpu):
    """Never run more parallel processes than there are CPUs (1 CPU per process minimum)."""
    return min(requested_parallel, cpu)


def check_qiime2_env_active(env_var='CONDA_DEFAULT_ENV'):
    """Raise if the active conda environment doesn't look like a QIIME2 environment.

    Matches both the legacy ``qiime2-*`` naming and the current ``rachis-qiime2-*`` naming.
    Returns the active environment name.
    """
    env_name = os.environ.get(env_var, '')
    if 'qiime2' not in env_name:
        raise EnvironmentError(
            'You must activate your QIIME2 conda environment to run this script, '
            'e.g. "conda activate rachis-qiime2-2026.7".'
        )
    return env_name
