import pytest

from qiime2_its import timing


@pytest.mark.parametrize('seconds,expected', [
    (0, '0s'),
    (45, '45s'),
    (90, '1m30s'),
    (3661, '1h1m1s'),
    (90000, '1d1h'),
])
def test_format_elapsed(seconds, expected):
    assert timing.format_elapsed(seconds) == expected
