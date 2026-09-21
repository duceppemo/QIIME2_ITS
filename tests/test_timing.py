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


def test_fractional_seconds_never_round_up_to_a_full_unit():
    """Regression test: each unit was rounded after the split, so 119.7 s
    came out as "1m60s"."""
    assert timing.format_elapsed(119.7) == '2m'
    assert timing.format_elapsed(59.6) == '1m'
