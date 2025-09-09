import pytest


def test_something():
    with pytest.raises(ZeroDivisionError):
        _ = 1 / 0
