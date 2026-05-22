"""Tests for gravtools.models.misc utility functions."""

import numpy as np
import pandas as pd
import datetime as dt
import pytest

from gravtools.models.misc import to_unix_seconds

# 2024-01-01 00:00:00 UTC
_REF_DT = dt.datetime(2024, 1, 1, tzinfo=dt.timezone.utc)
_REF_UNIX = 1704067200.0

# 2024-01-01 01:00:00 UTC  (+3600 s)
_REF_DT_PLUS1H = dt.datetime(2024, 1, 1, 1, 0, 0, tzinfo=dt.timezone.utc)
_REF_UNIX_PLUS1H = _REF_UNIX + 3600.0


class TestToUnixSeconds:
    """Tests for to_unix_seconds()."""

    def test_series_ns(self):
        """pd.Series with datetime64[ns, UTC] (old pandas default)."""
        s = pd.Series([_REF_DT]).dt.tz_localize(None).dt.tz_localize('UTC').astype('datetime64[ns, UTC]')
        result = to_unix_seconds(s)
        assert result == pytest.approx([_REF_UNIX])

    def test_series_us(self):
        """pd.Series with datetime64[us, UTC] (pandas 3.0 default)."""
        s = pd.to_datetime([_REF_DT]).tz_convert('UTC').astype('datetime64[us, UTC]')
        s = pd.Series(s)
        result = to_unix_seconds(s)
        assert result == pytest.approx([_REF_UNIX])

    def test_ndarray_datetime64s(self):
        """np.ndarray with datetime64[s]."""
        arr = np.array(['2024-01-01T00:00:00'], dtype='datetime64[s]')
        result = to_unix_seconds(arr)
        assert result == pytest.approx([_REF_UNIX])

    def test_returns_float64(self):
        """Return dtype is always float64."""
        s = pd.Series([_REF_DT]).dt.tz_localize(None).dt.tz_localize('UTC')
        result = to_unix_seconds(s)
        assert result.dtype == np.float64

    def test_ndarray_returns_float64(self):
        """Return dtype is float64 also for ndarray input."""
        arr = np.array(['2024-01-01T00:00:00'], dtype='datetime64[s]')
        result = to_unix_seconds(arr)
        assert result.dtype == np.float64

    def test_multiple_values(self):
        """Multiple timestamps: order and values are correct."""
        times = [_REF_DT, _REF_DT_PLUS1H]
        s = pd.Series(times).dt.tz_localize(None).dt.tz_localize('UTC')
        result = to_unix_seconds(s)
        assert len(result) == 2
        assert result == pytest.approx([_REF_UNIX, _REF_UNIX_PLUS1H])

    def test_series_nonzero_offset(self):
        """Non-epoch timestamp: difference to epoch is correct."""
        s = pd.Series([_REF_DT_PLUS1H]).dt.tz_localize(None).dt.tz_localize('UTC')
        result = to_unix_seconds(s)
        assert result == pytest.approx([_REF_UNIX_PLUS1H])
        assert result[0] - _REF_UNIX == pytest.approx(3600.0)
