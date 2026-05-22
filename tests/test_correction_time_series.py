"""Tests for gravtools.tides.correction_time_series."""

import numpy as np
import pytest

from gravtools.tides.correction_time_series import TimeSeries

# 2024-01-01 00:00:00 UTC = 1704067200 s
_REF_UNIX = 1704067200.0
_REF_DT = np.array(['2024-01-01T00:00:00'], dtype='datetime64[s]')

# 2024-01-01 01:00:00 UTC = 1704070800 s
_REF_DT_PLUS1H = np.array(['2024-01-01T01:00:00'], dtype='datetime64[s]')
_REF_UNIX_PLUS1H = _REF_UNIX + 3600.0


def _make_ts(ref_time_dt: np.ndarray, data: np.ndarray) -> TimeSeries:
    """Create a minimal TimeSeries for testing."""
    return TimeSeries(
        ref_time_dt=ref_time_dt,
        data=data,
        unit='µGal',
        data_source='test',
    )


class TestTimeSeriesRefTimeUnix:
    """Tests for TimeSeries.ref_time_unix property."""

    def test_single_timestamp(self):
        """Single known timestamp returns correct Unix value."""
        ts = _make_ts(_REF_DT, np.array([0.0]))
        result = ts.ref_time_unix
        assert result == pytest.approx([_REF_UNIX])

    def test_multiple_timestamps(self):
        """Multiple timestamps: length and values are correct."""
        ref = np.concatenate([_REF_DT, _REF_DT_PLUS1H])
        ts = _make_ts(ref, np.array([0.0, 1.0]))
        result = ts.ref_time_unix
        assert len(result) == 2
        assert result == pytest.approx([_REF_UNIX, _REF_UNIX_PLUS1H])

    def test_returns_float64(self):
        """Return dtype is float64."""
        ts = _make_ts(_REF_DT, np.array([0.0]))
        result = ts.ref_time_unix
        assert result.dtype == np.float64
