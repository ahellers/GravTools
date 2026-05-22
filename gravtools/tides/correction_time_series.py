"""Class for correction data provided as time series per station and survey.

Copyright (C) 2023  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from dataclasses import dataclass
import typing
from datetime import datetime, timezone
import warnings

import pandas as pd
import scipy
import numpy as np

import gravtools.tides.correction_time_series
from gravtools.tides.tide_data_tfs import TSF
from gravtools import settings
from gravtools.models.misc import to_unix_seconds

class CorrectionTimeSeries:
    """Class for correction data provided as time series per station and survey.

    Attributes
    ----------
    surveys: dict of `gravtools.tides.correction_time_series.SurveyCorrections`
        Dict of `SurveyCorrection` objects with ths survey names as keys.
    """

    def __init__(self):
        """Default initializer."""
        self.surveys = {}

    def load_tfs_file(self, survey_name: str, filename_tsf: str, location: str ='', instrument: str ='',
                      data_type: str ='', overwrite_channel: bool=True, is_correction: bool='False'):
        """Load time series data from a TSoft TSF file for one survey.

        Parameters
        ----------
        survey_name: str
            Name of the survey for which the data is loaded.
        filename_tsf: str
            Name and path of the TSF file from which the correction time series data is loaded.
        location: str, optional (default='')
            If not empty, only channels with matching locations are loaded.
        instrument: str, optional (default='')
            If not empty, only channels with matching instruments are loaded.
        data_type: str, optional (default='')
            If not empty, only channels with matching data types are loaded.
        overwrite_channel: bool, optional (default=`True`)
            `True` implies that an existing StationCorrection object of a station will be overwritten.
        is_correction: bool, optional (default=`False`)
            `True` implies that the channel data loaded from the TSF file model the gravity effect of phenomena, rather
            than corrections for gravity observations. Both options have the opposite sign: while corrections have to be
            added to gravity observations, effects have to be subtracted, in order to reduce observations.
            See TSoft manual (version 2.2.4 Release date 2015-09-09), p. 15: "The loading calculation by Tsoft will
            result in the correction, not the effect.
            On the other hand the prediction of the solid Earth tides using the WDD parameter set
            provided by Tsoft directly will result in the effect. Therefore both must be treated with
            different sign in order to reduce a gravity time series correctly."
        """
        # Load TSF file:
        tsf_data = TSF.from_tfs_file(filename_tsf)

        # Loop over channel and retrieve required data:
        retrieved_channel_locations = []
        retrieved_channel_numbers = []
        stations_corrections_dict = {}
        for channel_metadata in tsf_data.channel_metadata:
            ch_channel_name = channel_metadata.channel_name
            ch_instrument = channel_metadata.instrument_name
            ch_data_type = channel_metadata.data_type
            ch_location = channel_metadata.location
            channel_number = channel_metadata.channel_number

            if location:
                if location != ch_location:
                    continue
            if instrument:
                if instrument != ch_instrument:
                    continue
            if data_type:
                if data_type != ch_data_type:
                    continue
            if ch_location in retrieved_channel_locations:
                raise RuntimeError(f'Ambiguous TSF channel matching! The channels '
                                   f'{retrieved_channel_numbers[retrieved_channel_locations.index(ch_location)]} and '
                                   f'{channel_number} share the same location ({ch_location}). '
                                   f'Please constrain the channel selection by providing additional arguments ('
                                   f'instrument, data type) or edit the TSF file.')
            retrieved_channel_locations.append(ch_location)
            retrieved_channel_numbers.append(channel_number)
            ch_epoch_dt, ch_data  = tsf_data.get_channel_np(channel=channel_number)
            time_series = TimeSeries(ref_time_dt=ch_epoch_dt, data=ch_data, unit=channel_metadata.unit,
                                     data_source=f'{tsf_data.filename_without_path} ({tsf_data.filetype})',
                                     description=ch_channel_name,
                                     is_correction=is_correction)


            stations_corrections_dict[ch_location] = StationCorrections(station_name=ch_location,
                                                                        tidal_correction=time_series)

        # Checks:
        if len(retrieved_channel_locations) == 0:
            tmp_str_list = []
            if location:
                tmp_str_list.append(f'station "{location}"')
            if instrument:
                tmp_str_list.append(f'instrument "{instrument}"')
            if data_type:
                tmp_str_list.append(f'data_type "{data_type}"')
            if tmp_str_list:
                raise RuntimeError(f'Could not retrieve any data from TSF file {filename_tsf} with the following '
                                   f'restrictions:' + ', '.join(tmp_str_list) + '.')
            else:
                raise RuntimeError(f'Could not retrieve any data from TSF file {filename_tsf}.')

        # Add loaded data:
        if survey_name in self.surveys.keys():
            # Add station corrections to existing survey correction object:
            for stat_name, stat_corr in stations_corrections_dict.items():
                _ = self.surveys[survey_name].add_station_correction(station_name=stat_name,
                                                                     station_correction=stat_corr,
                                                                     overwrite=overwrite_channel)
        else:
            # Create survey correction object and add the station correction data:
            self.surveys[survey_name] = SurveyCorrections(survey_name=survey_name, stations=stations_corrections_dict)

    def delete_survey_correction(self, survey_name):
        """Removes a survey correction object.

        Parameters
        ----------
        survey_name: str
            Survey name.
        """
        del self.surveys[survey_name]

@dataclass
class TimeSeries:
    """Uniformly or irregularly sampled correction or effect time series for a single channel.

    This dataclass stores a paired array of reference epochs and scalar values together with
    metadata describing the source and physical meaning of the data.  It is the lowest-level
    building block used by :class:`StationCorrections` and :class:`SurveyCorrections` to hold,
    e.g., tidal corrections or atmospheric loading effects that are later interpolated onto
    individual gravity observations.

    The sign convention follows the TSoft convention (see TSoft manual v2.2.4, p. 15):

    * **Correction** (``is_correction=True``): the value is *added* to the raw gravity
      observation (e.g. the tidal loading correction as exported by TSoft).
    * **Effect** (``is_correction=False``): the value represents the *gravitational effect*
      of a phenomenon and must be *subtracted* from the observation to reduce it.

    Attributes
    ----------
    ref_time_dt : numpy.ndarray of datetime64
        Reference epochs of the time series, in UTC.
    data : numpy.ndarray of float
        Scalar values associated with each epoch.  Must have the same length as
        ``ref_time_dt``.
    unit : str
        Physical unit of ``data`` (e.g. ``'nm/s^2'``, ``'µGal'``).  Used for unit
        conversion via :data:`gravtools.settings.UNIT_CONVERSION_TO_MUGAL`.
    data_source : str
        Human-readable identification of the data origin, e.g. the file name and
        file type from which the series was loaded.
    description : str, optional
        Short free-text label for the time series (e.g. the TSF channel name).
        Defaults to ``''``.
    is_correction : bool, optional
        ``True`` if ``data`` contains *corrections* (to be added to observations);
        ``False`` if ``data`` contains *effects* (to be subtracted).
        Defaults to ``True``.
    created_utc_dt : datetime
        UTC timestamp set automatically in :meth:`__post_init__` when the object is
        created.  Not a constructor argument.
    """
    ref_time_dt: np.array  # np.array of DateTime objects
    data: np.array
    unit: str
    data_source: str
    description: str = ''
    is_correction: bool = True
    # created_utc_dt: datetime = None  # Is initialized in __post_init__

    def __post_init__(self):
        """Executed at the very end of __init__"""
        self.created_utc_dt = datetime.now(tz=timezone.utc)

    @property
    def created_datetime_utc(self):
        """Returns the time and date (in UTC) when the times series was created as datetime object."""
        return self.created_utc_dt

    @property
    def created_datetime_utc_str(self):
        """Returns the time and date (in UTC) when the times series was created as string."""
        return self.created_utc_dt.strftime('%Y-%m-%d %H:%M:%S.%f')

    @property
    def start_datetime(self):
        """Returns the start date and time."""
        return min(self.ref_time_dt)

    @property
    def start_datetime_str(self):
        """Returns the start date and time as string."""
        return pd.to_datetime(self.start_datetime).strftime('%Y-%m-%d %H:%M:%S.%f')

    @property
    def end_datetime(self):
        """Returns the end date and time."""
        return max(self.ref_time_dt)

    @property
    def end_datetime_str(self):
        """Returns the end date and time as string."""
        return pd.to_datetime(self.end_datetime).strftime('%Y-%m-%d %H:%M:%S.%f')

    @property
    def duration_timedelta64_ns(self):
        """Returns the duration as `numpy.timedelta64[ns]` object."""
        return max(self.ref_time_dt) - min(self.ref_time_dt)

    def duration_dhms(self):
        """Returns the duration as days, hours, minutes and seconds."""
        dt = max(self.ref_time_dt) - min(self.ref_time_dt)
        days = dt.astype('timedelta64[D]').astype(int)
        hours = dt.astype('timedelta64[h]').astype(int) - days * 24
        minutes = dt.astype('timedelta64[m]').astype(int) - days*24*60 - hours*60
        seconds = dt.astype('timedelta64[s]').astype(float) - days*24*60*60 - hours*60*60 - minutes*60
        return days, hours, minutes, seconds

    @property
    def duration_dhms_str(self):
        """Returns the duration as string (days, hours, minutes and seconds)"""
        days, hours, minutes, seconds = self.duration_dhms()
        return f'{days} days, {hours} hours, {minutes} minutes, {seconds} seconds'

    @property
    def number_of_datapoints(self):
        """Returns the number of datapoints."""
        return len(self.data)

    @property
    def model_type(self):
        """Returns the model type, i.e. whether effects or corrections are modeled."""
        if self.is_correction:
            return 'Correction'
        else:
            return 'Effect'

    def interpolate(self, interp_times: [np.array, typing.List[datetime]], kind: str = 'quadratic',
                    return_correction: bool = False) -> np.ndarray :
        """Return an interpolated value for the given interpolation epoch.

        Parameters
        ----------
        interp_times: list(datetime) or np.array(datetime)
            Interpolation epoch given as `datetime` object w.r.t. UTC! The interpolation times have to bin within the
            time range of the data series.
        kind: str or int, optional (default='quadratic')
            Specifies, which interpolation method is used. This argument is directly passed to
            `scipy.interpolate.inter1`. For options see scipy reference.
        return_correction: bool, optional (default=`False`)
            If `True`, corrections are returned considering the `self.is_correction` attribute.

        Returns
        -------
        `numpy.ndarray`: Interpolated value for the given interpolation times.
        """
        # Convert datetimes to UNIX timestamps, because numerical values are needed for interpolation:
        interp_times_unix = np.fromiter((t.replace(tzinfo=timezone.utc).timestamp() for t in interp_times), float)
        x = self.ref_time_unix
        y = self.data
        interp_func = scipy.interpolate.interp1d(x, y, kind=kind)
        interp_values = interp_func(interp_times_unix)
        if return_correction:
            if not self.is_correction:
                interp_values = interp_values * -1  # Convert from effect to correction!
        return interp_values

    def to_df(self):
        """Returns the time series as pandas DataFrame with the time reference as index (sorted!)."""
        df = pd.DataFrame({'epoch_dt': self.ref_time_dt,'data': self.data})
        df.set_index('epoch_dt', inplace=True)
        return df

    @property
    def ref_time_unix(self):
        """Returns the reference times as UNIX timestamps (seconds since Jan 1, 1970)."""
        return to_unix_seconds(self.ref_time_dt)


@dataclass
class StationCorrections:
    """Correction time series data for one station.

    Notes
    -----
    For each station and correction type (e.g. tidal correction) only ONE time series can be added!
    """
    station_name: str
    tidal_correction: TimeSeries


@dataclass
class SurveyCorrections:
    """Correction time series data for one survey with multiple stations."""
    survey_name: str
    stations: typing.Dict[str, StationCorrections]

    def add_station_correction(self, station_name: str, station_correction: StationCorrections, overwrite: bool=True,
                               warn: bool=True):
        """Add a StationCorrection object to `self.stations`.

        Parameters
        ----------
        station_name: str
            Station name.
        station_correction: `StationCorrections`
            StationCorrections object containing correction time series data for a specific station.
        overwrite: bool, optional (default=`True`)
            `True` implies that an existing StationCorrection object of a station will be overwritten.
        warn: bool, optional (default=`True`)
            `True` implies that a warning is issued if an existing StationCorrection will or would bne be overwritten.

        Returns
        -------
        bool: `True`, if the station correction data was added. Otherwise, `False`.
        """
        # Check, if data for this station already exists:
        if station_name in self.stations.keys():
            if not overwrite:
                if warn:
                    warnings.warn(f'Station correction data for station {station_name} already exists! Overwriting '
                                  f'data is not permitted!')
                return False
            else:
                if warn:
                    warnings.warn(f'Station correction data for station {station_name} already exists! Overwriting '
                                  f'data is permitted and the existing data will be overwritten!')
                self.stations[station_name] = station_correction
                return True
        else:
            self.stations[station_name] = station_correction
            return True

    def delete_station_correction(self, station_name):
        """Removes a station correction object.

        Parameters
        ----------
        station_name: str
            Station name.
        """
        del self.stations[station_name]

    @property
    def number_of_stations(self):
        """Returns the number of station correction objects."""
        return len(self.stations)

    @property
    def station_names(self):
        """Returns the names of all stations."""
        return list(self.stations.keys())

@staticmethod
def convert_to_mugal(data, data_unit: str):
    """Converts data to from a given unit to µGal."""
    if data_unit not in settings.UNIT_CONVERSION_TO_MUGAL:
        raise RuntimeError(f'Conversion factor from {data_unit} to µGal not defined. Please add the factor to'
                           f'gravtools.settings.UNIT_CONVERSION_TO_MUGAL.')
    return data * settings.UNIT_CONVERSION_TO_MUGAL[data_unit]








