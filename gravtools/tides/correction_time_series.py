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
                      data_type: str =''):
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
                                     description=ch_channel_name)


            stations_corrections_dict[ch_location] = StationCorrections(station_name=ch_location, tidal_correction=time_series)

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
                                                                     station_correction=stat_corr)
        else:
            # Create survey correction object and add the station correction data:
            self.surveys[survey_name] = SurveyCorrections(survey_name=survey_name, stations=stations_corrections_dict)


@dataclass
class TimeSeries:
    """Time series data."""
    ref_time_dt: np.array  # np.array of DateTime objects
    data: np.array
    unit: str
    data_source: str
    description: str = ''

    def interpolate(self, interp_times: [np.array, typing.List[datetime]], kind: str = 'quadratic') -> np.ndarray :
        """Return an interpolated value for the given interpolation epoch.

        Parameters
        ----------
        interp_times: list(datetime) or np.array(datetime)
            Interpolation epoch given as `datetime` object w.r.t. UTC! The interpolation times have to bin within the
            time range of the data series.
        kind: str or int, optional (default='quadratic')
            Specifies, which interpolation method is used. This argument is directly passed to
            `scipy.interpolate.inter1`. For options see scipy reference.

        Returns
        -------
        `numpy.ndarray`: Interpolated value for the given interpolation times.
        """
        # Convert datetimes to UNIX timestamps, because numerical values are needed for interpolation:
        interp_times_unix = np.fromiter((t.replace(tzinfo=timezone.utc).timestamp() for t in interp_times), float)
        x = self.ref_time_unix
        y = self.data
        interp_func = scipy.interpolate.interp1d(x, y, kind=kind)
        return interp_func(interp_times_unix)

    def to_df(self):
        """Returns the time series as pandas DataFrame with the time reference as index (sorted!)."""
        df = pd.DataFrame({'epoch_dt': self.ref_time_dt,'data': self.data})
        df.set_index('epoch_dt', inplace=True)
        return df

    @property
    def ref_time_unix(self):
        """Returns the reference times as UNIX timestamps (seconds since Jan 1, 1970)."""
        return self.ref_time_dt.astype('int64')/1e9

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
        bool: `True`, if the the station correction data was added. Otherwise, `False`.
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
                self.stations['station_name'] = station_correction
                return True
        else:
            self.stations['station_name'] = station_correction
            return True






