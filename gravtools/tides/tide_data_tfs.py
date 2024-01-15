"""Class for TSF-formatted tidal data.

TSF Files are created e.g. by TSoft

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

import re
import datetime as dt
import io
from dataclasses import dataclass
import os

import pandas as pd
from typing import Dict, Literal

from gravtools.tides.tide_data import AbstractTideData


@dataclass
class ChannelMetadata:
    """Class for channel metadata in TSF files."""
    location: str
    instrument_name: str
    data_type: str
    unit: str
    channel_number: int

    @property
    def channel_name(self):
        """Returns the channels name as read from the TSF file"""
        return self.location + ':' + self.instrument_name + ':' + self.data_type


class TSF(AbstractTideData):
    """Class for tidal data loaded from TSF files.

    TSF files are the native format for the TSoft which allows to calculate time-series of synthetic tide for a
    specific station defined by geographical coordinates and height.

    TSoft:
    Download TSoft and its manual: http://seismologie.oma.be/en/downloads/tsoft
    Reference: Michel Van Camp, Paul Vauterin: Tsoft: graphical and interactive software for the analysis of time series
    and Earth tides. Comput. Geosci. 31(5): 631-640 (2005)
    DOI: https://doi.org/10.1016/j.cageo.2004.11.015

    Attributes
    ----------
    _filename : str
        Path and name of the TSF file.
    tfs_format : str
        Format version of the TSF file (from block: [TSF-file]). Only version v01.0 is supported!
    timeformat : str
        Format of the time reference (from block: [TIMEFORMAT]). Has to be "DATETIME" or "DATETIMEFRAC"!
    undetval : float
        Determines which number is used to denote undetermined or unknown values in the [DATA] block (from block:
        [UNDETVAL]).
    increment : float
        Determines the time delay (in seconds) between two consecutive data points (from block: [INCREMENT]).
    channel_metadata : list of ChannelMetadata objects
        List of ChannelMetadata objects with one item per channel. Contains metadata such as location, instrument name,
        data type und unit.
    data_df : pandas.DataFrame
        Contains all columns from the [DATA] block. Additionally, the column `epoch_dt` describes the reference time and
        date in as DateTime object.
    countinfo : int, optiona (default=-1)
        Specifies the number of data points in the file (from block: [COUNTINFO]). If this block is missing the in the
        TSF file the default value -1 is set.
    comment : str, optional (default='')
        Optional multiline comment (from block: [COMMENT]).
    """
    def __init__(self, filename: str, tfs_format: str, timeformat: str, undetval: float, increment: float,
                 channel_metadata: list, data_df: pd.DataFrame, countinfo: int = -1, comment: str = ''):
        """Default constructor."""
        self._filename = filename
        self.tfs_format = tfs_format
        self.timeformat = timeformat
        self.undetval = undetval
        self.increment = increment
        self.countinfo = countinfo
        self.comment = comment
        self.channel_metadata = channel_metadata
        self.data_df = data_df


    @classmethod
    def from_tfs_file(cls, filename):
        """ Create class instance from TSF file.

        Parameters
        ----------
        filename : str
            Path and name of TSF file.

        Returns
        -------
        Class instance
        """
        with open(filename, 'r') as f:
            tsf_str = f.read()

            # [TSF-file]: single attribute, single occurrence
            expr = r'\[TSF-file\]\s*(?P<tsf_format>\S+)\s*\n'
            expr_count = 0
            for expr_result in re.finditer(expr, tsf_str):
                expr_result_dict = expr_result.groupdict()
                expr_count += 1
            if expr_count == 1:  # OK
                tfs_format = expr_result_dict['tsf_format']
                if tfs_format != 'v01.0':
                    raise RuntimeError(f'Invalid TSF file format: {tfs_format} (valid: v01.0)')
            else:
                raise RuntimeError(f'TSF-file format tag is missing at the begin of the file, or it occurs more than '
                                   f'once!')

            # [TIMEFORMAT]: single attribute, single occurrence
            expr = r'\[TIMEFORMAT\]\s*(?P<timeformat>\S+)\s*\n'
            expr_count = 0
            for expr_result in re.finditer(expr, tsf_str):
                expr_result_dict = expr_result.groupdict()
                expr_count += 1
            if expr_count == 1:  # OK
                timeformat = expr_result_dict['timeformat']
            elif expr_count == 0:
                raise RuntimeError(f'TSF block [TIMEFORMAT] not found in input TSF file!')
            else:
                raise RuntimeError(f'TSF block [TIMEFORMAT] occurs {expr_count} time! Only once allowed.')
            if timeformat not in ('DATETIME', 'DATETIMEFRAC'):
                raise RuntimeError(f'Invalif [TIMEFORMAT] in TSF File: {timeformat}')

            # [UNDETVAL]: single attribute, single occurrence
            expr = r'\[UNDETVAL\]\s*(?P<undetval>\S+)\s*\n'
            expr_count = 0
            for expr_result in re.finditer(expr, tsf_str):
                expr_result_dict = expr_result.groupdict()
                expr_count += 1
            if expr_count == 1:  # OK
                undetval = expr_result_dict['undetval']
            elif expr_count == 0:
                raise RuntimeError(f'TSF block [UNDETVAL] not found in input TSF file!')
            else:
                raise RuntimeError(f'TSF block [UNDETVAL] occurs {expr_count} time! Only once allowed.')
            undetval = float(undetval)

            # [INCREMENT]: single attribute, single occurrence
            expr = r'\[INCREMENT\]\s*(?P<increment>\S+)\s*\n'
            expr_count = 0
            for expr_result in re.finditer(expr, tsf_str):
                expr_result_dict = expr_result.groupdict()
                expr_count += 1
            if expr_count == 1:  # OK
                increment = expr_result_dict['increment']
            elif expr_count == 0:
                raise RuntimeError(f'TSF block [INCREMENT] not found in input TSF file!')
            else:
                raise RuntimeError(f'TSF block [INCREMENT] occurs {expr_count} time! Only once allowed.')
            increment = float(increment)

            # [COUNTINFO]: optional, single attribute, single occurrence
            expr = r'\[COUNTINFO\]\s*(?P<countinfo>\S+)\s*\n'
            expr_count = 0
            for expr_result in re.finditer(expr, tsf_str):
                expr_result_dict = expr_result.groupdict()
                expr_count += 1
            if expr_count == 1:  # OK
                countinfo = int(expr_result_dict['countinfo'])
            else:
                countinfo = -1  # Error code, if the block wasn't found!

            # [COMMENT]: optional, single attribute, single occurrence
            expr = r'\[COMMENT\]\n(?P<comment>[\s\S]*?(?=\n+\[))'
            expr_count = 0
            for expr_result in re.finditer(expr, tsf_str):
                expr_result_dict = expr_result.groupdict()
                expr_count += 1
            if expr_count == 1:  # OK
                comment = expr_result_dict['comment']
                comment = comment.strip()
            else:
                comment = ''

            # [CHANNELS]: multiple lines, single occurrence
            expr = r'\[CHANNELS\]\n(?P<channels>(?: *.+[\n\r]{0,1})+)\s*'
            expr_count = 0
            for expr_result in re.finditer(expr, tsf_str):
                expr_result_dict = expr_result.groupdict()
                expr_count += 1
            if expr_count == 1:  # OK
                channels = expr_result_dict['channels']
            elif expr_count == 0:
                raise RuntimeError(f'TSF block [CHANNELS] not found in input TSF file!')
            else:
                raise RuntimeError(f'TSF block [CHANNELS] occurs {expr_count} time! Only once allowed.')
            channels = channels.splitlines()
            channels = [s.strip() for s in channels]

            # [UNITS]: multiple lines, single occurrence
            expr = r'\[UNITS\]\n(?P<units>(?: *.+[\n\r]{0,1})+)\s*'
            expr_count = 0
            for expr_result in re.finditer(expr, tsf_str):
                expr_result_dict = expr_result.groupdict()
                expr_count += 1
            if expr_count == 1:  # OK
                units = expr_result_dict['units']
            elif expr_count == 0:
                raise RuntimeError(f'TSF block [UNITS] not found in input TSF file!')
            else:
                raise RuntimeError(f'TSF block [UNITS] occurs {expr_count} time! Only once allowed.')
            units = units.splitlines()
            units = [s.strip() for s in units]

            # Create channel metadata:
            channel_metadata = []
            for count, channel in enumerate(channels):
                channel_metadata.append(ChannelMetadata(location=channel.split(':')[0],
                                instrument_name=channel.split(':')[1],
                                data_type=channel.split(':')[2],
                                unit=units[count],
                                channel_number=count+1))

            # [DATA] multiple blocks possible, multiple lines per block
            data_df_list = []
            if timeformat == 'DATETIME':  # Multiple of seconds
                expr = r'\[DATA\]\n(?P<obs_data>(?:[0-9]{4} [0-9]{2} [0-9]{2}  [0-9]{2} [0-9]{2} [0-9]{2}.+[\n\r]{0,1})+)\s*'
                for expr_result in re.finditer(expr, tsf_str):
                    expr_result_dict = expr_result.groupdict()
                    obs_data_str = expr_result_dict['obs_data']
                    data_df_list.append(pd.read_csv(io.StringIO(obs_data_str), delim_whitespace=True, header=None))
                    columns = ['year', 'month', 'day', 'hour', 'minute', 'second'] + [ch.channel_name for ch in
                                                                                channel_metadata]
                    num_channels_in_data_block = data_df_list[0].shape[1] - 6

            elif timeformat == 'DATETIMEFRAC':
                expr = r'\[DATA\]\n(?P<obs_data>(?:[0-9]{4} [0-9]{2} [0-9]{2}  [0-9]{2} [0-9]{2} [0-9]{2} [0-9]{3}.+[\n\r]{0,1})+)\s*'
                for expr_result in re.finditer(expr, tsf_str):
                    expr_result_dict = expr_result.groupdict()
                    obs_data_str = expr_result_dict['obs_data']
                    data_df_list.append(pd.read_csv(io.StringIO(obs_data_str), delim_whitespace=True, header=None))
                    columns = ['year', 'month', 'day', 'hour', 'minute', 'second','ms'] + [ch.channel_name for ch in
                                                                                           channel_metadata]
                    num_channels_in_data_block = data_df_list[0].shape[1] - 7

            # Concat DATA blocks and add column names:
            if len(data_df_list) == 0:
                raise RuntimeError(f'No [DATA] block found in TSF file!')
            elif len(data_df_list) > 0:
                data_df = pd.concat(data_df_list, ignore_index=True)
                data_df.columns = columns

            # Check, if metadata is available for each channel:
            if num_channels_in_data_block != len(channel_metadata):
                raise RuntimeError(f'Number of data channels in [DATA] block ({num_channels_in_data_block}) not equal '
                                   f'to number of channel metadata '
                                   f'sets in the [CHANNELS] and [UNITS] blocks ({len(channel_metadata)})!')

            # Create datetime column:
            if timeformat == 'DATETIME':
                data_df['epoch_dt'] = pd.to_datetime(data_df[['year', 'month', 'day', 'hour', 'minute', 'second']])
            elif timeformat == 'DATETIMEFRAC':
                data_df['epoch_dt'] = pd.to_datetime(data_df[['year', 'month', 'day', 'hour', 'minute', 'second',
                                                              'ms']])

        return cls(filename=filename, tfs_format=tfs_format, timeformat=timeformat, undetval=undetval,
                   increment=increment, countinfo=countinfo, comment=comment, channel_metadata=channel_metadata,
                   data_df=data_df)

        # super().__init__(filename, filetype='tfs')

    @property
    def get_tide_corr_df(self):
        """Returns the complete tidal gravity time-series as pandas DataFrame"""
        return self.data_df


    @property
    def filename(self) -> str:
        """Returns the name of the data file."""
        return self._filename

    @property
    def filetype(self) -> str:
        """Returns the type of the data file."""
        return 'TSoft TSF file'

    @property
    def number_channels(self) -> str:
        """Returns the number of data channels."""
        return len(self.channel_metadata)

    @property
    def number_records(self) -> str:
        """Returns the number of data records."""
        return len(self.data_df)

    @property
    def starttime(self) -> dt.datetime:
        """Returns the start time of the tidal timeseries."""
        return self.data_df['epoch_dt'].min()

    @property
    def endtime(self) -> dt.datetime:
        """Returns the end time of the tidal timeseries."""
        return self.data_df['epoch_dt'].max()

    def get_channel_df(self, channel: [int, str]):
        """Return the specified channel and the time reference (DateTime) as pandas DataFrame.

        Parameters
        ----------
        channel : str or int or 'last'
            If `channel` is as string the argument is interpreted as channel name and the according channel data is
            returned. If `channel` is an integer, it is interpreted as the number of the channel to be returned. `last`
            indicates that the last channel in the TSF file should be returned.

        Returns
        -------
        pandas.DataFrame containing the data records of the selected channel and the time reference as DateTime.
        """
        if isinstance(channel, str):
            if channel == 'last':
                channel = max([num.channel_number for num in self.channel_metadata])
            elif channel in self.data_df.columns and channel in self.channel_names:
                return self.data_df.loc[:, ['epoch_dt', channel]]
        if isinstance(channel, int):
            for ch in self.channel_metadata:
                if ch.channel_number == channel:
                    return self.data_df.loc[:, ['epoch_dt', ch.channel_name]]
        raise RuntimeError(f'"{channel}" is not a valid TSF data record channel name or number!')

    def get_channel_np(self, channel: [int, str]):
        """Return the time reference and the data of the specified channel as numpy array.

        Parameters
        ----------
        channel : str or int or 'last'
            If `channel` is as string the argument is interpreted as channel name and the according channel data is
            returned. If `channel` is an integer, it is interpreted as the number of the channel to be returned. `last`
            indicates that the last channel in the TSF file should be returned.

        Returns
        -------
        numpy.array of reference epochs (DateTime), numpy.array of channel data.
        """
        if isinstance(channel, str):
            if channel == 'last':
                channel = max([num.channel_number for num in self.channel_metadata])
            elif channel in self.data_df.columns and channel in self.channel_names:
                return self.data_df.loc[:, 'epoch_dt'].to_numpy(), self.data_df.loc[:, channel].to_numpy()
        if isinstance(channel, int):
            for ch in self.channel_metadata:
                if ch.channel_number == channel:
                    return self.data_df.loc[:, 'epoch_dt'].to_numpy(), self.data_df.loc[:, ch.channel_name].to_numpy()
        raise RuntimeError(f'"{channel}" is not a valid TSF data record channel name or number!')


    @property
    def channel_names(self):
        """Returns a list of the names of all data channels."""
        return [ch.channel_name for ch in self.channel_metadata]

    @property
    def channel_names(self):
        """Returns a list of the names of all data channels."""
        return [ch.channel_name for ch in self.channel_metadata]

    @property
    def locations(self):
        """Returns a list of the locations of all stations."""
        return [ch.location for ch in self.channel_metadata]

    @property
    def instruments(self):
        """Returns a list of the instruments of all stations."""
        return [ch.instrument_name for ch in self.channel_metadata]

    @property
    def data_types(self):
        """Returns a list of the data types of all stations."""
        return [ch.data_type for ch in self.channel_metadata]

    @property
    def units(self):
        """Returns a list of the units of all stations."""
        return [ch.unit for ch in self.channel_metadata]

    @property
    def filename_without_path(self):
        """Returns the TSF filename without path."""
        return os.path.basename(self.filename)






