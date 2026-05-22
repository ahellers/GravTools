"""Abstract base class for tide data.

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
from abc import ABC, abstractmethod
import datetime as dt

import pandas as pd


class AbstractTideData(ABC):
    """Abstract base class for tidal data."""

    def get_tidal_gravity_correction(self, epoch_dt):
        pass

    @property
    @abstractmethod
    def get_tide_corr_df(self):
        """Returns the complete tidal gravity time-series as pandas DataFrame"""
        pass

    @property
    @abstractmethod
    def filename(self) -> str:
        """Returns the name of the data file."""
        pass

    @property
    @abstractmethod
    def filetype(self) -> str:
        """Returns the type of the data file."""
        pass

    @property
    @abstractmethod
    def number_channels(self) -> str:
        """Returns the number of data channels."""
        pass

    @property
    @abstractmethod
    def number_records(self) -> str:
        """Returns the number of data records."""
        pass

    @property
    @abstractmethod
    def starttime(self) -> dt.datetime:
        """Returns the start time of the tidal timeseries."""
        pass

    @property
    def starttime_str(self) -> str:
        """Returns the start time of the tidal timeseries as string."""
        return self.starttime.strftime('%Y-%d-%d %H:%M:%S')

    @property
    @abstractmethod
    def endtime(self) -> dt.datetime:
        """Returns the end time of the tidal timeseries."""
        pass

    @abstractmethod
    def get_channel_df(self, channel: [int, str]) -> pd.DataFrame:
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
        pass

    @property
    def endtime_str(self) -> str:
        """Returns the end time of the tidal timeseries as string."""
        return self.endtime.strftime('%Y-%d-%d %H:%M:%S')

    @property
    def filename_without_path(self) -> str:
        """Returns the filename without path."""
        pass

    def __str__(self):
        return f'Tidal data loaded from {self.filename} ({self.filetype}) with {self.number_channels} channels ' \
               f'and {self.number_records} datasets ({self.starttime_str} to {self.endtime_str})'




