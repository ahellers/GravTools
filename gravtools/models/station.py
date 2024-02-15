"""Modelling gravity stations.

Copyright (C) 2021  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

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

import os
import pandas as pd
from gravtools.models.exceptions import FileTypeError
from gravtools.settings import STATION_DATA_SOURCE_TYPES


class Station:
    """Holds information on stations with known attributes.

    Attributes
    ----------
    station_files : dict
        The keys are the filenames (valid filename and path) and the values are the types of the station file
        (see: STATION_DATA_SOURCE_TYPES.keys()).
    stat_df : :py:obj:`pandas.core.frame.DataFrame`
        Pandas Dataframe that holds all station data. The columns are specified in :py:obj:`.Station._STAT_DF_COLUMNS`:

        - station_name : str
            Name of station.
        - long_deg : float
            Geographical longitude of the station [°] as optioned from the data source. Be aware that details on the
            used reference fame are not handled here!
        - lat_deg : float
            Geographical longitude of the station [°] as optioned from the data source. Be aware that details on the
            used reference fame are not handled here!
        - height_deg : float
            Height of the station [m] as optioned from the data source. Be aware that details on the
            used reference fame are not handled here!
        - g_mugal : float
            Gravity at the station [µGal] as optioned from the data source.
        - sd_g_mugal : float
            Standard deviation of the gravity at the station [µGal] as optioned from the data source.
        - vg_mugalm : float, optional (default=np.nan)
            Vertical gradient at the station. If `NaN`, ths value is not available (and will probably be replaced by
            the standard value defined in :py:obj:`gravtools.const.VG_DEFAULT`).
        - is_observed : bool
            Flag, that indicates whether the station was observed at least once in at least one survey in the current
            campaign. If `True`, the station was observed.
        - source_type : str
            The data source type indicates where the station data originates. It hast to be listed in
            :py:obj:`gravtools.settings.STATION_DATA_SOURCE_TYPES`.
        - is_datum : bool
            This flag indicates, whether this station is used as datum station in the analysis. Hence, this is an
            estimation setting. For stations from observation files the default is `False`, for ÖSGN stations the
            default is `True`.
        - in_survey : str
            String that lists the names of the surveys in which the station was observed in (separated by ';')
    """

    _STAT_DF_COLUMNS_DTYPES = {
        'station_name': 'str',  # Station name, str
        'long_deg': 'float',  # longitude [deg], float
        'lat_deg': 'float',  # latitude [deg], float
        'height_m': 'float',  # Height [m]
        'g_mugal': 'float',  # gravity [µGal]
        'sd_g_mugal': 'float',  # standard deviation of the gravity [µGal]
        'vg_mugalm': 'float',  # vertical gradient [µGal/m]
        'is_observed': 'bool',  # flag: True, if station was observed in at least one survey in the current campaign.
        'source_type': 'str',  # Source of the station data, str
        'source_name': 'str',  # Name of the data source
        'is_datum': 'bool',  # flag, that indicates whether this is a datum station, bool
        'in_survey': 'str',  # Names of surveys in which ths station was observed (separator = ;), str
    }
    _STAT_DF_COLUMNS = list(_STAT_DF_COLUMNS_DTYPES.keys())

    # _STAT_DF_COLUMNS = (
    #     'station_name',  # Station name, str
    #     'long_deg',  # longitude [deg], float
    #     'lat_deg',  # latitude [deg], float
    #     'height_m',  # Height [m]
    #     'g_mugal',  # gravity [µGal]
    #     'sd_g_mugal',  # standard deviation of the gravity [µGal]
    #     'vg_mugalm',  # vertical gradient [µGal/m]
    #     'is_observed',  # flag: True, if station was observed in at least one survey in the current campaign.
    #     'source_type',  # Source of the station data, str
    #     'source_name',  # Name of the data source
    #     'is_datum',  # flag, that indicates whether this is a datum station, bool
    #     'in_survey',  # Names of surveys in which ths station was observed (separator = ;), str
    # )

    def __init__(self, station_files=None):
        """
        Parameters
        ----------
        station_files : dict, optional
            The keys are the filenames (valid filename and path) and the values are the types of the station file
            (see: STATION_DATA_SOURCE_TYPES.keys()). Default=None which implies that no stations are added to this
            object.

        Raises
        ------
        FileTypeError
            Station file format not defined in :py:const:`gravtools.settings.STATION_DATA_SOURCE_TYPES`.
        TypeError
            Constructor arguments of wrong type.
        """
        if station_files is None:
            station_files = {}
        if not isinstance(station_files, dict):
            raise TypeError('"station_files" has to be a dict!')

        self.station_files = station_files

        # Init. station dataframe:
        self.stat_df = pd.DataFrame(columns=self._STAT_DF_COLUMNS)
        self.stat_df = self.stat_df.astype(self._STAT_DF_COLUMNS_DTYPES)

        # Loop over station_files and append the content to stat_df:
        for file_name, file_type in station_files.items():
            if file_type not in STATION_DATA_SOURCE_TYPES.keys():
                raise FileTypeError(f'File type of {file_name} not valid.',
                                    valid_file_types=[*STATION_DATA_SOURCE_TYPES])
            # Load file dependent on the file type:
            if file_type == 'oesgn_table':
                self.add_stations_from_oesgn_table(file_name)
            else:
                raise NotImplementedError(f'Loading station files of type "{file_type}" not implemented yet.')

    def _read_oesgn_table(self, filename, is_datum=False):
        """Reads an OESGN table file and returns a dataframe containing this data.

        Notes
        -----
        Currently it is only possible to load station data from the OESGN Table data file!

        Parameters
        ----------
        filename : str
            Name (and path) of the OESGN table file.
        is_datum : bool, optional (default = False)
            `True` indicates that all loaded OESGN stations are initially selected as datum stations (`is_datum`=True)

        Returns
        -------
        :py:obj:`pandas.core.frame.DataFrame`
            Pandas dataframe with all columns specified in :py:attr:`.Station._STAT_DF_COLUMNS` containing the data of
            the OESGN table given by the input argument `filename`.
        """
        widths = (
            10,  # Punktnummer im ÖSGN
            24,  # Anmerkung
            8,  # Geograph. Breite [deg] 34
            8,  # Geograph. Länge [deg] 42
            8,  # Höhe [m] 50
            7,  # g - 9.8e8 [µGal] 58
            3,  # mittlerer Fehler von g [µGal] 65
            4,  # Vertikalgradient der Schwere [µGal/m] 68
            6,  # Datum 73
            12,  # Punktidentität 79
        )
        column_names = (
            'station_name',
            'notes',
            'lat_deg',
            'long_deg',
            'height_m',
            'g_mugal',
            'sd_g_mugal',
            'vg_mugalm',
            'date',
            'identity'
        )

        columns_to_be_dropped = []
        for col_name in column_names:
            if col_name not in self._STAT_DF_COLUMNS:
                columns_to_be_dropped.append(col_name)

        try:
            stat_df_oesgn = pd.read_fwf(filename, widths=widths, header=None)
        except UnicodeDecodeError:
            # To catch encoding errors with files created/converted on Windows.
            # See: https://docs.python.org/3/library/codecs.html#standard-encodings
            # See: https://stackoverflow.com/questions/45529507/unicodedecodeerror-utf-8-codec-cant-decode-byte-0x96-in-position-35-invalid
            stat_df_oesgn = pd.read_fwf(filename, widths=widths, header=None, encoding='ISO-8859-1')

        stat_df_oesgn.columns = column_names
        stat_df_oesgn['is_observed'] = False  # Default = False
        if is_datum:
            stat_df_oesgn['is_datum'] = True
        else:
            stat_df_oesgn['is_datum'] = False
        stat_df_oesgn['source_type'] = 'oesgn_table'
        stat_df_oesgn['source_name'] = os.path.basename(filename)  # filename

        stat_df_oesgn.drop(columns=columns_to_be_dropped, inplace=True)

        stat_df_oesgn['height_m'] = stat_df_oesgn['height_m'] * 1e-3  # Conversion: [mm] => [m]

        return stat_df_oesgn

    def add_stations_from_oesgn_table(self, filename, is_datum=False, verbose=False):
        """Adds (appends) all stations of the input OESGN table to the station dataframe.

        Parameters
        ----------
        filename : str
            Name (and path) of the OESGN table file.
        is_datum : bool, optional (default = False)
            `True` indicates that all loaded OESGN stations are initially selected as datum stations (is_datum=`True`)
        verbose : bool, optional
            Print notifications, if `True` (default=`False`)
        """
        stat_df_new = self._read_oesgn_table(filename, is_datum=is_datum)
        stat_df_new = self._stat_df_add_columns(stat_df_new)
        stat_df_new = self._stat_df_reorder_columns(stat_df_new)
        self.add_stations(stat_df_new, data_source_type='oesgn_table', verbose=verbose)

    def delete_station(self, station_names: list, verbose=False):
        """
        Deletes all rows where the station name occurs in the `station_names` list.

        Parameters
        ----------
        station_names : list of str
            List of station names. The according station records will be deleted from `self.stat_df`.
        verbose : bool
            If `True`, print notification on deleted items.
        """
        idx_series = self.stat_df.station_name.isin(station_names)
        if verbose:
            deleted_df = self.stat_df[idx_series]
            print(f'Deleted {len(deleted_df)} rows in station dataframe:')
            for index, row in deleted_df.iterrows():
                print(f' - Station {row["station_name"]:10s} (row index {index:d})')
        self.stat_df = self.stat_df[~idx_series]

    @property
    def get_number_of_stations(self) -> int:
        """int : Returns the number of stations."""
        return len(self.stat_df)

    def __str__(self):
        return f"Station DataFrame containing {self.get_number_of_stations} stations."

    @property
    def get_all_stations(self):
        """:py:obj:`pandas.core.frame.DataFrame` : Returns all available stations."""
        return self.stat_df

    def set_observed_info_from_survey(self, survey):
        """Set the `is_observed` flags and the `in_survey` strings in the py:obj:`Station.stat_df` DataFrame.

        Parameters
        ----------
        survey : py:obj:`Survey`
            Survey object which is used as reference to set the `is_observed` flags in the stations Dataframe.
        """
        self.stat_df.loc[self.stat_df['station_name'].isin(survey.obs_df['station_name']), 'is_observed'] = True

        # Delete None entries:
        self.stat_df.loc[self.stat_df['station_name'].isin(survey.obs_df['station_name']) & self.stat_df[
            'in_survey'].isnull(), 'in_survey'] = ''
        # self.stat_df.loc[self.stat_df['station_name'].isin(survey.obs_df['station_name']), 'in_survey'] = ''

        # Add survey names:
        num_entries = len(
            self.stat_df.loc[self.stat_df['station_name'].isin(survey.obs_df['station_name']), 'in_survey'])
        self.stat_df.loc[self.stat_df['station_name'].isin(survey.obs_df['station_name']), 'in_survey'] = \
            self.stat_df.loc[self.stat_df['station_name'].isin(survey.obs_df['station_name']),
                             'in_survey'].str.cat([survey.name + '; '] * num_entries)

    @classmethod
    def _stat_df_add_columns(cls, stat_df):
        """Add and initialize (as None) all columns that are listed in cls._STAT_DF_COLUMNS and not present in input
        stat_df."""
        columns_to_be_initialized_as_none = list(set(cls._STAT_DF_COLUMNS) - set(stat_df.columns))
        stat_df[columns_to_be_initialized_as_none] = None
        return stat_df

    @classmethod
    def _stat_df_reorder_columns(cls, stat_df):
        """Change order of columns of stat_df to the order specified in cls._STAT_DF_COLUMNS.

        See: https://erikrood.com/Python_References/change_order_dataframe_columns_final.html
        """
        stat_df = stat_df[list(cls._STAT_DF_COLUMNS)]
        return stat_df

    @classmethod
    def _stat_df_check_columns(cls, stat_df, verbose=False) -> bool:
        """Check if all columns specified in cls._STAT_DF_COLUMNS are present."""
        is_valid = True

        invalid_cols = list(set(stat_df) - set(cls._STAT_DF_COLUMNS))
        if len(invalid_cols) > 0:
            is_valid = False
            if verbose:
                print(f'The following columns are not valid: {", ".join(invalid_cols)}')

        invalid_cols = list(set(cls._STAT_DF_COLUMNS) - set(stat_df))
        if len(invalid_cols) > 0:
            is_valid = False
            if verbose:
                print(f'The following columns are missing: {", ".join(invalid_cols)}')

        return is_valid

    def add_stations_from_survey(self, survey, verbose=False):
        """Add stations from a survey dataset, if they are not included in the station data yet.

        Notes
        -----
        The location information (longitude, latitude and height) are taken from the first appearance of a station in
        the survey dataframe, although this information may differ between setups if the location was determined (e.g.
        by GPS) for each individual measurement setup.

        Parameters
        ----------
        survey : py:obj:`Survey`
            Survey object which is search for stations that are not included in the station dataframe yet.
        verbose : bool, defualt=False
            If `True`, print notifications to command line.
        """
        # Get the essential information from the input survey:
        tmp_df = survey.obs_df[['station_name', 'lon_deg', 'lat_deg', 'alt_m', 'obs_epoch']].copy(deep=True)
        tmp_df = tmp_df.sort_values(by=['obs_epoch']).drop_duplicates(subset=['station_name']).drop(
            columns=['obs_epoch'])

        # Add missing columns (defined in self._STAT_DF_COLUMNS) to temp_df and initialize data:
        tmp_df.rename(columns={'lon_deg': 'long_deg', 'alt_m': 'height_m'}, inplace=True)
        tmp_df = self._stat_df_add_columns(tmp_df)
        tmp_df = self._stat_df_reorder_columns(tmp_df)
        tmp_df['is_datum'] = False
        tmp_df['is_observed'] = True
        tmp_df['source_type'] = 'obs_file'
        tmp_df['source_name'] = survey.name  # survey name
        tmp_df['in_survey'] = ''
        tmp_df.reset_index(drop=True, inplace=True)
        tmp_df = tmp_df.astype(self._STAT_DF_COLUMNS_DTYPES)

        self.add_stations(tmp_df, data_source_type='obs_file', verbose=verbose)

        # # Add missing stations to stat_df by matching the station names only (location parameters of observed stations
        # # may differ!):
        # self.stat_df = pd.concat([self.stat_df, tmp_df]).drop_duplicates(
        #     subset=['station_name'],
        #     keep='first'
        # ).reset_index(drop=True)

    def add_stations(self, stat_df_add, data_source_type, verbose=False):
        """Adds stations from various sources (e.g. station and observation files) to the station dataframe.

        Notes
        -----

        - Station names (column 'station_name') have to be unique in the resulting station dataframe
        - Keep the added station and drop existing entries, if the station names match.
        - If a station (with a unique station_name) is already available in the station dataframe and originates from a
            station file (e.g. ÖSGN file), the new station (from any source) is not added.
        - If a station originating from an observation file is already available in the station dataframe an the same
            station (same name) should be added from a station file (e.g. ÖSGN), the old version is overwritten be the
            entry from the station file.
        - If a station originating from an observation file is already available in the station dataframe and the same
            station (same name) should be added from an observation file (e.g. CG5 observation file), the old version
            is kept and the new version is discarded.

        Parameters
        ----------
        stat_df_add : py:oby:`.Station.stat_df`
            Station dataframe, that should be added so `self.stat_df`. The new station data has to come from one source
            (specified in `data_source_type`).
        data_source_type : str
            Specifies the type of the data source. Valid data source types are defined in
            :py:obj:`gravtools.settings.STATION_DATA_SOURCE_TYPES`.
        verbose : bool, default = `False`
            If `True`, the status messages are printed to the command line.

        """
        number_of_existing_stations = len(self.stat_df)

        # check, if the input station dataframe has valid columns:
        if not self._stat_df_check_columns(stat_df_add, verbose=verbose):
            raise AssertionError('The columns of the input dataframe are not valid!')

        # Check if data source type is valid:
        if data_source_type not in STATION_DATA_SOURCE_TYPES:
            raise AssertionError(f'The station data source type "{data_source_type}" is not valid!')

        # Check if the data in the input dataframe comes from the (valid) source given as input parameter:
        if not all(stat_df_add['source_type'] == data_source_type):
            raise AssertionError(f'All entries in the input station dataframe need to have the same data source '
                                 f'type ("{data_source_type}")!')

        # Add stations to a temp. copy of self.stat_df:
        stat_df_old = self.stat_df.copy(deep=True)

        # Add missing stations to stat_df by matching the station names only (location parameters of observed stations
        # may differ!):
        # - From station files => Drop existing entries with the same name in any case!
        if data_source_type == 'oesgn_table' or data_source_type == 'csv_stat_file':
            stat_df_new = pd.concat([stat_df_old, stat_df_add]).drop_duplicates(
                subset=['station_name'],
                keep='last'
            ).reset_index(drop=True)
        # From observation file  => Drop existing entries with the same name!
        elif data_source_type == 'obs_file':
            # => Replace entries with data_source = obs_file:
            stat_df_new = pd.concat([stat_df_old, stat_df_add]).drop_duplicates(
                subset=['station_name', 'source_type'],
                keep='last'
            ).reset_index(drop=True)
            # => Drop new entries from obs files, if entries from another data source type (station file) with the same
            # name already exist:
            stat_df_new = pd.concat([stat_df_new, stat_df_add]).drop_duplicates(
                subset=['station_name'],
                keep='first'
            ).reset_index(drop=True)

        self.stat_df = stat_df_new

        if verbose:
            number_of_new_stations = len(stat_df_add)
            number_of_stations = len(self.stat_df)
            stations_added = number_of_stations - number_of_existing_stations
            print(f"{number_of_new_stations} stations loaded. "
                  f"{stations_added} stations added ({number_of_new_stations - stations_added} already listed).")

    def add_stations_from_csv_file(self, filename, verbose=False):
        """Adds (appends) all stations of the input OESGN table to the station dataframe.

        Parameters
        ----------
        filename : str
            Name (and path) of the OESGN table file.
        verbose : bool, optional
            Print notifications, if `True` (default=`False`)
        """
        if verbose:
            print(f'Read station csv file: {filename}')
        stat_df_new = pd.read_csv(filename, comment='#')
        stat_df_new = self._stat_df_add_columns(stat_df_new)
        stat_df_new = self._stat_df_reorder_columns(stat_df_new)
        stat_df_new['is_observed'] = False  # Default = False
        stat_df_new['is_datum'] = False
        stat_df_new['source_type'] = 'csv_stat_file'
        stat_df_new['source_name'] = os.path.basename(filename)
        self.add_stations(stat_df_new, data_source_type='csv_stat_file', verbose=verbose)

    def set_datum_stations(self, station_names: list, is_datum: bool, verbose=False):
        """
        Dependent on the `is_datum` flag set the datum status of all stations in the list

        Parameters
        ----------
        station_names : list of str
            List of station names. The according station records will be deleted from `self.stat_df`.
        is_datum : bool
            `True` implies that all stations in the `station_names` list are set to datum stations. `False` implies
            that all stations in the list are no datum stations.
        verbose : bool
            If `True`, print notification to command line.
        """
        idx_series = self.stat_df.station_name.isin(station_names)
        if verbose:
            station_list = self.stat_df.loc[idx_series, 'station_name'].to_list()
            if is_datum:
                print(f'The following stations are defined as datum stations:')
            else:
                print(f'The following stations are NO datum stations:')
            for station_name in station_list:
                print(f' - {station_name}')
        self.stat_df.loc[idx_series, 'is_datum'] = is_datum
