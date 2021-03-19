"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
Contains classes for modeling gravity campaigns.
"""

import pandas as pd
import numpy as np
import datetime as dt
import os

from gravtools.settings import SURVEY_DATA_SOURCE_TYPES, STATION_DATA_SOURCE_TYPES, GRAVIMETER_ID_BEV, \
    TIDE_CORRECTION_TYPES, DEFAULT_GRAVIMETER_ID_CG5_SURVEY, REFERENCE_HEIGHT_TYPE, NAME_OBS_FILE_BEV, \
    PATH_OBS_FILE_BEV, BEV_GRAVIMETER_TIDE_CORR_LOOKUP, GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m
from gravtools.const import VG_DEFAULT
from gravtools.models.exceptions import FileTypeError
from gravtools.CG5_utils.cg5_survey import CG5Survey
from gravtools.models.lsm import LSM, LSMDiff


class Station:
    """
    Holds information on stations with known attributes.

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

    _STAT_DF_COLUMNS = (
        'station_name',  # Station name, str
        'long_deg',  # longitude [deg], float
        'lat_deg',  # latitude [deg], float
        'height_m',  # Height [m]
        'g_mugal',  # gravity [µGal]
        'sd_g_mugal',  # standard deviation of the gravity [µGal]
        'vg_mugalm',  # vertical gradient [µGal/m]
        'is_observed',  # flag: True, if station was observed in at least one survey in the current campaign.
        'source_type',  # Source of the station data, str
        'source_name',  # Name of the data source
        'is_datum',  # flag, that indicates whether this is a datum station, bool
        'in_survey',  # Names of surveys in which ths station was observed (separator = ;), str
    )

    # _STAT_DF_COLUMNS_SAME_STATION_INDICATORS = (
    #     'station_name',  # Station name, str
    #     'long_deg',  # longitude [deg], float
    #     'lat_deg',  # latitude [deg], float
    #     'height_m',  # Height [m]
    #     'source_type',
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

        # Loop over station_files and append the content to stat_df:
        for file_name, file_type in station_files.items():
            if file_type not in STATION_DATA_SOURCE_TYPES.keys():
                raise FileTypeError(f'File type of {file_name} not valid.',
                                    valid_file_types=[*STATION_DATA_SOURCE_TYPES])
            # Load file dependent on the file type:
            if file_type == 'oesgn_table':
                self.add_stations_from_oesgn_table(file_name)
            else:
                raise NotImplementedError(f'Loading station files of type "{file_type}" not implemneted yet.')

    def _read_oesgn_table(self, filename):
        """Reads an OESGN table file and returns a dataframe containing this data.

        Notes
        -----
        Currently it is only possible to load station data from the OESGN Table data file!

        Parameters
        ----------
        filename : str
            Name (and path) of the OESGN table file.

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
            7,  # g [µGal] 58
            3,  # mittlerer Fehler von g [?] 65
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

        stat_df_oesgn = pd.read_fwf(filename, widths=widths, header=None)
        stat_df_oesgn.columns = column_names
        stat_df_oesgn['is_observed'] = False  # Default = False
        stat_df_oesgn['is_datum'] = True
        stat_df_oesgn['source_type'] = 'oesgn_table'
        stat_df_oesgn['source_name'] = os.path.basename(filename)  # filename

        stat_df_oesgn.drop(columns=columns_to_be_dropped, inplace=True)

        stat_df_oesgn['height_m'] = stat_df_oesgn['height_m'] * 1e-3  # Conversion: [mm] => [m]

        return stat_df_oesgn

    def add_stations_from_oesgn_table(self, filename, verbose=False):
        """Adds (appends) all stations of the input OESGN table to the station dataframe.

        Parameters
        ----------
        filename : str
            Name (and path) of the OESGN table file.
        verbose : bool, optional
            Print notifications, if `True` (default=`False`)
        """
        stat_df_new = self._read_oesgn_table(filename)
        stat_df_new = self._stat_df_add_columns(stat_df_new)
        stat_df_new = self._stat_df_reorder_columns(stat_df_new)
        self.add_stations(stat_df_new, data_source_type='oesgn_table', verbose=verbose)

        # number_of_existing_stations = len(self.stat_df)
        #
        # # Concatenate multiple df without duplicate rows and a unique indices:
        # # https://stackoverflow.com/questions/21317384/pandas-python-how-to-concatenate-two-dataframes-without-duplicates
        # # https://learndataanalysis.org/concatenate-pandas-dataframes-without-duplicates/
        # # - pd.concat([self.stat_df, stat_df_new]) => simply concat both df, while it may result in duplicate indices
        # #   and rows.
        # # - drop_duplicates() => Drop duplicate rows
        # # - reset_index(drop=True) => Reset index creates a new indec col with unique indices. With drop=True the old
        # #   index column is dropped.
        #
        # # Concatenate existing with new stations and drop duplicates, based on the columns in specified in
        # # _STAT_DF_COLUMNS_SAME_STATION_INDICATORS:
        # self.stat_df = pd.concat([self.stat_df, stat_df_new]).drop_duplicates(
        #     subset=self._STAT_DF_COLUMNS_SAME_STATION_INDICATORS,
        #     keep='first'
        # ).reset_index(drop=True)
        #
        # # Check, if station names still have multiple occurrences in the station DataFrame:
        # number_of_occurrences_per_station_name_series = self.stat_df.value_counts(subset=['station_name'])
        # if sum(number_of_occurrences_per_station_name_series > 1):
        #     duplicate_stations = number_of_occurrences_per_station_name_series[
        #         number_of_occurrences_per_station_name_series > 1
        #     ].index.values  # The index contains the names of stations with multiple occurrences
        #     duplicate_stations_list = []
        #     for item in duplicate_stations:  # Convert list of tuples of strings to list of strings
        #         duplicate_stations_list.append(item[0])
        #     duplicate_stations_str = ', '.join(duplicate_stations_list)
        #     raise AssertionError(f'When adding news stations the following station already exist with differing '
        #                          f'location parameters and/or a differing data source:' + duplicate_stations_str + f'.'
        #                          f'\nPlease load station files before adding observations to the campaign!')
        #
        # if verbose:
        #     number_of_new_stations = len(stat_df_new)
        #     number_of_stations = len(self.stat_df)
        #     stations_added = number_of_stations - number_of_existing_stations
        #     print(f"{number_of_new_stations} stations loaded from {filename}. "
        #           f"{stations_added} stations added ({number_of_new_stations - stations_added} already existed).")

    def delete_station(self, station_names, verbose=False):
        """
        Deletes all rows with the observed station's name occurs in the `station_names` list.

        Parameters
        ----------
        station_names : list of str
            List of station names. The related station records will be deleted from `self.stat_df`.
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

    def set_observed_info_from_survey(self, survey, verbose=False):
        """Set the `is_observed` flags and the `in_survey` strings in the py:obj:`Station.stat_df` DataFrame.

        Parameters
        ----------
        survey : py:obj:`Survey`
            Survey object which is used as reference to set the `is_observed` flags in the stations Dataframe.
        verbose : bool, optional
            Print notifications, if `True` (default=`False`)
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
        tmp_df.reset_index(drop=True, inplace=True)

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
        - Keep the added station and drop existing entries, if the sation names match.
        - If a station (with a unique station_name) is already available in the station dataframe and originates from a
            station file (e.g. ÖSGN file), the new station (from any source) is not added.
        - If a station originating from an observation file is already available in the station dataframe an the same
            station (same name) should be added from a station file (e.g. ÖSGN), the old version is overwritten be the
            entry from the staton file.
        - If a station originating from an observation file is already available in the station dataframe an the same
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
        verbose : bool, delauft = `False`
            If `True`, the status mesages are printed to the command line.

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
        if data_source_type == 'oesgn_table':
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


class Survey:
    """Gravity survey (instrument-independent).

    A gravity survey object contains all data that belongs to a single field observation job, carried out under the
    same circumstances (same instrument, measurement area, drift- and datum-point planing, same observer, etc.). A
    survey is usually observed on the same day, under similar conditions by the same operator.

    All analytical observation-level reductions and corrections are applied here. The reduced g values are stored in the
    columns `g_red_mugal` and the according standard deviation in `sd_g_red_mugal` of `obs_df`. The following
    corrections/reductions are supported:

    - Tidal corrections of observations

      - Bulit-in corrections of the CG-5 (Longman, 1959)

    - Reduction of the observed gravity to various height levels

      - Control point level
      - Ground level
      - Level of the gravimeter top
      - Sensor level

    Notes
    -----
    Basically it is possible to initialize an empty survey, by just defining the survey's name on the instantiation
    step. In this case all other (class and instance) attributes are initialized with the default values listed below
    in the Attributes section.

    Attributes
    ----------
    name : str
        Name of the survey. This parameter is mandatory!
    date :
        Date of te survey.
    operator : str, optional (default='')
        Name of the responsible operator that carried out the observations of this survey.
    gravimeter_id : str, optional (default='')
        ID of the used gravimeter. Valid gravimeter IDs have to be listed in
        :py:obj:`gravtools.settings.GRAVIMETER_ID_BEV`.
    data_file_name : str, optional (default='')
        Name of the data source file (observation file). Without file path!
    data_file_type : str, optional (default='')
        Type of the data source file. Since gravtools allows to load data from different sources, it is important to
        track the data source type.
    obs_tide_correction_type : str, (default='')
        Type of the tidal corrections applied on the observations (gravimeter readings as obtained from the data
        source; column `g_obs_mugal` in `obs_df`). Valid entries have to be listed in
        :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`.
    obs_reference_height_type : str, (default='')
        Reference level type of the observations (gravimeter readings as obtained from the data
        source; column `g_obs_mugal` in `obs_df`). Valid entries have to be listed in
        :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`.
    red_tide_correction_type : str, optional (default='')
        Type of the tidal corrections applied on the reduced observations (column `g_red_mugal` in `obs_df`). Valid
        entries have to be listed in :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`. `Empty string`, if corrected
        observations are not available.
    red_reference_height_type : str, optional (default='')
        Reference level type of the reduced observations (column `g_red_mugal` in `obs_df`). Valid entries have to be
        listed in :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`. `Empty string`, if corrected observations are not
        available.
    keep_survey : bool (default=True)
        Flag that indicates whether this survey will be considered in the data analysis (adjustment).`True` is the
        default and implies that this survey is considered.
    obs_df : :py:obj:`pandas.core.frame.DataFrame`, optional (default=None)
        Contains all observation data that belongs to this survey. One observation per line. If `None`, no observations
        have been assigned to the survey. All Columns are listed in :py:obj:`.Survey._OBS_DF_COLUMNS`:

        - station_name : str
            Name of observed station
        - setup_id : int (UNIX timestamp of the first observation reference epoch of this setup)
            Each setup (one or more observations at a station without moving the instrument) gets a unique ID in order
            to distinguish between independent groups of observations at a station.
        - loop_id: int, optional (default=None)
            Unique ID of a line. A survey can be split up into multiple lines. A line needs to have at least one station
            that was observed at least twice for drift control (preferably at the begin and end of the line). If `None`,
            loops were not defined. The purpose of splitting surveys into loops is to carry out drift correction for
            individual loops (shorter time period) rather than fo the complete survey.
        -  lon_deg : float, optional (default=None)
            Geographical longitude of the station [°]. When loading observation data from the CG-5 observation files,
            geographical coordinates are obtained from measurements of the built-in GPS device of teh instrument.
        -  lat_deg : float, optional (default=None)
            Geographical latitude of the station [°]. When loading observation data from the CG-5 observation files,
            geographical coordinates are obtained from measurements of the built-in GPS device of teh instrument.
        -  alt_m : float, optional (default=None)
            Altitude of the station [m]. When loading observation data from the CG-5 observation files,
            altitudes are obtained from measurements of the built-in GPS device of teh instrument.
        - obs_epoch : :py:obj:`datetime.datetime`; timezone aware, if possible
            Reference epoch of the observation. Per default the start epoch of an observation. Be aware that for the
            determination of tidal corrections by the CG-5 built-in model (Longman, 1959) the middle of the observation
            with the duration dur_sec [sec] is used (obs_epoch + dur_sec/2)!
        - g_obs_mugal : float
            Observed gravity value (instrument reading) [µGal], as obtained from the data source (observation files).
            Be aware that, depending on the instrument settings and the observation data source, different corrections
            and/or reductions may have been applied already! The reference height type and the applied tidal correction
            have to be concise with the statements in :py:obj:`.Survey.obs_reference_height_type` and
            :py:obj:`.Survey.obs_tide_correction_type`, respectively.
        - sd_g_obs_mugal : float, optional (default=None)
            Standard deviation of `g_obs_mugal`. If `None`, the standard deviation is not available.
        - g_red_mugal : float, optional (default=None)
            Reduced gravity observation. This value is derived from the observed gravimeter reading (`g_obs_mugal`) by
            applying reductions and corrections, e.g. from analytical models for tidal effects, or by a reduction to a
            different height level (using the vertical gravity gradient). Take care, that the same  reductions and
            corrections are not applied more than once! The reference height type and the applied tidal correction
            have to be concise with the statements in :py:obj:`.Survey.red_reference_height_type` and
            :py:obj:`.Survey.red_tide_correction_type`, respectively. If `None`, no reductions and/or corrections have
            been applied so far.
        - sd_g_red_mugal : float, optional (default=None)
            Standard deviation of the reduced gravity reading (`g_red_mugal`). This value is derived from the standard
            deviation of the gravity reading (`sd_g_obs_mugal`) by applying proper error propagation. If `None`, the SD
            of the gravity reading is not available, or no reductions/corrections have been applied so far, or an error
            propagation model is still missing. In general, this value should not be `None`, if `g_red_mugal` is not
            `None`.
        - corr_terrain : float, optional (default=None)
            Terrain correction [??] as determined by the built-in model of the Scintrex CG-5. If `None`, this
            correction is not available in the observation data.
        - corr_tide_mugal : float, optional (default=None)
            Tidal correction [µGal] as determined by the built-in model of the Scintrex CG-5 (Longman, 1959). Be aware
            that the tidal corrections by the CG-5 model is determined fot the middle of the observation (also see
            `obs_epoch`). If `None`, this correction is not available in the observation data.
        - temp : float, optional (default=None)
            Temperature [mK] as determined by the Scintrex CG-5. Be aware that this is not the ambient temperature!
            This temperature is obtained from the GC-5 observation file and indicates internal temperature variations
            that are measured and used for the instrumental temperature compensation. If `None`, this information is
            not available in the observation data.
        - tiltx : float, optional (default=None)
            Tilt in x-direction [arcsec] of the gravimeter logged during the observation. If `None`, this
            information is not available in the observation data.
        - tilty : float, optional (default=None)
            Tilt in y-direction [arcsec] of the gravimeter logged during the observation. If `None`, this
            information is not available in the observation data.
        - dhf_m : float
            Vertical distance between instrument top and physical reference point [m]. This information is required
            for reducing the observed gravity to the reference point level (vertical gravity gradient also required).
        - dhb_m : float
            Vertical distance between instrument top and the ground [m]. This information is required
            for reducing the observed gravity to the ground level (vertical gravity gradient also required).
        - keep_obs : bool (default=True)
            Flag, that indicates whether this observation should be considered in the data analysis (adjustment).
            If `True`, the observations takes part in the analysis.
        - vg_mugalm : float, optional (default=None)
            Vertical gravity gradient at the station [µGal/m], obtained from an external source (e.g. station info
            file). The vertical gradient is required for reducing the observed gravity to different reference heights.
        - corr_tide_red_mugal : float, optional (default=None)
            Tidal correction [µGal] that is applied to `g_red_mugal`.

    setup_df : :py:obj:`pandas.core.frame.DataFrame`, optional (default=None)
        Contains a single pseudo observation per setup calculated as variance weighted mean of all active
        (flag `keep_obs = True`) reduced observations of a setup and other information that is required for the
        subsequent parameter adjustment. If `None`, setup data ha not been determined yet. All Columns are listed in
        :py:obj:`.Survey._SETUP_DF_COLUMNS`:

        - station_name : str
            Name of the station that is observed in the setup.
        - setup_id : int
            Same as in `obs_df`.
        - g_mugal : float
            Variance weighted mean of all active reduced observations the setup [µGal].
        - sd_g_red_mugal : float
            Standard deviation of `g_mugal` [µGal]
        - obs_epoch : :py:obj:`datetime.datetime`; timezone aware, if possible
            Reference epoch of `g_mugal`.
    """

    _OBS_DF_COLUMNS = (
        'station_name',  # Name of station (str)
        'setup_id',  # unique ID of setup (int)
        'loop_id',  # Line ID, optional (default=None) (int)
        'lon_deg',  # Longitude [deg], optional (float)
        'lat_deg',  # Latitude [deg], optional (float)
        'alt_m',  # Altitude [m], optional (float)
        'obs_epoch',  # Observation epoch (datetime object, TZ=<UTC>)
        'g_obs_mugal',  # observed g from obs file [µGal] (float)
        'sd_g_obs_mugal',  # Standard deviation of g observed from obs file (float) [µGal]
        'g_red_mugal',  # Reduced gravity observation at station (float) [µGal]
        'sd_g_red_mugal',  # Standard deviation of the reduced gravity (float) [µGal]
        'corr_terrain',  # Terrain correction [??]
        'corr_tide_mugal',  # Tidal correction loaded from input file [mGal], optional (e.g. from CG5 built-in model)
        'temp',  # Temperature [mK], optional
        'tiltx',  # [arcsec], optional
        'tilty',  # [arcsec], optional
        'dhf_m',  # Distance between instrument top and physical reference point (float) [m]
        'dhb_m',  # Distance between instrument top and ground (float) [m]
        'keep_obs',  # Remove observation, if false (bool)
        'vg_mugalm',  # vertical gradient [µGal/m]
        'corr_tide_red_mugal',  # Alternative tidal correction [mGal], optional
    )

    _SETUP_DF_COLUMNS = (
        'station_name',  # Name of station (str)
        'setup_id',  # Unique ID of setup (int)
        'g_mugal',  # Variance weighted mean of all active observations in setup [µGal] (float)
        'sd_g_mugal',  # Standard deviation of `g_mugal` (float) [µGal]
        'epoch_unix',  # Reference epoch of `g_mugal` (unix time [sec])
        'epoch_dt',  # Reference epoch of `g_mugal` (datetime obj)
    )

    # TODO: Get missing infos on columns in CG5 obs file!
    # How is the reference epoch of the observations defined?
    # - begin, middle or end of the observation?
    # - Does this depend on the obs-file type?
    # => Enter the missing information in the docstring above!

    def __init__(self,
                 name,
                 date=None,
                 operator='',
                 institution='',
                 gravimeter_id='',
                 data_file_name='',
                 data_file_type='',
                 obs_df=None,
                 obs_tide_correction_type='',  # of "g_obs_mugal"
                 obs_reference_height_type='',  # of "g_obs_mugal"
                 red_tide_correction_type='',  # of "g_red_mugal"
                 red_reference_height_type='',  # of "g_red_mugal"
                 keep_survey=True,  # Flag
                 setup_df=None
                 ):
        """Default constructor of class Survey."""

        # Check input arguments:
        # name:
        if name is not None:
            if isinstance(name, str):
                if len(name) > 0:
                    self.name = name
                else:
                    raise ValueError('"name" needs to be a non-empty string')
            else:
                raise TypeError('"name" needs to be a non-empty string')

        # date:
        if date is not None:
            if isinstance(date, dt.date):
                self.date = date
            else:
                raise TypeError('"date" needs to be a datetime object')
        else:
            self.date = date  # None

        # operator:
        if isinstance(operator, str):
            self.operator = operator
        else:
            raise TypeError('"operator" needs to be a string')

        # institution:
        if isinstance(institution, str):
            self.institution = institution
        else:
            raise TypeError('"institution" needs to be a string')

        # gravimeter_id:
        if isinstance(gravimeter_id, str):
            if gravimeter_id:
                if gravimeter_id in GRAVIMETER_ID_BEV.keys():
                    self.gravimeter_id = gravimeter_id
                else:
                    raise ValueError('"gravimeter_id" needs to be a key in GRAVIMETER_ID_BEV')
            else:
                self.gravimeter_id = gravimeter_id  # None
        else:
            raise TypeError('"gravimeter_id" needs to be a string')

        # data_file_name:
        if isinstance(data_file_name, str):
            self.data_file_name = data_file_name
        else:
            raise TypeError('"data_file_name" needs to be a string')

        # data_file_type:
        if isinstance(data_file_type, str):
            if data_file_type:
                if data_file_type in SURVEY_DATA_SOURCE_TYPES.keys():
                    self.data_file_type = data_file_type
                else:
                    raise ValueError('"gravimeter_id" needs to be a key in GRAVIMETER_ID_BEV '
                                     '({})'.format(', '.join(GRAVIMETER_ID_BEV.keys())))
            else:
                if self.data_file_name:  # Not empty
                    raise ValueError('If a data file is specified ("data_file_name" not empty, "data_file_type" has '
                                     'to be specified also.')
                self.data_file_type = data_file_type
        else:
            raise TypeError('"data_file_type" needs to be a string')

        # obs_df:
        if obs_df is not None:
            if isinstance(obs_df, pd.DataFrame):
                # Check if obs_df contains exactly all columns defined by self._OBS_DF_COLUMNS:
                if all([item for item in obs_df.columns.isin(self._OBS_DF_COLUMNS)]) and \
                        obs_df.shape[1] == len(self._OBS_DF_COLUMNS):
                    self.obs_df = obs_df
                else:
                    raise ValueError('"obs_df" needs the following columns:{}'.format(', '.join(self._OBS_DF_COLUMNS)))
            else:
                raise TypeError('"obs_df" needs to be a pandas DataFrame.')
        else:
            self.obs_df = obs_df  # None

        # obs_tide_correction_type:
        if isinstance(obs_tide_correction_type, str):
            if obs_tide_correction_type:
                if obs_tide_correction_type in TIDE_CORRECTION_TYPES.keys():
                    self.obs_tide_correction_type = obs_tide_correction_type
                else:
                    raise ValueError('"obs_tide_correction_type" needs to be a key in TIDE_CORRECTION_TYPES '
                                     '({})'.format(', '.join(TIDE_CORRECTION_TYPES.keys())))
            else:
                self.obs_tide_correction_type = obs_tide_correction_type  # None
        else:
            raise TypeError('"obs_tide_correction_type" needs to be a string')

        # red_tide_correction_type
        if isinstance(red_tide_correction_type, str):
            if red_tide_correction_type:
                if red_tide_correction_type in TIDE_CORRECTION_TYPES.keys():
                    self.red_tide_correction_type = red_tide_correction_type
                else:
                    raise ValueError('"red_tide_correction_type" needs to be a key in TIDE_CORRECTION_TYPES '
                                     '({})'.format(', '.join(TIDE_CORRECTION_TYPES.keys())))
            else:
                self.red_tide_correction_type = red_tide_correction_type  # None
        else:
            raise TypeError('"red_tide_correction_type" needs to be a string')

        # obs_reference_height_type:
        if isinstance(obs_reference_height_type, str):
            if obs_reference_height_type:
                if obs_reference_height_type in REFERENCE_HEIGHT_TYPE.keys():
                    self.obs_reference_height_type = obs_reference_height_type
                else:
                    raise ValueError('"obs_reference_height_type" needs to be a key in REFERENCE_HEIGHT_TYPE '
                                     '({})'.format(', '.join(REFERENCE_HEIGHT_TYPE.keys())))
            else:
                self.obs_reference_height_type = obs_reference_height_type  # None
        else:
            raise TypeError('"obs_reference_height_type" needs to be a string')

        # red_reference_height_type:
        if isinstance(red_tide_correction_type, str):
            if red_reference_height_type:
                if red_reference_height_type in REFERENCE_HEIGHT_TYPE.keys():
                    self.red_reference_height_type = red_reference_height_type
                else:
                    raise ValueError('"red_reference_height_type" needs to be a key in REFERENCE_HEIGHT_TYPE '
                                     '({})'.format(', '.join(REFERENCE_HEIGHT_TYPE.keys())))
            else:
                self.red_reference_height_type = red_reference_height_type  # None
        else:
            raise TypeError('"red_tide_correction_type" needs to be a string')

        # keep_survey
        if isinstance(keep_survey, bool):
            self.keep_survey = keep_survey
        else:
            raise TypeError('"keep_survey" needs to be a bool type.')

        # setup_df:
        if setup_df is not None:
            if isinstance(setup_df, pd.DataFrame):
                # Check if setup_df contains exactly all columns defined by self._SETUP_DF_COLUMNS:
                if all([item for item in obs_df.columns.isin(self._SETUP_DF_COLUMNS)]) and \
                        obs_df.shape[1] == len(self._SETUP_DF_COLUMNS):
                    self.setup_df = setup_df
                else:
                    raise ValueError(
                        '"setup_df" needs the following columns:{}'.format(', '.join(self._SETUP_DF_COLUMNS)))
            else:
                raise TypeError('"setup_df" needs to be a pandas DataFrame.')
        else:
            self.setup_df = setup_df  # None

    @classmethod
    def from_cg5_survey(cls, cg5_survey, keep_survey=True):
        """Constructor that generates and populates the survey object from a CG5Survey class object.

        Notes
        -----
        The observation epochs are represented by timezone aware datetime objects with TZ=<UTC>. The TZ is changed
        if necessary when loading data from any source.

        Parameters
        ----------
        keep_survey : bool, optional (default=True)
            If False, this survey is excluded from further processing.
        cg5_survey : :py:obj:`gravtools.CG5_utils.survey.CG5Survey`
            Objects of the class CG5Survey contain all data from CG-5 observation files.

        Returns
        -------
        :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        """

        # Check input arguments:
        if not isinstance(cg5_survey, CG5Survey):
            raise TypeError('"cg5_survey" has to be a CG5Survey object.')

        # Parse and check data from CG5Survey:
        # Get date from date and time:
        if cg5_survey.survey_parameters.date_time is not None:
            survey_date = cg5_survey.survey_parameters.date_time.date()
        else:
            survey_date = None

        # obs_df:
        if cg5_survey.obs_df is not None:
            # Refactor observation dataframe:
            obs_df: object = cg5_survey.obs_df.sort_values('obs_epoch').copy(
                deep=True)  # deep copy => No struggles with references

            # Add missing columns (initialized with default values):
            obs_df['g_obs_mugal'] = obs_df['g_mgal'] * 1e3
            obs_df['sd_g_obs_mugal'] = obs_df['sd_mgal'] * 1e3
            obs_df['sd_g_obs_mugal'] = obs_df['sd_mgal'] * 1e3
            obs_df['tide'] = obs_df['tide'] * 1e3
            obs_df['keep_obs'] = True

            # Check timezone of observation epoch and convert it to UTC, if necessary:
            if obs_df['obs_epoch'].dt.tz is None:  # TZ unaware => set TZ to <UTC>
                obs_df.loc[:, 'obs_epoch'] = obs_df['obs_epoch'].dt.tz_localize('UTC')
            else:
                if obs_df['obs_epoch'].dt.tz.zone is not 'UTC':  # Change TZ to <UTC>
                    obs_df.loc[:, 'obs_epoch'] = obs_df['obs_epoch'].dt.tz_convert('UTC')

            # Rename columns:
            obs_df.rename(columns={'terrain': 'corr_terrain',
                                   'tide': 'corr_tide_mugal', },
                          inplace=True)

            # Drop columns that are not needed any more:
            # - all columns that are not in _OBS_DF_COLUMNS
            obs_df = cls._obs_df_drop_columns(obs_df)

            # Add all missing columns (init as None):
            obs_df = cls._obs_df_add_columns(obs_df)

            # Change column order:
            obs_df = cls._obs_df_reorder_columns(obs_df)

            # Check, if all columns are there:
            # if not all(obs_df.columns.values == cls._OBS_DF_COLUMNS):
            if not all([item for item in obs_df.columns.isin(cls._OBS_DF_COLUMNS)]):
                raise AssertionError('Columns missing in "obs_df"')
        else:
            obs_df = None  # If no observations are available, initialize as None

        if cg5_survey.options.tide_correction is None:
            obs_tide_correction_type = 'unknown'  # e.g. "CG-5 OPTIONS" block in observation file missing
        else:
            if cg5_survey.options.tide_correction:
                obs_tide_correction_type = 'cg5_longman1959'  # built-in tide correction of the CG5
            else:
                obs_tide_correction_type = 'no_tide_corr'

        return cls(name=cg5_survey.survey_parameters.survey_name,
                   date=survey_date,
                   operator=cg5_survey.survey_parameters.operator,
                   institution=cg5_survey.survey_parameters.client,
                   gravimeter_id=DEFAULT_GRAVIMETER_ID_CG5_SURVEY,
                   data_file_name=os.path.split(cg5_survey.obs_filename)[1],  # Filename only, without path
                   data_file_type='cg5_obs_file_txt',
                   obs_df=obs_df,
                   obs_tide_correction_type=obs_tide_correction_type,
                   obs_reference_height_type='sensor_height',
                   red_tide_correction_type='',  # Not specified
                   red_reference_height_type='',  # Not specified
                   keep_survey=keep_survey,
                   )

    @classmethod
    def from_cg5_obs_file(cls, filename, keep_survey=True):
        """Constructor that generates and populates the survey object directly from a CG5 observation file.

        Parameters
        ----------
        filename : str
            Name (and path) of a CG-5 observation file (text format).
        keep_survey : bool, optional (default=True)
            If False, this survey is excluded from further processing.

        Returns
        -------
        :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        """
        cg5_survey = CG5Survey(filename)
        return cls.from_cg5_survey(cg5_survey, keep_survey=keep_survey)

    @classmethod
    def from_bev_obs_file(cls, filename, keep_survey=True, verbose=False):
        """Constructor that generates and populates the survey object from an observation file in the legacy BEV format.

        Notes
        -----
        Format of the legacy BEV observation files:
        - Line 1: Scaled values? ('Y'=True, 'N'=False)
        - Line 2: Instrument ID (e.g. 5 fpr CG5) and institution
        - Line 3: Degree of the drift polynomial to be fitted
        - Line 4: Date (YYYY MM DD)
        - Line 5: Timezone (UTC, MEZ, OEZ)
        - Line 6 bis n: Station name(x-xxx-xx), epoch (hh.mm), g [mGal], dhb [cm] (hh.h), dhf [cm] (hh.h)
        - Line n+1: 'end'

        The observation epochs are represented by timezone aware datetime objects with TZ=<UTC>. The TZ is changed
        if necessary when loading data from any source.

        Parameters
        ----------
        filename : str
            Name (and path) of an observation in the legacy BEV format.
        keep_survey : bool, optional (default=True)
            If False, this survey is excluded from further processing.
        verbose : bool, optional (default=False)
            If True, status messages are printed.

        Returns
        -------
        :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        """

        data_file_name = os.path.split(filename)[1]

        if verbose:
            print(f'Read observations from file: {filename}.')

        # Read header lines:
        num_of_header_lines = 5
        with open(filename) as myfile:
            head = [next(myfile) for x in range(num_of_header_lines)]
        scaling = head[0][:-1] == 'Y'
        gravimeter_id = head[1].split()[0]
        institution = head[1].split()[1]
        polynomial_degree = int(head[2][:-1])
        date_str = head[3][:-1]
        timezone_str = head[4][:-1]

        survey_date = dt.datetime.strptime(date_str, '%Y %m %d').date()

        # Read observations (fixed width file):
        widths = (
            11,  # Station name
            5,  # Time
            9,  # g [mGal]
            7,  # dhb [cm]
            6,  # dhf [cm]
        )
        column_names = (
            'station_name',
            'time_hh.mm',
            'g_mgal',
            'dhb_cm',
            'dhf_cm',
        )
        df = pd.read_fwf(filename, widths=widths, header=None, names=column_names, skiprows=5, skipfooter=1,
                         dtype={'time_hh.mm': object})

        # Prepare and initialize DataFrame:
        df['obs_epoch'] = pd.to_datetime(date_str + ' ' + df['time_hh.mm'] + ' ' + timezone_str,
                                         format='%Y %m %d %H.%M %Z')

        # Check timezone of observation epoch and convert it to UTC, if necessary:
        if df['obs_epoch'].dt.tz is None:  # TZ unaware => set TZ to <UTC>
            df.loc[:, 'obs_epoch'] = df['obs_epoch'].dt.tz_localize('UTC')
        else:
            if df['obs_epoch'].dt.tz.zone is not 'UTC':  # Change TZ to <UTC>
                df.loc[:, 'obs_epoch'] = df['obs_epoch'].dt.tz_convert('UTC')

        # Timestamp: https://stackoverflow.com/questions/40881876/python-pandas-convert-datetime-to-timestamp-effectively-through-dt-accessor
        df['setup_id'] = df['obs_epoch'].values.astype(np.int64) // 10 ** 9

        df['g_obs_mugal'] = df['g_mgal'] * 1e3
        df['dhb_m'] = df['dhb_cm'] * 1e-2
        df['dhf_m'] = df['dhf_cm'] * 1e-2

        df['keep_obs'] = True

        # Drop columns:
        obs_df = cls._obs_df_drop_columns(df)

        # Add all missing columns (init as None):
        obs_df = cls._obs_df_add_columns(obs_df)

        # Change column order:
        obs_df = cls._obs_df_reorder_columns(obs_df)

        # Get tide correction type:
        try:
            obs_tide_correction_type = BEV_GRAVIMETER_TIDE_CORR_LOOKUP[gravimeter_id]
        except KeyError:
            obs_tide_correction_type = 'unknown'
            if VERBOSE:
                print(f'Warning: For gravimeter ID "{gravimeter_id}" the tide correction type is unknown '
                      f'(not specified in settings.BEV_GRAVIMETER_TIDE_CORR_LOOKUP). '
                      f'It is set to "{obs_tide_correction_type}".')

        return cls(name=os.path.split(filename)[1],
                   date=survey_date,
                   operator='',
                   institution=institution,
                   gravimeter_id=gravimeter_id,
                   data_file_name=os.path.split(filename)[1],  # Filename only, without path
                   data_file_type='bev_obs_file',
                   obs_df=obs_df,
                   obs_tide_correction_type=obs_tide_correction_type,
                   obs_reference_height_type='sensor_height',
                   red_tide_correction_type='',  # Not specified
                   red_reference_height_type='',  # Not specified
                   keep_survey=keep_survey,
                   )

    @classmethod
    def _obs_df_drop_columns(cls, obs_df):
        """Drop all columns of obs_df that are not listed in cls._OBS_DF_COLUMNS"""
        columns_to_be_dropped = list(set(obs_df.columns) - set(cls._OBS_DF_COLUMNS))
        obs_df.drop(columns=columns_to_be_dropped, inplace=True)
        return obs_df

    @classmethod
    def _obs_df_add_columns(cls, obs_df):
        """Add and initialize (as None) all columns that are listed in cls._OBS_DF_COLUMNS and not present in input
        obs_df."""
        columns_to_be_initialized_as_none = list(set(cls._OBS_DF_COLUMNS) - set(obs_df.columns))
        obs_df[columns_to_be_initialized_as_none] = None
        return obs_df

    @classmethod
    def _obs_df_reorder_columns(cls, obs_df):
        """Change order of columns of obs_df to the order specified in cls._OBS_DF_COLUMNS.

        See: https://erikrood.com/Python_References/change_order_dataframe_columns_final.html
        """
        obs_df = obs_df[list(cls._OBS_DF_COLUMNS)]
        return obs_df

    def activate_setup(self, setup_id, flag_activate):
        """Activate or deactivate instrument setup with the specified ID.

        Parameters
        ----------
        setup_id : int
            ID of the setup that will be activated or deactivated.
        flag_activate : bool
            Specified whether the setup with the ID `setup_id` will be activated (`TRUE`) or deactivated (`False`).
        """
        if not isinstance(flag_activate, bool):
            raise TypeError('"flag_activate" has to be a boolean variable!')
        self.obs_df.loc[self.obs_df['setup_id'] == setup_id, 'keep_obs'] = flag_activate

    def activate_observation(self, obs_idx, flag_activate):
        """Activate or deactivate teh observation with the specified dataframe index.

        Parameters
        ----------
        obs_idx : int
            Index of the observation within the observation dataframe (`obs_df`) that will be activated or deactivated.
        flag_activate : bool
            Specified whether the observation with the index `obs_idx` will be activated (`TRUE`) or deactivated
            (`False`).
        """
        self.obs_df.at[obs_idx, 'keep_obs'] = flag_activate

    def is_valid_obs_df(self, verbose=False) -> bool:
        """Check, whether the observations DataFrame (`obs_df`) is valid.

        This method carried out the following checks:

        - Check the columns as specified in :py:obj:`.Survey._OBS_DF_COLUMNS`?

          - Does `obs_df` have all required columns?

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If True, messages are printed that indicate why `obs_df` is not valid in case.

        Returns
        -------
        bool
            True, if `obs_df` is valid.
        """
        is_valid = True

        invalid_cols = list(set(self.obs_df.columns) - set(self._OBS_DF_COLUMNS))
        if len(invalid_cols) > 0:
            is_valid = False
            if verbose:
                print(f'The following columns are not valid: {", ".join(invalid_cols)}')

        invalid_cols = list(set(self._OBS_DF_COLUMNS) - set(self.obs_df.columns))
        if len(invalid_cols) > 0:
            is_valid = False
            if verbose:
                print(f'The following columns are missing: {", ".join(invalid_cols)}')

        return is_valid

    def obs_df_drop_redundant_columns(self):
        """Drop all columns of obs_df that are not listed in self._OBS_DF_COLUMNS"""
        columns_to_be_dropped = list(set(self.obs_df.columns) - set(self._OBS_DF_COLUMNS))
        self.obs_df.drop(columns=columns_to_be_dropped, inplace=True)

    @classmethod
    def get_obs_df_column_name(cls, col_index):
        """Return the name of the observation dataframe column with the specified index.

        Parameters
        ----------
        col_index : int
            Column index for the observation dataframe.

        Returns
        -------
        str : Column name.
            Name of the column with index `col_index`
        """
        return cls._OBS_DF_COLUMNS[col_index]

    @classmethod
    def get_setup_df_column_name(cls, col_index):
        """Return the name of the setup dataframe column with the specified index.

        Parameters
        ----------
        col_index : int
            Column index for the setup dataframe.

        Returns
        -------
        str : Column name.
            Name of the column with index `col_index`
        """
        return cls._SETUP_DF_COLUMNS[col_index]

    @classmethod
    def get_obs_df_column_index(cls, column_name: str) -> int:
        """Returns the column index for specific column name for the obs_df dataframe.

        Parameters
        ----------
        column_name: str
            Name of the columns for which the column index is returned.

        Returns
        -------
        int : Column index
            Index of the column with the name `column_name`.
        """
        return cls._OBS_DF_COLUMNS.index(column_name)

    def get_number_of_observations(self) -> int:
        """Returns the number of observations of the current survey.

        Returns
        -------
        int
            Number of observations.
        """
        if self.obs_df is None:
            return 0
        else:
            return len(self.obs_df)

    def obs_df_populate_vg_from_stations(self, stations, verbose=False):
        """Populates the vertical gradient columns of the observation DataFrame with values from a Station object.

        Parameters
        ----------
        stations : :py:obj:`.Station` object
            Station data (datum- and non-datum-stations).

        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Notes
        -----
        - If an observed station is not present in the station object, or the station has no vertical gradient in the
          station object, the default vertical gradient is assigned to the observation.
        """
        # Merge stations Dataframe (`stat_df`) and observations DataFrame (`obs_df`) by the station names:
        self.obs_df = self.obs_df.merge(stations.stat_df[['station_name', 'vg_mugalm']], on='station_name',
                                        how='left', suffixes=('_x', ''))

        # Populate all missing VGs with the default value:
        self.obs_df.loc[self.obs_df['vg_mugalm'].isna(), 'vg_mugalm'] = VG_DEFAULT

        # Drop columns that are not required and check for validity:
        self.obs_df_drop_redundant_columns()
        if not self.is_valid_obs_df():
            raise AssertionError('The DataFrame "obs_df" is not valid.')

    def reduce_observations(self, target_ref_height: str = None, target_tide_corr: str = None,
                            verbose: bool = False) -> [bool, str]:
        """Reduce the observed gravity values by applying the selected corrections.

        The following corrections can be applied:

        - Reference height: Using the vertical gravity gradient at the measurement points (from station data file) and
          the vertical distances between instrument top, sensor height, ground and control point, the observations are
          reduced to the selected reference height. Information on the reference height of the input data (usually
          sensor height, as observed) and the target reference height are required.

        - Tidal correction: Observations are reduced due tidal gravitational attraction (caused by sun and moon).
          The applied model to determine the corrections for each observation, dependent on the location and time, can
          be selected. Information on the tidal corrections that were already applied on the input data is required

        Notes
        -----
        - For the reduction of the reference height vertical gravity gradients are required!

        Parameters
        ----------
        target_ref_height : string, specifying the target reference height type (default = `None`).
            The target reference height type has to be listed in :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`.
            Default is `None` indicating that the reference heights of the input data are not changed.
        target_tide_corr : str, specifying the tidal correction type to be applied (default = `None`).
            The target tidal correction type specifies what kind of tidal correction will be applied. Valid types have
            to be listed in :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`. Default is `None` indicating that the
            tidal corrections are not considered here (tidal corrections are inherited from input data).
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        Returns
        -------
        flag_corrections_applied_correctly : bool
            `False` indicates that an error occurred when applying the observation corrections.
        error_msg : str, default = ''
            Message that describes the error in case an error occurred.
        """
        # Init.:
        flag_corrections_applied = True
        error_msg = ''

        # Create a copy auf obs_df in order prevent problems with manipulation assigned py reference variables:
        obs_df = self.obs_df.copy(deep=True)

        # Initialize pandas series for reduced data by copying the observation data as loaded from the input file:
        g_red_mugal = obs_df['g_obs_mugal']
        sd_g_red_mugal = obs_df['sd_g_obs_mugal']

        # Check whether field for reduced data are initialized correctly:
        flag_all_field_are_nan = all([obs_df['g_red_mugal'].isna().all(),
                                      obs_df['sd_g_red_mugal'].isna().all(),
                                      obs_df['corr_tide_red_mugal'].isna().all()
                                      ])
        flag_no_field_is_nan = all([
            all(~obs_df['g_red_mugal'].isna().values),
            all(~obs_df['sd_g_red_mugal'].isna().values),
            all(~obs_df['corr_tide_red_mugal'].isna().values)])
        if not (flag_all_field_are_nan or flag_no_field_is_nan):
            # raise AssertionError('In "obs_df" the columns "g_red_mugal", "sd_g_red_mugal" and "corr_tide_red_mugal" '
            #                      'are not initialized consistently!'
            error_msg = 'In "obs_df" the columns "g_red_mugal", "sd_g_red_mugal" and "corr_tide_red_mugal" are not ' \
                        'initialized consistently! '
            flag_corrections_applied = False
            if verbose:
                print(error_msg)
            return flag_corrections_applied, error_msg

        if verbose:
            print(f'## Calculate reduced observations by applying the specified corrections:')

        # 1.) Check reference heights:
        if verbose:
            print(
                f'Reduction of reference heights to "{target_ref_height}" ({REFERENCE_HEIGHT_TYPE[target_ref_height]}):')

        if target_ref_height is None:  # Do nothing
            if verbose:
                print(f' - Reference heights are not changed!')
            flag_corrections_applied = True
        else:
            # Check if the target reference height type is valid:
            if target_ref_height not in REFERENCE_HEIGHT_TYPE:
                raise ValueError(f'"{target_ref_height}" is an unknown reference height type.')

            # Check if reduction is necessary:
            if self.obs_reference_height_type == target_ref_height:
                if verbose:
                    print(f' - Observations are already referenced to: {target_ref_height}')
                flag_corrections_applied = True
            else:
                # Check, if VG are available in obs_df:
                if obs_df.vg_mugalm.isna().any():
                    error_msg = f'For the following stations the vertical gravity gradient is not available: ' \
                                f'{", ".join(obs_df[obs_df.vg_mugalm.isna()].station_name.unique())}. ' \
                                f'Without vertical gradient the observed gravity cannot be reduced to another height ' \
                                f'level! '
                    flag_corrections_applied = False
                    if verbose:
                        print(error_msg)
                    return flag_corrections_applied, error_msg

                # Reduction:
                # Distance between instrument top and sensor level:
                dst_m = GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[GRAVIMETER_ID_BEV[self.gravimeter_id]]
                if self.obs_reference_height_type == 'sensor_height':
                    if target_ref_height == 'control_point':
                        # + dst_m + dhf_m
                        g_red_mugal = g_red_mugal + (dst_m + obs_df['dhf_m']) * \
                                      obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                    elif target_ref_height == 'instrument_top':
                        # + dst_m
                        g_red_mugal = g_red_mugal + dst_m * obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                    elif target_ref_height == 'ground':
                        # + dst_m + dhb_m
                        g_red_mugal = g_red_mugal + (dst_m + obs_df['dhb_m']) * \
                                      obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                elif self.obs_reference_height_type == 'instrument_top':
                    if target_ref_height == 'sensor_height':
                        # - dst_m
                        g_red_mugal = g_red_mugal - dst_m * obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                    elif target_ref_height == 'ground':
                        # + dhb_m
                        g_red_mugal = g_red_mugal + obs_df['dhb_m'] * \
                                      obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                    elif target_ref_height == 'control_point':
                        # + dhf_m
                        g_red_mugal = g_red_mugal + obs_df['dhf_m'] * \
                                      obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                elif self.obs_reference_height_type == 'ground':
                    if target_ref_height == 'instrument_top':
                        # - dhb_m - dst_m + dst_m = - dhb_m
                        g_red_mugal = g_red_mugal - obs_df['dhb_m'] * \
                                      obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                    elif target_ref_height == 'sensor_height':
                        # - dhb_m - dst_m
                        g_red_mugal = g_red_mugal - (obs_df['dhb_m'] + dst_m) * \
                                      obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                    elif target_ref_height == 'control_point':
                        # - dhb_m - dst_m + dhf_m + dst_m = dhf_m - dhb_m
                        g_red_mugal = g_red_mugal + \
                                      (obs_df['dhf_m'] - obs_df['dhb_m']) * obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                elif self.obs_reference_height_type == 'control_point':
                    if target_ref_height == 'instrument_top':
                        # - dhf_m - dst_m + dst_m = - dhf_m
                        g_red_mugal = g_red_mugal - obs_df['dhf_m'] * \
                                      obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                    elif target_ref_height == 'ground':
                        # - dhf_m - dst_m + dst_m + dhb_m = dhb_m - dhf_m
                        g_red_mugal = g_red_mugal + \
                                      (obs_df['dhb_m'] - obs_df['dhf_m']) * obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                    elif target_ref_height == 'sensor_height':
                        # - dhf_m - dst_m
                        g_red_mugal = g_red_mugal - (obs_df['dhf_m'] + dst_m) * \
                                      obs_df['vg_mugalm']
                        # sd_g_red_mugal = sd_g_red_mugal  # Keep SD
                if verbose:
                    print('...done!')
                flag_corrections_applied = True

        # 2.) Check tidal corrections:
        if verbose:
            print(f'Tidal corrections "{target_tide_corr}" ({TIDE_CORRECTION_TYPES[target_tide_corr]}):')

        if target_tide_corr is None:  # Do nothing
            if verbose:
                print(f' - Tidal corrections are not changed!')
            corr_tide_red_mugal = obs_df['corr_tide_mugal']
            flag_corrections_applied = True
        else:
            # Check if the target tidal correction type is valid:
            if target_tide_corr not in TIDE_CORRECTION_TYPES:
                # raise ValueError(f'"{target_tide_corr}" is unknown ab invalid!')
                error_msg = f'"{target_tide_corr}" is unknown ab invalid!'
                flag_corrections_applied = False
                if verbose:
                    print(error_msg)
                return flag_corrections_applied, error_msg

            # Check if reduction is necessary or possible:
            if self.obs_tide_correction_type == target_tide_corr:
                if verbose:
                    print(f' - Observations are already reduced by: {target_tide_corr}')
                corr_tide_red_mugal = obs_df['corr_tide_mugal']
                flag_corrections_applied = True

            elif self.obs_tide_correction_type == 'unknown':
                error_msg = 'Unknown tidal corrections at input data!'
                flag_corrections_applied = False
                if verbose:
                    print(error_msg)
                return flag_corrections_applied, error_msg
            else:
                # Reduction:
                if self.obs_tide_correction_type == 'no_tide_corr':  # status of unreduced observations
                    if target_tide_corr == 'cg5_longman1959':  # Add instrumental corrections
                        g_red_mugal = g_red_mugal + obs_df['corr_tide_mugal']
                        # sd_g_red_mugal = sd_g_red_mugal # Keep SD
                        corr_tide_red_mugal = obs_df['corr_tide_mugal']
                    # elif target_tide_corr == '......':
                elif self.obs_tide_correction_type == 'cg5_longman1959':
                    if target_tide_corr == 'no_tide_corr':  # Subtract instrumental corrections
                        g_red_mugal = g_red_mugal - obs_df['corr_tide_mugal']
                        # sd_g_red_mugal = sd_g_red_mugal # Keep SD
                        corr_tide_red_mugal = obs_df['corr_tide_mugal']
                        corr_tide_red_mugal.values[:] = 0
                    # elif target_tide_corr == '......':
                if verbose:
                    print('...done!')
                flag_corrections_applied = True

        # Apply changes, if everything worked out in here:
        if flag_corrections_applied:
            self.red_tide_correction_type = target_tide_corr
            self.red_reference_height_type = target_ref_height
            self.obs_df['g_red_mugal'] = g_red_mugal
            self.obs_df['sd_g_red_mugal'] = sd_g_red_mugal
            self.obs_df['corr_tide_red_mugal'] = corr_tide_red_mugal
        return flag_corrections_applied, error_msg

    def autselect_tilt(self, threshold_arcsec: int, setup_id: int = None):
        """Deactivate all observations of the survey or of a setup with a tilt larger than the defined threshold.

        Parameters
        ----------
        threshold_arcsec : int
            Observations in this survey or the specified setup are deactivated, if their tilt in X or Y direction
            (columns `tiltx` and `tilty` in :py:obj:`.Survey.obs_df`) exceeds the given threshold [arcsec].
        setup_id : int (default=None)
            `None` implies that this autoselection function is applied on all observations of this survey. Otherwise,
            the autoselection function is only applied on observations of the setup wirth the provided ID (`setup_id`)
        """
        filter_tilt = (abs(self.obs_df['tiltx']) > threshold_arcsec) | (abs(self.obs_df['tilty']) > threshold_arcsec)
        if setup_id is not None:  # Apply on whole survey
            filter_tilt = filter_tilt & (self.obs_df['setup_id'] == setup_id)
        self.obs_df.loc[filter_tilt, 'keep_obs'] = False

    def autselect_g_sd(self, threshold_mugal: int, obs_type: str = 'reduced', setup_id: int = None,
                       verbose: bool = False):
        """Deactivate observations with a gravity standard deviation larger than the defined threshold.

        Parameters
        ----------
        threshold_mugal : int
            Observations in this survey or the specified setup are deactivated, if the standard deviation of the
            observed or the reduced gravity (columns `sd_g_obs_mugal` or `sd_g_red_mugal` in
            :py:obj:`.Survey.obs_df`) exceed the given threshold [µGal].
        obs_type : str, 'observed' or 'reduced' (default)
            Defines whether the observed (as loaded from an observation file) or the reduced observations are used as
            reference for this autoselection function.
        setup_id : int (default=None)
            `None` implies that this autoselection function is applied on all observations of this survey. Otherwise,
            the autoselection function is only applied on observations of the setup wirth the provided ID (`setup_id`).
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        Returns
        -------
        bool : Error indicator.
            `False`, if at least one standard deviation value is `None` is case, reduced observations are usd as
            reference.
        """
        if obs_type == 'reduced':
            if self.obs_df['sd_g_red_mugal'].isna().any():
                if verbose:
                    print('ERROR: At least one value in column "sd_g_red_mugal" is None.')
                return False
            filter_g_sd = (abs(self.obs_df['sd_g_red_mugal']) > threshold_mugal)
        elif obs_type == 'observed':
            filter_g_sd = (abs(self.obs_df['sd_g_obs_mugal']) > threshold_mugal)
        else:
            raise ValueError(f'Invalid value assigned to the parameter "obs_type": {obs_type}')
        if setup_id is not None:  # Apply on whole survey
            filter_g_sd = filter_g_sd & (self.obs_df['setup_id'] == setup_id)
        self.obs_df.loc[filter_g_sd, 'keep_obs'] = False
        return True

    def autselect_delta_g(self, threshold_mugal: int, n_obs: int = 3, obs_type: str = 'reduced', setup_id: int = None,
                          verbose: bool = False):
        """Deactivate observations based ob the deviation of the stabilized gravity at a setup.

        Notes
        -----
        If a setup consists of less than `n_obs` observations, this autosection function is not applied.


        Parameters
        ----------
        threshold_mugal : int
            Observations in this survey or the specified setup are deactivated, if the observed or the reduced gravity
            (columns `sd_g_obs_mugal` or `sd_g_red_mugal` in :py:obj:`.Survey.obs_df`) deviates from the stabilized
            gravity at the setup by more than the given threshold. The stabilized gravity at a setup is calculated as
            the mean (observed or reduced) gravity of the last `n_obs` observations in a setup.
        n_obs : int, optional (default=3)
            Number of observations at the end of a setup that are used to calculate the stabilized gravity as reference
            for deactivating observations.
        obs_type : str, 'observed' or 'reduced' (default)
            Defines whether the observed (as loaded from an observation file) or the reduced observations are used as
            reference for this autoselection function.
        setup_id : int, optional (default=None)
            `None` implies that this autoselection function is applied on all setup of this survey. Otherwise,
            the autoselection function is only applied on observations of the setup with the provided ID (`setup_id`).
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        Returns
        -------
        bool : Error indicator.
            `False`, if at least one standard deviation value is `None` is case, reduced observations are usd as
            reference.
        """
        # TODO: Check, if it makes sense to keept the last n_obs observations anyway (even if they do not meet the
        #  conditions here)! Check if reduced observations are available, if required:
        if obs_type == 'reduced':
            if self.obs_df['g_red_mugal'].isna().any():
                if verbose:
                    print('ERROR: At least one value in column "g_red_mugal" is None.')
                return False

        # Get list of IDs of all setups that will be handled:
        if setup_id is None:
            setup_ids = self.get_setup_ids()
        else:
            setup_ids = [setup_id]

        # Loop over all setups:
        for id in setup_ids:
            filter_id = self.obs_df['setup_id'] == id
            if verbose:
                print(f' - setup ID: {id}')
            # Check number of observations in setup:
            if len(self.obs_df.loc[filter_id]) < (n_obs + 1):
                if verbose:
                    print(f'   - Less than {n_obs} observations available: {len(self.obs_df.loc[filter_id])}')
            else:  # Enough observations in this setup
                # Calculate stabilized gravity and set up filter:
                if obs_type == 'reduced':
                    g_stabilized_mugal = self.obs_df.loc[filter_id, 'g_red_mugal'].tail(n_obs).mean()
                    filter_delta_g = (self.obs_df['g_red_mugal'] < (g_stabilized_mugal - threshold_mugal)) | (
                            self.obs_df['g_red_mugal'] > (g_stabilized_mugal + threshold_mugal))
                elif obs_type == 'observed':
                    g_stabilized_mugal = self.obs_df.loc[filter_id, 'g_obs_mugal'].tail(n_obs).mean()
                    filter_delta_g = (self.obs_df['g_obs_mugal'] < (g_stabilized_mugal - threshold_mugal)) | (
                            self.obs_df['g_obs_mugal'] > (g_stabilized_mugal + threshold_mugal))
                else:
                    raise ValueError(f'Invalid value assigned to the parameter "obs_type": {obs_type}')

                # Apply filter:
                filter_all = filter_delta_g & filter_id
                if verbose:
                    print(f'   - Number of removed observations: {len(self.obs_df.loc[filter_all, "keep_obs"])}')
                self.obs_df.loc[filter_all, 'keep_obs'] = False
        return True

    def get_setup_ids(self) -> list:
        """Return the IDs of all setups in this survey

        Returns
        list : List of all setup IDs.
        """
        return self.obs_df.setup_id.unique().tolist()

    def __str__(self):
        if self.obs_df is not None:
            return f'Survey "{self.name}" with {len(self.obs_df)} observations.'
        else:
            return f'Survey "{self.name}"'

    def reset_setup_data(self, verbose=False):
        """Deletes the setup data.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        if self.setup_df is None:
            if verbose:
                print('Nothing to delete. Setup data is empty.')
        else:
            if verbose:
                print('Setup data deleted.')
            self.setup_df = None

    def calculate_setup_data(self, obs_type='reduced', verbose=False):
        """Accumulate all active observation within each setup and calculate a single representative pseudo observation.

        Parameters
        ----------
        obs_type : str, 'observed' or 'reduced' (default)
            Defines whether the observed (as loaded from an observation file) or the reduced observations from
            `self.obs_df` are used to determine the weighted mean values per setup.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """

        _VALID_OBS_TYPES = ('observed', 'reduced',)

        self.reset_setup_data(verbose)

        # Initial checks:
        if obs_type not in _VALID_OBS_TYPES:
            raise AssertionError(f'Invalid observation type: {obs_type} (valid: "reduced" or "observed").')
        if self.obs_df is None:
            raise AssertionError('Observation dataframe is empty!')

        # Get all active observations:
        tmp_filter = self.obs_df['keep_obs'] == True
        active_obs_df = self.obs_df[tmp_filter].copy(deep=True)

        # Check, if at least one observation is active:
        if len(active_obs_df) == 0:
            # raise AssertionError(f'No active observations in survey {self.name}')
            if verbose:
                print(f'No active observations in survey {self.name}')
        else:

            # Check, if reduced data is available:
            if obs_type == 'reduced':
                if active_obs_df['g_red_mugal'].isnull().any() or active_obs_df['sd_g_red_mugal'].isnull().any():
                    raise AssertionError('Reduced observations (g and sd) are missing!')

            # Initialize colums lists for creating dataframe:
            station_name_list = []
            setup_id_list = []
            g_mugal_list = []
            sd_g_red_mugal_list = []
            obs_epoch_list_unix = []
            obs_epoch_list_dt = []

            # Loop over setups:
            setup_ids = active_obs_df['setup_id'].unique()
            for setup_id in setup_ids:
                tmp_filter = active_obs_df['setup_id'] == setup_id
                if obs_type == 'observed':
                    g_mugal = active_obs_df.loc[tmp_filter, 'g_obs_mugal'].to_numpy()
                    sd_g_mugal = active_obs_df.loc[tmp_filter, 'sd_g_obs_mugal'].to_numpy()
                elif obs_type == 'reduced':
                    g_mugal = active_obs_df.loc[tmp_filter, 'g_red_mugal'].to_numpy()
                    sd_g_mugal = active_obs_df.loc[tmp_filter, 'sd_g_red_mugal'].to_numpy()
                weights = 1 / sd_g_mugal ** 2
                g_setup_mugal = np.sum(g_mugal * weights) / np.sum(weights)
                sd_g_setup_mugal = np.sqrt(1 / np.sum(weights))
                if len(active_obs_df.loc[tmp_filter, 'station_name'].unique()) > 1:
                    raise AssertionError(
                        f'Setup with ID "{setup_id}" in survey "{self.name}" contains '
                        f'{len(active_obs_df.loc[tmp_filter, "station_name"].unique())} stations (only 1 allowed)!'
                    )
                # observation epoch (UNIX timestamps in full seconds):
                obs_epochs_series = active_obs_df.loc[tmp_filter, 'obs_epoch']
                unix_obs_epochs = obs_epochs_series.values.astype(np.int64) // 10 ** 9
                unix_setup_epoch = np.sum(unix_obs_epochs * weights) / np.sum(weights)

                # Reference epoch as datetime object (TZ=<UTC>):
                dt_setup_epoch = dt.datetime.fromtimestamp(unix_setup_epoch)
                dt_setup_epoch = dt_setup_epoch.replace(tzinfo=dt.timezone.utc)  # TZ = <UTC>

                station_name_list.append(active_obs_df.loc[tmp_filter, 'station_name'].unique()[0])
                setup_id_list.append(setup_id)
                g_mugal_list.append(g_setup_mugal)
                sd_g_red_mugal_list.append(sd_g_setup_mugal)
                obs_epoch_list_unix.append(unix_setup_epoch)
                obs_epoch_list_dt.append(dt_setup_epoch)

            # convert to pd dataframe:
            self.setup_df = pd.DataFrame(list(zip(station_name_list, setup_id_list, g_mugal_list, sd_g_red_mugal_list,
                                                  obs_epoch_list_unix, obs_epoch_list_dt)),
                                         columns=self._SETUP_DF_COLUMNS)


class Campaign:
    """Gravity Campaign dataset.

    A gravity campaign datasets consists of:

    - Campaign name
    - One or more gravity surveys that belong together and
      - Each survey was observed with one gravimeter on a single day
    - Station data (datum and non-datum stations)
    - Reductions and corrections
      - All observations (from surveys) are corrected and reduced in the same way

    Attributes
    ----------
    campaign_name : str
        Name of the campaign.
    output_directory : str
        Path to output directory (all output files are stored there).
    surveys: dict of :py:obj:`.Survey` objects
        Arbitrary number of survey objects.
    stations : :py:obj:`.Station` object
        Data of known stations (datum- and non-datum-stations).
    lsm_runs : list of objects inherited from :py:obj:`gravtools.models.lsm.LSM`
        Each item in the list contains one enclosed LSM object. Each LSM object reflects one dedicated run of an
        least-squares adjustment in order to estimate target parameters.
    """

    def __init__(self,
                 campaign_name,
                 output_directory,
                 surveys=None,  # Always use non-mutable default arguments!
                 stations=None,  # Always use non-mutable default arguments!
                 lsm_runs=None  # Always use non-mutable default arguments!
                 ):
        """
        Parameters
        ----------
        campaign_name : str
            Name of the campaign.
        surveys: dict of :py:obj:`.Survey` objects, optional
            Arbitrary number of survey data objects. Default=None which implies that the campaign will be initialized
            without surveys.
        stations: :py:obj:`.Station` object, optional
            Station data (datum- and non-datum-stations). Default=None implies that the campaign will be
            initialized without station data.
        lsm_runs : list of objects inherited from :py:obj:`gravtools.models.lsm.LSM`
            Each item in the list contains one enclosed LSM object. Each LSM object reflects one dedicated run of an
            least-squares adjustment in order to estimate target parameters.

        Raises
        ------
        TypeError
            Wrong input argument type.
        """

        # Check campaign_name:
        if not isinstance(campaign_name, str):
            raise TypeError('The argument "campaign_name" needs to be a string.')
        else:
            if not campaign_name:
                raise ValueError('"campaign_name" should not be empty!')
        self.campaign_name = campaign_name

        # Check output directory:
        if not isinstance(output_directory, str):
            raise TypeError('The argument "output_directory" needs to be a string.')
        else:
            if not output_directory:
                raise ValueError('"output_directory" should not be empty!')
        self.output_directory = output_directory

        # Check surveys:
        if surveys is None:
            surveys = {}
        else:
            if not isinstance(surveys, dict):
                raise TypeError('The argument "survey" needs to be a dict of Survey objects.')
            else:
                for survey_name, survey_obj in surveys.items():
                    if not isinstance(survey_name, str):
                        raise TypeError('The argument "survey" needs to be a string.')
        self.surveys = surveys  # dict: key=Name of Survey, value=Survey object

        # Check stations:
        if stations is None:
            stations = Station()
        else:
            if not isinstance(stations, Station):
                raise TypeError('The argument "stations" needs to be a Station object.')
        self.stations = stations

        # Check lsm_runs:
        if lsm_runs is None:
            lsm_runs = []  # Empty list
        else:
            if not isinstance(lsm_runs, list):
                raise TypeError('The argument "lsm_runs" needs to be a list of LSM-objects.')
            else:
                for items in lsm_runs:
                    if not isinstance(lsm_runs, LSM):
                        raise TypeError('The argument "lsm_runs" needs to be a list of LSM-objects.')
        self.lsm_runs = lsm_runs

    def add_survey(self, survey_add: Survey, verbose=False) -> bool:
        """Add a survey to campaign and specify whether to use it for ths analysis.

        Notes
        -----
        A survey can only be added, f the survey's name is unique within the campaign.

        Parameters
        ----------
        survey_add : :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        bool
            True, if the survey was successfully added; False, if not.
        """
        # Check if a survey with the dame name ("survey_add.name") already exists in this campaign:
        # - Raise warning:
        if survey_add.name in self.surveys.keys():
            if verbose:
                print(f'Warnung: the current campaign already contains a survey named {survey_add.name}.')
                print(' - Survey names need to be unique within a campaign.')
            return False
        else:
            # Add survey:
            self.surveys[survey_add.name] = survey_add
            if verbose:
                print(f"Survey {survey_add.name} added to the campaign.")
            return True

    def remove_survey(self, survey_name: str, verbose=False) -> bool:
        """Remove survey with the specified name from the campaign

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be removed from the campaign.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        bool
            True, if the survey was successfully removed; False, if not.
        """
        try:
            del self.surveys[survey_name]
            if verbose:
                print(f'Survey "{survey_name}" removed from campaign.')
        except KeyError:
            if verbose:
                print(f'Survey "{survey_name}" does not exist.')
            return False
        except:
            if verbose:
                print(f'Failed to remove survey {survey_name}.')
            return False
        else:
            return True

    def activate_survey(self, survey_name: str, verbose=False) -> bool:
        """Set the survey with the specified name active.

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be set active.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        bool
            True, if the survey was successfully activated; False, if not.
        """
        try:
            if self.surveys[survey_name].keep_survey:
                if verbose:
                    print(f'Survey "{survey_name}" already active.')
                return True
            else:
                self.surveys[survey_name].keep_survey = True
                if verbose:
                    print(f'Survey "{survey_name}" activated.')
        except KeyError:
            if verbose:
                print(f'Survey "{survey_name}" does not exist.')
            return False
        except:
            if verbose:
                print(f'Failed to activate survey "{survey_name}"')
            return False
        else:
            return True

    def deactivate_survey(self, survey_name: str, verbose=False) -> bool:
        """Set the survey with the specified name inactive.

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be set inactive.
        verbose : bool, optional (default=False)
                If True, status messages are printed to the command line.

        Returns
        -------
        bool
            True, if the survey was successfully deactivated; False, if not.
        """
        try:
            if not self.surveys[survey_name].keep_survey:
                if verbose:
                    print(f'Survey "{survey_name}" already inactive.')
                return True
            else:
                self.surveys[survey_name].keep_survey = False
                if verbose:
                    print(f'Survey "{survey_name}" deactivated.')
        except KeyError:
            if verbose:
                print(f'Survey "{survey_name}" does not exist.')
            return False
        except:
            if verbose:
                print(f'Failed to deactivate survey "{survey_name}"')
            return False
        else:
            return True

    def get_survey_names_and_status(self, verbose: bool = False) -> dict:
        """Return list with all survey names and information whether the survey is set active.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If True, survey names and status are printed to the command line.

        Returns
        -------
        dict
            The keys are the survey names and the values represent the respective status (active=True, inactive=False).
        """
        if verbose:
            print('Surveys and their status:')
        info_dict = {}
        lookup_dict = {True: 'active', False: 'inactive'}
        for surv_name, surv_obj in self.surveys.items():
            if verbose:
                activity_str = lookup_dict[surv_obj.keep_survey]
                print(f' - {surv_name:12s} ({activity_str:8s}): {surv_obj.get_number_of_observations()} observations')
            info_dict[surv_name] = surv_obj.keep_survey
        return info_dict

    @property
    def number_of_surveys(self) -> int:
        """int : Returns the number of surveys in this campaign."""
        return len(self.surveys)

    @property
    def number_of_stations(self) -> int:
        """int : Returns the number of stations in this campaign."""
        return self.stations.get_number_of_stations

    def reduce_observations_in_all_surveys(self, target_ref_height=None, target_tide_corr=None, verbose=False):
        """Reduce the observed gravity by applying the specified corrections.

        Notes
        -----
        - For this reduction vertical gravity gradients are required. They are obtained from the `Station` object.
          Hence, a `Station` object has to be attached to the Campaign object beforehand.

        - All corrections are applied on the survey-level. See py:obj:`.Survey.reduce_observations` for more details.

        Parameters
        ----------
        target_ref_height : string, specifying the target reference height type (default = `None`).
            The target reference height type has to be listed in :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`.
            Default is `None` indicating that the reference heights of the input data are not changed.
        target_tide_corr : str, specifying the tidal correction type to be applied (default = `None`).
            The target tidal correction type specifies what kind of tidal correction will be applied. Valid types have
            to be listed in :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`. Default is `None` indicating that the
            tidal corrections are not considered here (tidal corrections are inherited from input data).
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        flag_corrections_applied_correctly : bool, default = True
            `False` indicates that an error occurred when applying the observation corrections.
        error_msg : str, default = ''
            Message that describes the error in case an error occurred.

        """
        flag_corrections_applied_correctly = True
        error_msg = ''
        if verbose:
            print(f'## Reduce all observation in this campaign:')
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f'Survey {survey_name}:')
            if target_ref_height is not None:
                if verbose:
                    print(f' - Get vertical gradients')
                survey.obs_df_populate_vg_from_stations(self.stations, verbose=verbose)
            flag_corrections_applied_correctly, error_msg = survey.reduce_observations(
                target_ref_height=target_ref_height,
                target_tide_corr=target_tide_corr,
                verbose=verbose)
            # Check, if corrections were applied correctly and return error message if an error occurred:
            if not flag_corrections_applied_correctly:
                error_msg = f'Survey: {survey_name}: ' + error_msg
                if verbose:
                    print(error_msg)
                return flag_corrections_applied_correctly, error_msg
                # raise AssertionError('Reduction to reference height failed!')
        return flag_corrections_applied_correctly, error_msg

    def add_stations_from_oesgn_table_file(self, oesgn_filename, verbose=False):
        """Add station from an OESGN table file.

        Parameters
        ----------
        oesgn_filename : string, specifying the path/file of the OESGN file
            Stations in the specified OESGN table file are added to the Campaign.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.
        """
        self.stations.add_stations_from_oesgn_table(filename=oesgn_filename, verbose=verbose)

    def __str__(self):
        return f'Campaign "{self.campaign_name}" with {self.number_of_surveys} surveys ' \
               f'and {self.stations.get_number_of_stations} stations.'

    def synchronize_stations_and_surveys(self, verbose=False):
        """Synchronize information between station and survey data in the campaign.

        The following information is synchronized:

        - The `is_observed` flags in the :py:obj:`.Campaign.stations.stat_df` are set according to the surveys in
          :py:obj:`.Campaign.surveys`. `True` indicated that the station as observed at least once.
        - Populates the vertical gradient columns (``) of the observation DataFrames (:py:obj:`.Campaign.surveys`) with
          values from a Station object (:py:obj:`.Campaign.stations`).
        - Add observed stations

        Notes
        -----
        It is recommended to run this method whenever new stations and/or new surveys are added to the campaign on order
        to synchronize the data.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        self.stations.stat_df['is_observed'] = False  # Reset to default.
        self.stations.stat_df['in_survey'] = None  # Reset to default.
        self.sync_observed_stations(verbose=verbose)  # Add stations from surveys.

        # Loop over all surveys to match and synchronize the survey data with the station data:
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f' - Survey: {survey_name}')
            survey.obs_df_populate_vg_from_stations(self.stations, verbose=verbose)
            self.stations.set_observed_info_from_survey(survey, verbose=verbose)

    def sync_observed_stations(self, verbose=False):
        """Adds all stations that were observed in at least one survey to this campaign's station dataframe.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Loop over surveys in this campaign:
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f' - Survey: {survey_name}')
            self.stations.add_stations_from_survey(survey, verbose)

    def calculate_setup_data(self, obs_type='reduced', verbose=False):
        """Calculate accumulated pseudo observations for each active setup in all active surveys.

        Parameters
        ----------
        obs_type : str, 'observed' or 'reduced' (default)
            Defines whether the observed (as loaded from an observation file) or the reduced observations from
            `self.obs_df` are used to determine the weighted mean values per setup.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Loop over all surveys in the campaign:
        if verbose:
            print(f'Calculate setup data:')
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f' - Survey: {survey_name}')
            if survey.keep_survey:
                survey.calculate_setup_data(obs_type=obs_type, verbose=verbose)
        pass

    def initialize_and_add_lsm_run(self, lsm_method, comment=''):
        """Initialize and add an least-squares adjustment run (object) to the campaign.

        Parameters
        ----------
        lsm_method : str
            Defines the adjustment method. Has to b listed in :py:obj:`gravtools.settings.ADJUSTMENT_METHODS`.
        comment : str, optional (default = '')
            Optional comment on the adjustment run.
        """
        # Initialize LSM object:
        if lsm_method == 'LSM_diff':
            lsm_run = LSMDiff.from_campaign(self, comment)
        else:
            raise AssertionError('Unknown LSM method: {lsm_method}')
        # Add LSM object to campaign:
        self.lsm_runs.append(lsm_run)


if __name__ == '__main__':
    """Main function, primarily for debugging and testing."""

    from gravtools.settings import NAME_OESGN_TABLE, PATH_OESGN_TABLE, VERBOSE, PATH_OBS_FILE_CG5, NAME_OBS_FILE_CG5

    station_files_dict = {PATH_OESGN_TABLE + NAME_OESGN_TABLE: 'oesgn_table'}
    stat = Station(station_files_dict)
    stat.add_stations_from_oesgn_table(PATH_OESGN_TABLE + NAME_OESGN_TABLE, verbose=VERBOSE)

    # delete stations:
    stat.delete_station(['2-008-02', '2-001-00'], verbose=VERBOSE)

    # ### Getters ###
    # Get ÖSGN stations:
    # print(stat.get_oesgn_stations)

    # Get all stations:
    print(stat.get_all_stations)

    cg5surv = CG5Survey(PATH_OBS_FILE_CG5 + NAME_OBS_FILE_CG5)

    surv = Survey.from_cg5_survey(cg5surv)
    print(surv)

    surv.activate_setup(1592449837, False)
    surv.activate_observation(18, False)

    surv2 = Survey.from_cg5_obs_file(PATH_OBS_FILE_CG5 + NAME_OBS_FILE_CG5)
    print(surv2)

    camp = Campaign('AD1', 'dir_name', stations=stat, surveys={surv2.name: surv2})

    surv_bev = Survey.from_bev_obs_file(PATH_OBS_FILE_BEV + NAME_OBS_FILE_BEV, verbose=VERBOSE)

    validity = surv_bev.is_valid_obs_df(verbose=True)

    camp.add_survey(survey_add=surv_bev, verbose=VERBOSE)

    surveys_info_dict = camp.get_survey_names_and_status(verbose=VERBOSE)

    surv.obs_df_populate_vg_from_stations(stat, verbose=True)

    print(surv.reduce_observations(target_ref_height='control_point', target_tide_corr='no_tide_corr', verbose=True))
    print(surv.reduce_observations(target_ref_height='control_point', target_tide_corr='cg5_longman1959', verbose=True))
    print(surv.reduce_observations(target_ref_height='ground', target_tide_corr='no_tide_corr', verbose=True))
    print(surv.reduce_observations(target_ref_height='instrument_top', target_tide_corr='cg5_longman1959',
                                   verbose=True))
    print(surv.reduce_observations(target_ref_height='sensor_height', target_tide_corr='no_tide_corr', verbose=True))

    print(camp.reduce_observations_in_all_surveys(target_ref_height='control_point', target_tide_corr='cg5_longman1959',
                                                  verbose=True))

    camp.synchronize_stations_and_surveys()

    surv.autselect_tilt(threshold_arcsec=5, setup_id=None)
    # surv.autselect_g_sd(threshold_mugal=10, obs_type='observed', setup_id=None, verbose=True)
    # surv.autselect_delta_g(threshold_mugal=10, obs_type='observed', setup_id=None, verbose=True)
    # surv.autselect_delta_g(threshold_mugal=10, obs_type='observed', setup_id=1599551889, verbose=True)

    surv.reset_setup_data(verbose=True)
    surv.calculate_setup_data(obs_type='observed', verbose=True)

    camp.deactivate_survey('n20200701_1')  # From BEV obsrvation file => SD of observations is missing!
    camp.calculate_setup_data(obs_type='observed', verbose=True)

    # Test adjustment:
    from gravtools.models import lsm

    lsm_diff = lsm.LSMDiff.from_campaign(camp)
    lsm_diff.adjust(drift_pol_degree=1, sig0_mugal=10, scaling_factor_datum_observations=1e-3, verbose=True)

    camp.initialize_and_add_lsm_run('LSM_diff', 'Test number one!')

    pass
