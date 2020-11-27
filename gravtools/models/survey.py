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
        - is_oesgn : bool
            Flag, that indicates whether the station is part auf the Austrian gravity base network (ÖSGN). If `True`,
            the station is part of the reference network.
    """

    _STAT_DF_COLUMNS = (
        'station_name',  # Station name, str
        'long_deg',  # longitude [deg], float
        'lat_deg',  # latitude [deg], float
        'height_m',  # Height [m]
        'g_mugal',  # gravity [µGal]
        'sd_g_mugal',  # standard deviation of the gravity [µGal]
        'vg_mugalm',  # vertical gradient [µGal/m]
        'is_oesgn',  # flag: True, if station is a OESGN station.
    )

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
        stat_df_oesgn['is_oesgn'] = True

        stat_df_oesgn.drop(columns=columns_to_be_dropped, inplace=True)

        return stat_df_oesgn

    def add_stations_from_oesgn_table(self, filename, verbose=False):
        """Adds (appends) all stations of the input OESGN table to the station dataframe.

        Parameters
        ----------
        filename : str
            Name (and path) of the OESGN table file.
        verbose : bool, optional
            Print notifications if True (default=False)
        """
        stat_df_new = self._read_oesgn_table(filename)
        number_of_existing_stations = len(self.stat_df)

        # Concatinate multiple df without duplicate rows and a unique indices:
        # https://stackoverflow.com/questions/21317384/pandas-python-how-to-concatenate-two-dataframes-without-duplicates
        # https://learndataanalysis.org/concatenate-pandas-dataframes-without-duplicates/
        # - pd.concat([self.stat_df, stat_df_new]) => simply concat both df, while it may result in duplicate indices
        #   and rows.
        # - drop_duplicates() => Drop duplicate rows
        # - reset_index(drop=True) => Reset index creates a new indec col with unique indices. With drop=True the old
        #   index column is dropped.

        self.stat_df = pd.concat([self.stat_df, stat_df_new]).drop_duplicates().reset_index(drop=True)

        if verbose:
            number_of_new_stations = len(stat_df_new)
            number_of_stations = len(self.stat_df)
            stations_added = number_of_stations - number_of_existing_stations
            print(f"{number_of_new_stations} stations loaded from {filename}. "
                  f"{stations_added} stations added ({number_of_new_stations - stations_added} already existed).")

    def delete_station(self, station_names, verbose=False):
        """
        Deletes all rows with the observed station's name occurs in the `station_names` list.

        Parameters
        ----------
        station_names : list of str
            List of station names. The related station records will be deleted from `self.stat_df`.
        verbose : bool
            If True, print notification on deleted items.
        """
        idx_series = self.stat_df.station_name.isin(station_names)
        if verbose:
            deleted_df = self.stat_df[idx_series]
            print(f'Deleted {len(deleted_df)} rows in station dataframe:')
            for index, row in deleted_df.iterrows():
                print(f' - Station {row["station_name"]:10s} (row index {index:d})')
        self.stat_df = self.stat_df[~idx_series]

    @property
    def get_number_of_stations(self):
        """int : Returns the number of stations."""
        return len(self.stat_df)

    def __str__(self):
        return f"Station DataFrame containing {self.get_number_of_stations} stations."

    @property
    def get_all_stations(self):
        """:py:obj:`pandas.core.frame.DataFrame` : Returns all available stations."""
        return self.stat_df

    @property
    def get_oesgn_stations(self):
        """:py:obj:`pandas.core.frame.DataFrame` : Returns all OESGN stations."""
        return self.stat_df[self.stat_df['is_oesgn']]


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
        Name of the survy. This parameter is mandatory!
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
        - corr_tide : float, optional (default=None)
            Tidal correction [µGal] as determined by the built-in model of the Scintrex CG-5 (Longman, 1959). Be aware
            that the tidal corrections by the CG-5 model is determined fot the middl of the observation (also see
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
    """

    _OBS_DF_COLUMNS = (
        'station_name',  # Name of station (str)
        'setup_id',  # unique ID of setup (int)
        'loop_id',  # Line ID, optional (default=None) (int)
        'lon_deg',  # Longitude [deg], optional (float)
        'lat_deg',  # Latitude [deg], optional (float)
        'alt_m',  # Altitude [m], optional (float)
        'obs_epoch',  # Observation epoch (datetime obj)
        'g_obs_mugal',  # observed g from obs file [µGal] (float)
        'sd_g_obs_mugal',  # Standard deviation of g observed from obs file (float) [µGal]
        'g_red_mugal',  # Reduced gravity observation at station (float) [µGal]
        'sd_g_red_mugal',  # Standard deviation of the reduced gravity (float) [µGal]
        'corr_terrain',  # Terrain correction [??]
        'corr_tide',  # Tidal correction [mGal], optional
        'temp',  # Temperature [mK], optional
        'tiltx',  # [arcsec], optional
        'tilty',  # [arcsec], optional
        'dhf_m',  # Distance between instrument top and physical reference point (float) [m]
        'dhb_m',  # Distance between instrument top and ground (float) [m]
        'keep_obs',  # Remove observation, if false (bool)
        'vg_mugalm',  # vertical gradient [µGal/m]
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

    @classmethod
    def from_cg5_survey(cls, cg5_survey, keep_survey=True):
        """Constructor that generates and populates the survey object from a CG5Survey class object.

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
            obs_df: object = cg5_survey.obs_df.copy(deep=True)  # deep copy => No struggles with references

            # Add missing columns (initialized with default values):
            obs_df['g_obs_mugal'] = obs_df['g_mgal'] * 1e3
            obs_df['sd_g_obs_mugal'] = obs_df['sd_mgal'] * 1e3
            obs_df['keep_obs'] = True

            # Rename columns:
            obs_df.rename(columns={'terrain': 'corr_terrain',
                                   'tide': 'corr_tide',},
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
        # Timestamp: https://stackoverflow.com/questions/40881876/python-pandas-convert-datetime-to-timestamp-effectively-through-dt-accessor
        df['setup_id'] = df['obs_epoch'].values.astype(np.int64) // 10 ** 9

        df['g_obs_mugal'] = df['g_mgal'] * 1e3
        df['dhb_m'] = df['dhb_cm'] * 1e-2
        df['dhf_m'] = df['dhf_cm'] * 1e-2

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
            Data of known stations (datum- and non-datum-stations).

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
            raise AssertionError('Th DataFrame "obs_df" is not valid.')

    def reduce_to_reference_height(self, target_ref_height: str, verbose: bool = False) -> bool:
        """Reduce the observed gravity to the specified target reference height.

        Notes
        -----
        - For this reduction vertical gravity gradients are required!

        Parameters
        ----------
        target_ref_height : string, specifying the target reference height type.
            The target reference height type has to be listed in :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        bool : True, if no errors occurred
            True implies that the reduction was successful, or that it was not necessary (already reduced).
        """
        # Initial checks:
        if target_ref_height not in REFERENCE_HEIGHT_TYPE:
            raise ValueError(f'"{target_ref_height}" is an unknown reference height type.')

        # Whether all or none of the data fields have to be None:
        if (self.obs_df['g_red_mugal'].isna().all() != self.obs_df['g_red_mugal'].isna().any()) or \
            (self.obs_df['sd_g_red_mugal'].isna().all() != self.obs_df['sd_g_red_mugal'].isna().any()):
            raise AssertionError('In "obs_df" the columns "g_red_mugal" and/or "sd_g_red_mugal" are not initialized '
                                 'consistently!')

        # Check, if VG are available in self.obs_df:
        if self.obs_df.vg_mugalm.isna().any():
            if verbose:
                print('For the following stations the vertical gravity gradient is not available:')
                print(f'{", ".join(self.obs_df[self.obs_df.vg_mugalm.isna()].station_name.unique())}')
                print('Without vertical gradient the observed gravity cannot be reduced to another height level!')
            return False

        # Initilalize fields for reduced data, if empty:
        if self.obs_df['g_red_mugal'].isna().all():
            self.obs_df['g_red_mugal'] = self.obs_df['g_obs_mugal']
            self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_obs_mugal']
            self.red_reference_height_type = self.obs_reference_height_type
            self.red_tide_correction_type = self.obs_tide_correction_type

        # Reduction not necessary:
        if self.red_reference_height_type == target_ref_height:
            if verbose:
                print(f'Reduced gravity values are already referenced to: {target_ref_height}')
                return True

        # Reduction:
        # # Distance between instrument top and sensor level:
        dst_m = GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[GRAVIMETER_ID_BEV[self.gravimeter_id]]
        if self.red_reference_height_type == 'sensor_height':
            if target_ref_height == 'control_point':
                # + dst_m + dhf_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] + (dst_m + self.obs_df['dhf_m']) * \
                                             self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
            elif target_ref_height == 'instrument_top':
                # + dst_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] + dst_m * self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
            elif target_ref_height == 'ground':
                # + dst_m + dhb_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] + (dst_m + self.obs_df['dhb_m']) * \
                                             self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
        if self.red_reference_height_type == 'instrument_top':
            if target_ref_height == 'sensor_height':
                # - dst_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] - dst_m * self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
            elif target_ref_height == 'ground':
                # + dhb_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] + self.obs_df['dhb_m'] * \
                                             self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
            elif target_ref_height == 'control_point':
                # + dhf_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] + self.obs_df['dhf_m'] * \
                                             self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
        if self.red_reference_height_type == 'ground':
            if target_ref_height == 'instrument_top':
                # - dhb_m - dst_m + dst_m = - dhb_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] - self.obs_df['dhb_m'] * \
                                             self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
            elif target_ref_height == 'sensor_height':
                # - dhb_m - dst_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] - (self.obs_df['dhb_m'] + dst_m) * \
                                             self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
            elif target_ref_height == 'control_point':
                # - dhb_m - dst_m + dhf_m + dst_m = dhf_m - dhb_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] + \
                                             (self.obs_df['dhf_m'] - self.obs_df['dhb_m']) * self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
        if self.red_reference_height_type == 'control_point':
            if target_ref_height == 'instrument_top':
                # - dhf_m - dst_m + dst_m = - dhf_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] - self.obs_df['dhf_m'] * \
                                             self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD
            elif target_ref_height == 'ground':
                # - dhf_m - dst_m + dst_m + dhb_m = dhb_m - dhf_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] + \
                                             (self.obs_df['dhb_m'] - self.obs_df['dhf_m']) * self.obs_df['vg_mugalm']
            elif target_ref_height == 'sensor_height':
                # - dhf_m - dst_m
                self.obs_df['g_red_mugal'] = self.obs_df['g_red_mugal'] - (self.obs_df['dhf_m'] + dst_m) * \
                                             self.obs_df['vg_mugalm']
                self.obs_df['sd_g_red_mugal'] = self.obs_df['sd_g_red_mugal']  # Keep SD

        self.red_reference_height_type = target_ref_height
        return True

    def __str__(self):
        if self.obs_df is not None:
            return f'Survey "{self.name}" with {len(self.obs_df)} observations.'
        else:
            return f'Survey "{self.name}"'


class Campaign:
    """Gravity Campaign dataset.

    A gravity campaign datasets consists of:

    - Campaign name
    - One or more gravity surveys that belong together and
      - Each survey was observed with one gravimeter on a single day
    - Station data (datum an non-datum stations)
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
    """

    def __init__(self,
                 campaign_name,
                 output_directory,
                 surveys=None,  # Always use non-mutable default arguments!
                 stations=None,  # Always use non-mutable default arguments!
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
            Data of known stations (datum- and non-datum-stations). Default=None implies that the campaign will be
            initialized without station data.

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

    def add_survey(self, survey_add: Survey, verbose=False) -> bool:
        """Add a survey to campaign ans specify whether to use it for ths analysis.

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
                print('Warnung: the current campaign already contains a survey named {survey_add.name}.')
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

    def reduce_to_reference_height(self, target_ref_height, verbose=False):
        """Reduce the observed gravity to the specified target reference heights.

        Notes
        -----
        - For this reduction vertical gravity gradients are required. They are obtained from the `Station` object.
          Hence, a `Station` object has to be attached to the Campaign object beforehand.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        target_ref_height : string, specifying the target reference height type.
            The target reference height type has to be listed in :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.
        """
        if verbose:
            print(f'Get VG from stations and reduce to reference height "{target_ref_height}":')
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f' - Survey: {survey_name}')
            survey.obs_df_populate_vg_from_stations(self.stations, verbose=verbose)
            if not survey.reduce_to_reference_height(target_ref_height, verbose=verbose):
                raise AssertionError('Reduction to reference height failed!')

    def __str__(self):
        return f'Campaign "{self.campaign_name}" with {self.number_of_surveys} surveys ' \
               f'and {self.stations.get_number_of_stations} stations.'


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

    surv2 = Survey.from_cg5_obs_file(PATH_OBS_FILE_CG5 + NAME_OBS_FILE_CG5)
    print(surv2)

    camp = Campaign('AD1', 'dir_name', stations=stat, surveys={surv2.name: surv2})
    camp.activate_survey(True)

    surv_bev = Survey.from_bev_obs_file(PATH_OBS_FILE_BEV + NAME_OBS_FILE_BEV, verbose=VERBOSE)

    validity = surv_bev.is_valid_obs_df(verbose=True)

    camp.add_survey(survey_add=surv_bev, verbose=VERBOSE)

    surveys_info_dict = camp.get_survey_names_and_status(verbose=VERBOSE)

    surv.obs_df_populate_vg_from_stations(stat, verbose=True)

    target_ref_h = 'control_point'
    surv.reduce_to_reference_height(target_ref_h, verbose=True)

    camp.reduce_to_reference_height(target_ref_height=target_ref_h, verbose=True)

    pass
