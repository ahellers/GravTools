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
import datetime as dt
import os

from gravtools.settings import SURVEY_DATA_SOURCE_TYPES, STATION_DATA_SOURCE_TYPES, GRAVIMETER_ID_BEV, \
    TIDE_CORRECTION_TYPES, DEFAULT_GRAVIMETER_ID_CG5_SURVEY, REFERENCE_HEIGHT_TYPE
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
        Pandas Dataframe that with the columns specified in `_STAT_DF_COLUMNS`. It holds all station information.
    """

    _STAT_DF_COLUMNS = (
        'name',  # Station name, str
        'long_deg',  # longitude [deg], float
        'lat_deg',  # latitude [deg], float
        'height_m',  # Height [m]
        'g_mugal',  # gravity [µGal]
        'g_sig_mugal',  # standard deviation of the gravity [µGal]
        'vg_mugal',  # vertical gradient [µGal/m]
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
            'name',
            'notes',
            'lat_deg',
            'long_deg',
            'height_m',
            'g_mugal',
            'g_sig_mugal',
            'vg_mugal',
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
        idx_series = self.stat_df.name.isin(station_names)
        if verbose:
            deleted_df = self.stat_df[idx_series]
            print(f'Deleted {len(deleted_df)} rows in station dataframe:')
            for index, row in deleted_df.iterrows():
                print(f' - Station {row["name"]:10s} (row index {index:d})')
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
    survey is usually observed on the same day.

    Attributes
    ----------
    # TODO: Add Attribute description!
    """

    _OBS_DF_COLUMNS = (
        'station_name',  # Name of station (str)
        'setup_id',  # unique ID of setup (int)
        'line_id',  # Line ID, optional (default=None) (int)
        'lon_deg',  # Longitude [deg], optional (float)
        'lat_deg',  # Latitude [deg], optional (float)
        'alt_m',  # Altitude [m], optional (float)
        'obs_epoch',  # Observation epoch (datetme obj)
        'g_obs_mugal',  # observed g from obs file [µGal] (float)
        'sd_g_obs_mugal',  # Standard deviation of g observed from obs file (float) [µGal]
        'g_red_mugal',  # Reduced gravity observation at station (float) [µGal]
        'corr_terrain',  # Terrain correction [??]
        'corr_tide',  # Tidal correction [??], optional
        'corr_temp',  # Temperature [??], optional
        'tiltx',  # [??], optional
        'tilty',  # [??], optional
        'dhf_m',  # Distance between instrument top and physical reference point (float) [m]
        'dhb_m',  # Distance between instrument top and ground (float) [m]
        'keep_obs',  # Remove observation, if false (bool)
    )

    # TODO: Get missing infos on columns in CG5 obs file!

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
                 ):

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

    @classmethod
    def from_cg5_survey(cls, cg5_survey):
        """Constructor that generates and populates the survey object from a CG5Survey class object.

        Parameters
        ----------
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
            obs_df = cg5_survey.obs_df.copy(deep=True)  # deep copy => No struggles with references

            # Add missing columns (initialized with default values):
            obs_df['line_id'] = None
            obs_df['g_obs_mugal'] = obs_df['g_mgal'] * 1e3
            obs_df['sd_g_obs_mugal'] = obs_df['sd_mgal'] * 1e3
            obs_df['g_red_mugal'] = None
            obs_df['keep_obs'] = True

            # Rename columns:
            obs_df.rename(columns={'terrain': 'corr_terrain',
                                   'tide': 'corr_tide',
                                   'temp': 'corr_temp'},
                          inplace=True)

            # Drop columns that are not needed any more:
            # - all columns that are not in _OBS_DF_COLUMNS
            columns_to_be_dropped = list(set(obs_df.columns) - set(cls._OBS_DF_COLUMNS))
            obs_df.drop(columns=columns_to_be_dropped, inplace=True)

            # Change column order: https://erikrood.com/Python_References/change_order_dataframe_columns_final.html
            obs_df = obs_df[list(cls._OBS_DF_COLUMNS)]

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
                   )

    @classmethod
    def from_cg5_obs_file(cls, filename):
        """Constructor that generates and populates the survey object directly from a CG5 observation file.

        Parameters
        ----------
        filename : str
            Name (and path) of a CG-5 observation file (text format).

        Returns
        -------
        :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        """
        cg5_survey = CG5Survey(filename)
        return cls.from_cg5_survey(cg5_survey)

    @classmethod
    def from_bev_obs_file(cls, filename):
        """Constructor that generates populates the survey object from an observation file in the legacy BEV format.

        Parameters
        ----------
        filename : str
            Name (and path) of an observation in the legacy BEV format.

        Returns
        -------
        :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        """
        pass
        # TODO

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
    surveys: dict of :py:obj:`.Survey` objects
        Arbitrary number of survey objects.
    stations : :py:obj:`.Station` object
        Data of known stations (datum- and non-datum-stations).
    """

    def __init__(self,
                 campaign_name,
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

        # Check surveys:
        # Format:
        # survey = {<survey_name>: {'survey': <survey_object>,
        #                          'keep_survey': <bool>}
        # }
        if surveys is None:
            surveys = {}
        else:
            if not isinstance(surveys, dict):
                raise TypeError('The argument "survey" needs to be a dict of Survey objects.')
            else:
                for survey_name, survey_data in surveys.items():
                    if not isinstance(survey_name, str):
                        raise TypeError('The argument "survey" needs to be a string.')
                    if not isinstance(survey_data, dict):
                        raise TypeError('The value in the dict "survey" needs to be a dict.')
                    if not isinstance(survey_data['survey'], Survey):
                        raise TypeError('"survey" needs to be a Survey object.')
                    if not isinstance(survey_data['keep_survey'], bool):
                        raise TypeError('"keep_survey" needs to be a boolean type.')
        self.surveys = surveys  # dict: key=Name of Survey, value=Survey object

        # Check stations:
        if stations is None:
            stations = Station()
        else:
            if not isinstance(stations, Station):
                raise TypeError('The argument "stations" needs to be a Station object.')
        self.stations = stations

    def add_survey(self, survey_add: Survey, keep_survey: bool = True, verbose=False):
        """Add a survey to campaign ans specify whether to use it for ths analysis.

        Notes
        -----
        A survey can only be added, f the survey's name is unique within the campaign.

        Parameters
        ----------
        survey_add : :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        keep_survey : bool, optional (default=True)
            Specifies whether the survey is used in the data analysis. 'True' implies that this survey will be used.
        verbose : bool, optional (default=False)
            If True, Warnings are printed to the command line.

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
            self.surveys[survey_add.name] = {'survey': survey_add, 'keep_survey': keep_survey}
            return True

    def remove_survey(self, survey_name: str):
        """Remove survey with the specified name from the campaign

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be removed from the campaign.
        """
        pass
        # TODO

    def activate_survey(self, survey_name: str):
        """Set the survey with the specified name active.

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be set active.
        """
        pass
        # TODO

    def deactivate_survey(self, survey_name: str):
        """Set the survey with the specified name inactive.

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be set inactive.
        """
        pass
        # TODO

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
        for surv_name, surv_data in self.surveys.items():
            info_dict[surv_name] = surv_data['keep_survey']
            if verbose:
                activity_str = lookup_dict[surv_data['keep_survey']]
                print(f' - {surv_name:12s} ({activity_str})')
        return info_dict

    @property
    def number_of_surveys(self) -> int:
        """int : Returns the number of surveys in this campaign."""
        return len(self.surveys)

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

    survey = {surv2.name: {'survey': surv2,
                           'keep_survey': True}
              }
    camp = Campaign('AD1', stations=stat, surveys=survey)

    surv.name += '_2'  # mock 2nd independent survey object
    camp.add_survey(survey_add=surv, keep_survey=True, verbose=VERBOSE)

    surveys_info_dict = camp.get_survey_names_and_status(verbose=VERBOSE)

    pass
