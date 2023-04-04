"""Classes for modelling relative gravity surveys.

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
import pandas as pd
import numpy as np
import datetime as dt
import os

import gravtools.models.lsm_diff
import gravtools.tides.longman1959
from gravtools.settings import SURVEY_DATA_SOURCE_TYPES, TIDE_CORRECTION_TYPES, DEFAULT_GRAVIMETER_TYPE_CG5_SURVEY, \
    REFERENCE_HEIGHT_TYPE, BEV_GRAVIMETER_TIDE_CORR_LOOKUP, GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m, \
    GRAVIMETER_SERIAL_NUMBERS, GRAVIMETER_TYPES, GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE
from gravtools.const import VG_DEFAULT
from gravtools.CG5_utils.cg5_survey import CG5Survey
from gravtools.models.misc import format_seconds_to_hhmmss


class Survey:
    """Gravity survey (instrument-independent).

    A gravity survey object contains all data that belongs to a single field observation job, carried out under the
    same circumstances (same instrument, measurement area, drift- and datum-point planing, same observer, etc.). A
    survey is usually observed on the same day, under similar conditions by the same operator.

    All analytical observation-level reductions and corrections are applied here. The reduced g values are stored in the
    columns `g_red_mugal` and the according standard deviation in `sd_g_red_mugal` in `obs_df`. The following
    corrections/reductions are supported:

    - Tidal corrections of observations

      - Corrections of the CG-5 instruments provided in the observation files (Longman, 1959)
      - No corrections
      - Longman (1959) model, evaluated in GravTools for the longitude, latitude , altitude and UTC time stamp of each observation. The correction is calculated for the middle of the readint time (obs_epoch + duration_sec/2)

    - Reduction of the observed gravity to different reference height levels

      - Control point level
      - Ground level
      - Level of the gravimeter top
      - Sensor level

    Notes
    -----
    Basically it is possible to initialize an empty survey, by just defining the survey's name on the instantiation
    step. In this case all other (class and instance) attributes are initialized with the default values listed below
    in the Attributes section.

    The observation reference time in GravTools is equal to the start time of an instrument reading - in accordance with
    the CG-5 observation files. This has to be taken into account when calcualting time dependent corrections, such as
    tidal corrections.

    Attributes
    ----------
    name : str
        Name of the survey. This parameter is mandatory!
    date :
        Date of te survey.
    operator : str, optional (default='')
        Name of the responsible operator that carried out the observations of this survey.
    gravimeter_serial_number : str, optional (default='')
        Valid gravimeter types have to be listed in :py:obj:`gravtools.settings.GRAVIMETER_SERIAL_NUMBERS`.
    gravimeter_type: str, optional (default='')
        Valid gravimeter types have to be listed in :py:obj:`gravtools.settings.GRAVIMETER_TYPES`.
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
    setup_tide_correction_type : str, optional (default='')
        Type of the tidal corrections applied on the observations that are used to calculate the setup data in the
        `setup_df` dataframe. Valid entries have to be listed in :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`.
        `Empty string`, if corrected setup data has not been calculated so far (`setup_df` = None).
    setup_reference_height_type : str, optional (default='')
        Reference height reduction type applied on the observations that are used to calculate the setup data in the
        `setup_df` dataframe. Valid entries have to be listed in :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`.
        `Empty string`, if corrected setup data has not been calculated so far (`setup_df` = None).
    keep_survey : bool (default=True)
        Flag that indicates whether this survey will be used to derive setup observations.`True` is the
        default and implies that this survey is considered. This flag is independent of the `keep_obs` flags
        in the `obs_df` dataframe (used to flag individual observations).
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
            that was observed at least twice for drift control (preferably at the beginning and end of the line). If
            `None`, loops were not defined. The purpose of splitting surveys into loops is to carry out drift correction
            for individual loops (shorter time period) rather than fo the complete survey.
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
            Reference epoch of the observation. Per default the start epoch (!) of an observation. Be aware that for the
            determination of tidal corrections by the CG-5 built-in model (Longman, 1959) the middle of the observation
            with the duration dur_sec [sec] is used (obs_epoch + duration_sec/2)!
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
        - duration_sec : int
            Duration of each observation from the CG-5 observation files. Given in seconds.

    ref_delta_t_dt : datetime, optional (default=None)
        Reference time for relative times (), e.g. reference time t0 the for drift polynomial adjustment.
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
        - sd_g_mugal : float
            Standard deviation of `g_mugal` [µGal]
        - epoch_unix : float
            Reference epoch of `g_mugal` in unix time [sec]
        - epoch_dt : :py:obj:`datetime.datetime`; timezone aware, if possible
            Reference epoch of `g_mugal`.
        - delta_t_h : float
            Time span since reference time :py:obj:`Survey.ref_delta_t_dt` in hours.
        - sd_setup_mugal : float
             Standard deviation of active observations in this setup [µGal].
        - number_obs : int
            Number of observations in a setup.
    setup_obs_list_df : :py:obj:`pandas.core.frame.DataFrame`, optional (default=None)
        List of observations in the `obs_df` dataframe that were used to calculate the setup observations in the
        `setup_df` dataframe.

        - station_name : str
            Name of the station.
        - obs_epoch : :py:obj:`datetime.datetime`; timezone aware, if possible
            Reference epoch of the observation
        - keep_obs : bool (default=True)


    """

    _OBS_DF_COLUMNS = (
        'station_name',  # Name of station (str)
        'setup_id',  # unique ID of setup (int)
        'loop_id',  # Line ID, optional (default=None) (int)
        'lon_deg',  # Longitude [deg], optional (float)
        'lat_deg',  # Latitude [deg], optional (float)
        'alt_m',  # Altitude [m], optional (float)
        'obs_epoch',  # Observation epoch (datetime object, TZ=<UTC>), start of instrument reading!
        'g_obs_mugal',  # observed g from obs file [µGal] (float)
        'sd_g_obs_mugal',  # Standard deviation of g observed from obs file (float) [µGal]
        'g_red_mugal',  # Reduced gravity observation at station (float) [µGal]
        'sd_g_red_mugal',  # Standard deviation of the reduced gravity (float) [µGal]
        'corr_terrain',  # Terrain correction [??]
        'corr_tide_mugal',  # Tidal correction loaded from input file [µGal], optional (e.g. from CG5 built-in model)
        'temp',  # Temperature [mK], optional
        'tiltx',  # [arcsec], optional
        'tilty',  # [arcsec], optional
        'dhf_m',  # Distance between instrument top and physical reference point (float) [m]
        'dhb_m',  # Distance between instrument top and ground (float) [m]
        'keep_obs',  # Remove observation, if false (bool)
        'vg_mugalm',  # vertical gradient [µGal/m]
        'corr_tide_red_mugal',  # Alternative tidal correction [µGal], optional
        'duration_sec',  # Duration [sec]
    )

    _SETUP_DF_COLUMNS = (
        'station_name',  # Name of station (str)
        'setup_id',  # Unique ID of setup (int)
        'g_mugal',  # Variance weighted mean of all active observations in setup [µGal] (float)
        'sd_g_mugal',  # Standard deviation of `g_mugal` (float) [µGal]
        'epoch_unix',  # Reference epoch of `g_mugal` (unix time [sec])
        'epoch_dt',  # Reference epoch of `g_mugal` (datetime obj)
        'delta_t_h',  # Time span since reference time [hours]
        'delta_t_campaign_h',  # Time span since reference time [hours]
        'sd_setup_mugal',  # Standard deviation of active observations in this setup [µGal]
        'number_obs',  # Number of observations in a setup
        'dhf_sensor_m',  # Vertical distance between control point and sensor height
    )

    _SETUP_OBS_LIST_DF_COLUMNS = (
        'station_name',  # Name of station (str)
        'obs_epoch',  # Observation epoch (datetime object, TZ=<UTC>), start of instrument reading!
        'keep_obs',  # If False, the observation is not used for calculating setup observations (bool)
    )

    def __init__(self,
                 name,
                 date=None,
                 operator='',
                 institution='',
                 gravimeter_serial_number='',
                 gravimeter_type='',
                 data_file_name='',
                 data_file_type='',
                 obs_df=None,
                 obs_tide_correction_type='',  # of "g_obs_mugal"
                 obs_reference_height_type='',  # of "g_obs_mugal"
                 red_tide_correction_type='',  # of "g_red_mugal"
                 red_reference_height_type='',  # of "g_red_mugal"
                 setup_tide_correction_type='',
                 setup_reference_height_type='',
                 keep_survey=True,  # Flag
                 setup_df=None,
                 ref_delta_t_dt=None,  # Datetime object (UTC)
                 setup_obs_list_df=None  #
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

        # gravimeter_serial_number:
        if isinstance(gravimeter_serial_number, str):
            if gravimeter_serial_number:
                if gravimeter_serial_number in GRAVIMETER_SERIAL_NUMBERS.keys():
                    self.gravimeter_serial_number = gravimeter_serial_number
                else:
                    raise ValueError('"gravimeter_serial_number" needs to be a key in GRAVIMETER_SERIAL_NUMBERS')
            else:
                self.gravimeter_serial_number = gravimeter_serial_number  # ''
        else:
            raise TypeError('"gravimeter_serial_number" needs to be a string')

        # gravimeter_type:
        if isinstance(gravimeter_type, str):
            if gravimeter_type:
                if gravimeter_type in GRAVIMETER_TYPES.keys():
                    self.gravimeter_type = gravimeter_type
                else:
                    raise ValueError('"gravimeter_type" needs to be a key in GRAVIMETER_TYPES')
            else:
                self.gravimeter_type = gravimeter_type  # ''
        else:
            raise TypeError('"gravimeter_type" needs to be a string')

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
                    raise ValueError('"data_file_type" needs to be a key in SURVEY_DATA_SOURCE_TYPES')
            else:
                if self.data_file_name:  # Not empty
                    raise ValueError('If a data file is specified ("data_file_name" not empty), "data_file_type" has '
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

        # setup_tide_correction_type:
        if isinstance(setup_tide_correction_type, str):
            if setup_tide_correction_type:
                if setup_tide_correction_type in REFERENCE_HEIGHT_TYPE.keys():
                    self.setup_tide_correction_type = setup_tide_correction_type
                else:
                    raise ValueError('"setup_tide_correction_type" needs to be a key in REFERENCE_HEIGHT_TYPE '
                                     '({})'.format(', '.join(REFERENCE_HEIGHT_TYPE.keys())))
            else:
                self.setup_tide_correction_type = setup_tide_correction_type  # ''
        else:
            raise TypeError('"setup_tide_correction_type" needs to be a string')

        # setup_reference_height_type:
        if isinstance(setup_reference_height_type, str):
            if setup_reference_height_type:
                if setup_reference_height_type in REFERENCE_HEIGHT_TYPE.keys():
                    self.setup_reference_height_type = setup_reference_height_type
                else:
                    raise ValueError('"setup_reference_height_type" needs to be a key in REFERENCE_HEIGHT_TYPE '
                                     '({})'.format(', '.join(REFERENCE_HEIGHT_TYPE.keys())))
            else:
                self.setup_reference_height_type = setup_reference_height_type  # ''
        else:
            raise TypeError('"setup_reference_height_type" needs to be a string')

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

        if ref_delta_t_dt is not None:
            if not isinstance(ref_delta_t_dt, dt.datetime):
                raise TypeError('`ref_delta_t_dt` needs to be a datetime object.')
        self.ref_delta_t_dt = ref_delta_t_dt

        # setup_obs_list_df:
        if setup_obs_list_df is not None:
            if isinstance(setup_obs_list_df, pd.DataFrame):
                # Check if setup_df contains exactly all columns defined by self._SETUP_DF_COLUMNS:
                if all([item for item in obs_df.columns.isin(self._SETUP_OBS_LIST_DF_COLUMNS)]) and \
                        obs_df.shape[1] == len(self._SETUP_OBS_LIST_DF_COLUMNS):
                    self.setup_obs_list_df = setup_obs_list_df
                else:
                    raise ValueError(
                        '"setup_obs_list_df" needs the following columns:{}'.format(', '.join(self._SETUP_OBS_LIST_DF_COLUMNS)))
            else:
                raise TypeError('"setup_obs_list_df" needs to be a pandas DataFrame.')
        else:
            self.setup_obs_list_df = setup_obs_list_df  # None

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
                if obs_df['obs_epoch'].dt.tz.zone != 'UTC':  # Change TZ to <UTC>
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
                   gravimeter_type=DEFAULT_GRAVIMETER_TYPE_CG5_SURVEY,
                   gravimeter_serial_number=cg5_survey.survey_parameters.instrument_sn,
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
            if df['obs_epoch'].dt.tz.zone != 'UTC':  # Change TZ to <UTC>
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
            if verbose:
                print(f'Warning: For gravimeter ID "{gravimeter_id}" the tide correction type is unknown '
                      f'(not specified in settings.BEV_GRAVIMETER_TIDE_CORR_LOOKUP). '
                      f'It is set to "{obs_tide_correction_type}".')

        gravimeter_serial_number = ''
        # Get gravimeter type and serial number from ID in obs file:
        for serial_number, id_tmp in GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
            if id_tmp == gravimeter_id:
                gravimeter_serial_number = serial_number
                break
        if not gravimeter_serial_number:
            raise AssertionError(
                f'No serial number found that matches the gravimeter id "{gravimeter_id}"in GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE.')
        try:
            gravimeter_type = GRAVIMETER_SERIAL_NUMBERS[gravimeter_serial_number]
        except KeyError:
            raise AssertionError(
                f'serial number "{gravimeter_serial_number}" not found in the lookuptable "GRAVIMETER_SERIAL_NUMBERS".')

        return cls(name=os.path.split(filename)[1],
                   date=survey_date,
                   operator='',
                   institution=institution,
                   gravimeter_type=gravimeter_type,
                   gravimeter_serial_number=gravimeter_serial_number,
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

    def set_reference_time(self, ref_delta_t_dt):
        """Set refernce time for the determination of relative time spans, e.g. for the drift polynomial.

        Parameters
        ----------
        ref_delta_t_dt : datetime object
            Reference time epoch w.r.t. UTC.
        """
        if isinstance(ref_delta_t_dt, dt.datetime):
            self.ref_delta_t_dt = ref_delta_t_dt
        else:
            raise ValueError('`ref_delta_t_dt` needs to be a datetime object.')

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
                dst_m = GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[self.gravimeter_type]
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
                    elif target_tide_corr == 'longman1959':
                        # Calculate corrections:
                        tmp_df = self.obs_df[['obs_epoch', 'duration_sec', 'lon_deg', 'lat_deg', 'alt_m']].copy(
                            deep=True)
                        # Shift obs reference epoch from start to the middle of the gravity reading
                        # (= evaluation time for the tide model):
                        tmp_df['obs_epoch'] = tmp_df['obs_epoch'] + pd.to_timedelta(tmp_df['duration_sec'], 'sec') / 2
                        # Remove TZ info for longman evaluation:
                        tmp_df['obs_epoch'] = tmp_df['obs_epoch'].dt.tz_localize(None)
                        tmp_df['tmp_index'] = tmp_df.index
                        tmp_df = tmp_df.set_index('obs_epoch')
                        tmp_df = gravtools.tides.longman1959.solve_tide_df(tmp_df, lat='lat_deg', lon='lon_deg',
                                                                           alt='alt_m')
                        tmp_df['longman_tide_corr_mugal'] = tmp_df['g0'] * 1e3  # Convert from mGal to µGal
                        tmp_df = tmp_df.set_index('tmp_index')  # Restore numerical index
                        # Apply correction:
                        g_red_mugal = g_red_mugal - tmp_df['longman_tide_corr_mugal']
                        # sd_g_red_mugal = sd_g_red_mugal # Keep SD
                        corr_tide_red_mugal = tmp_df['longman_tide_corr_mugal']
                    # elif target_tide_corr == '......':
                elif self.obs_tide_correction_type == 'cg5_longman1959':
                    if target_tide_corr == 'no_tide_corr':  # Subtract instrumental corrections
                        g_red_mugal = g_red_mugal - obs_df['corr_tide_mugal']
                        # sd_g_red_mugal = sd_g_red_mugal # Keep SD
                        corr_tide_red_mugal = obs_df['corr_tide_mugal']
                        corr_tide_red_mugal.values[:] = 0
                    elif target_tide_corr == 'longman1959':
                        # Calculate corrections:
                        tmp_df = self.obs_df[['obs_epoch', 'duration_sec', 'lon_deg', 'lat_deg', 'alt_m']].copy(
                            deep=True)
                        # Shift obs reference epoch from start to the middle of the gravity reading
                        # (= evaluation time for the tide model):
                        tmp_df['obs_epoch'] = tmp_df['obs_epoch'] + pd.to_timedelta(tmp_df['duration_sec'], 'sec') / 2
                        # Remove TZ info for longman evaluation:
                        tmp_df['obs_epoch'] = tmp_df['obs_epoch'].dt.tz_localize(None)
                        tmp_df['tmp_index'] = tmp_df.index
                        tmp_df = tmp_df.set_index('obs_epoch')
                        tmp_df = gravtools.tides.longman1959.solve_tide_df(tmp_df, lat='lat_deg', lon='lon_deg', alt='alt_m')
                        tmp_df['longman_tide_corr_mugal'] = tmp_df['g0']*1e3  # Convert from mGal to µGal
                        tmp_df = tmp_df.set_index('tmp_index')  # Restore numerical index
                        # Apply correction:
                        g_red_mugal = g_red_mugal - tmp_df['longman_tide_corr_mugal']
                        # sd_g_red_mugal = sd_g_red_mugal # Keep SD
                        corr_tide_red_mugal = tmp_df['longman_tide_corr_mugal']
                    # elif target_tide_corr == '......':
                elif self.obs_tide_correction_type == 'longman1959':
                    if target_tide_corr == 'no_tide_corr':  # Subtract instrumental corrections
                        g_red_mugal = g_red_mugal - obs_df['corr_tide_mugal']
                        # sd_g_red_mugal = sd_g_red_mugal # Keep SD
                        corr_tide_red_mugal = obs_df['corr_tide_mugal']
                        corr_tide_red_mugal.values[:] = 0
                    elif target_tide_corr == 'cg5_longman1959':  # Add instrumental corrections
                        g_red_mugal = g_red_mugal + obs_df['corr_tide_mugal']
                        # sd_g_red_mugal = sd_g_red_mugal # Keep SD
                        corr_tide_red_mugal = obs_df['corr_tide_mugal']
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

    def autselect_tilt(self, threshold_arcsec: int, setup_id: int = None, verbose: bool = False):
        """Deactivate all observations of the survey or of a setup with a tilt larger than the defined threshold.

        Parameters
        ----------
        threshold_arcsec : int
            Observations in this survey or the specified setup are deactivated, if their tilt in X or Y direction
            (columns `tiltx` and `tilty` in :py:obj:`.Survey.obs_df`) exceeds the given threshold [arcsec].
        setup_id : int (default=None)
            `None` implies that this autoselection function is applied on all observations of this survey. Otherwise,
            the autoselection function is only applied on observations of the setup wirth the provided ID (`setup_id`)
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        filter_tilt = (abs(self.obs_df['tiltx']) > threshold_arcsec) | (abs(self.obs_df['tilty']) > threshold_arcsec)
        if setup_id is not None:  # Apply on whole survey
            filter_tilt = filter_tilt & (self.obs_df['setup_id'] == setup_id)
        self.obs_df.loc[filter_tilt, 'keep_obs'] = False
        if verbose:
            print(f'Removed observations due to tilt threshold ({threshold_arcsec} asec): {filter_tilt[filter_tilt].count()}')

    def autselect_duration(self, threshold_sec: int, setup_id: int = None, verbose: bool = False):
        """Detect and deactivate all observations with a measurement duration smaller than the threshold.

        Parameters
        ----------
        threshold_sec : int
            Threshold for the measurement duration.
        setup_id : int (default=None)
            `None` implies that this autoselection function is applied on all observations of this survey. Otherwise,
            the autoselection function is only applied on observations of the setup wirth the provided ID (`setup_id`)
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        filter_duration = self.obs_df['duration_sec'] < threshold_sec
        if setup_id is not None:  # Apply on whole survey
            filter_duration = filter_duration & (self.obs_df['setup_id'] == setup_id)
        self.obs_df.loc[filter_duration, 'keep_obs'] = False
        if verbose:
            print(f'Removed observations due to duration threshold ({threshold_sec} sec): {filter_duration[filter_duration].count()}')

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
        for setup_id in setup_ids:
            filter_id = self.obs_df['setup_id'] == setup_id
            if verbose:
                print(f' - setup ID: {setup_id}')
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
            self.setup_reference_height_type = ''
            self.setup_tide_correction_type = ''
            self.setup_obs_list_df = None

    def calculate_setup_data(self, obs_type='reduced',
                             ref_delta_t_campaign_dt=None,
                             active_obs_only_for_ref_epoch=True,
                             verbose=False):
        """Accumulate all active observation within each setup and calculate a single representative pseudo observation.

        Parameters
        ----------
        obs_type : str, 'observed' or 'reduced' (default)
            Defines whether the observed (as loaded from an observation file) or the reduced observations from
            `self.obs_df` are used to determine the weighted mean values per setup.
        ref_delta_t_campaign_dt : datetime object, optional (default = None)
            Reference time for calculation of `delta_t_campaign_h`. `None` implies that the reference epoch is not
            available/defined. In the latter case `delta_t_campaign_h` is `None` for all setup observations.
        active_obs_only_for_ref_epoch: bool, optional (default=True)
            `True` implies that the relative reference epochs are determined by considering active observations only.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """

        if verbose:
            print(f'Calculate setup data for survey {self.name}')

        _VALID_OBS_TYPES = ('observed', 'reduced',)
        flag_calculate_delta_t_campaign_h = False

        self.reset_setup_data(verbose)

        # Initial checks:
        if obs_type not in _VALID_OBS_TYPES:
            raise AssertionError(f'Invalid observation type: {obs_type} (valid: "reduced" or "observed").')
        if self.obs_df is None:
            raise AssertionError('Observation dataframe is empty!')

        # Get all active observations:
        tmp_filter = self.obs_df['keep_obs']
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

            # Check input
            if ref_delta_t_campaign_dt is None:  # No reference time for the campaign defined
                flag_calculate_delta_t_campaign_h = False
            elif isinstance(ref_delta_t_campaign_dt, dt.datetime):  # Input is OK!
                flag_calculate_delta_t_campaign_h = True
            else:
                raise TypeError('`ref_delta_t_campaign_dt` needs to be a datetime object!')
            # - SD == 0? => This would create a divided by zero error!
            if obs_type == 'reduced':
                tmp_filter = active_obs_df['sd_g_red_mugal'] <= 0.0
            elif obs_type == 'observed':
                tmp_filter = active_obs_df['sd_g_red_mugal'] <= 0.0
            if tmp_filter.any():
                error_str = active_obs_df.loc[tmp_filter, ['station_name', 'obs_epoch']].to_string()
                raise AssertionError(f'The SD of the following observations ({obs_type}) in the survey {self.name}'
                                     f' is <= 0.0 µGal (invalid!):\n{error_str}')

            # Determine reference time for the survey:
            if active_obs_only_for_ref_epoch:
                ref_delta_t_dt = active_obs_df['obs_epoch'].min()  # First observation epoch in survey (active only)
            else:
                ref_delta_t_dt = self.obs_df['obs_epoch'].min()  # First observation epoch (also inactive obs)

            # Initialize columns lists for creating dataframe:
            station_name_list = []
            setup_id_list = []
            g_mugal_list = []
            sd_g_red_mugal_list = []
            obs_epoch_list_unix = []
            obs_epoch_list_dt = []
            delta_t_h_list = []
            delta_t_campaign_h_list = []
            sd_setup_mugal_list = []
            number_obs_list = []
            dhf_sensor_m_list = []

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

                # Standard deviation of active observations within setup:
                # - If less than 2 active observations in setup => Calculation not possible => NaN
                if len(g_mugal) >= 2:
                    sd_setup_mugal_list.append(g_mugal.std(ddof=1))  # degree of freedom = (len(g_mugal) - 1)
                else:
                    sd_setup_mugal_list.append(np.nan)

                number_obs_list.append(len(g_mugal))  # Number of observations

                # Get vertical distance between sensor height and control point:
                dhf_sensor_m_list.append(active_obs_df.loc[tmp_filter, 'dhf_m'].values[0] + GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[self.gravimeter_type])

                # observation epoch (UNIX timestamps in full seconds):
                obs_epochs_series = active_obs_df.loc[tmp_filter, 'obs_epoch']
                unix_obs_epochs = obs_epochs_series.values.astype(np.int64) / 10 ** 9
                unix_setup_epoch = np.sum(unix_obs_epochs * weights) / np.sum(weights)

                # Reference epoch as datetime object (TZ=<UTC>):
                dt_setup_epoch = dt.datetime.utcfromtimestamp(unix_setup_epoch)
                dt_setup_epoch = dt_setup_epoch.replace(tzinfo=dt.timezone.utc)  # TZ = <UTC>

                station_name_list.append(active_obs_df.loc[tmp_filter, 'station_name'].unique()[0])
                setup_id_list.append(setup_id)
                g_mugal_list.append(g_setup_mugal)
                sd_g_red_mugal_list.append(sd_g_setup_mugal)
                obs_epoch_list_unix.append(unix_setup_epoch)
                obs_epoch_list_dt.append(dt_setup_epoch)
                delta_t_h_list.append((unix_setup_epoch - (ref_delta_t_dt.value / 10 ** 9)) / 3600.0)
                if flag_calculate_delta_t_campaign_h:
                    delta_t_campaign_h_list.append((unix_setup_epoch - (ref_delta_t_campaign_dt.value / 10 ** 9)) / 3600.0)
                else:
                    delta_t_campaign_h_list.append(None)

            if obs_type == 'observed':
                self.setup_tide_correction_type = self.obs_tide_correction_type
                self.setup_reference_height_type = self.obs_reference_height_type
            elif obs_type == 'reduced':
                self.setup_tide_correction_type = self.red_tide_correction_type
                self.setup_reference_height_type = self.red_reference_height_type
            self.create_setup_obs_list()

            # convert to pd dataframe:
            self.setup_df = pd.DataFrame(list(zip(station_name_list,
                                                  setup_id_list,
                                                  g_mugal_list,
                                                  sd_g_red_mugal_list,
                                                  obs_epoch_list_unix,
                                                  obs_epoch_list_dt,
                                                  delta_t_h_list,
                                                  delta_t_campaign_h_list,
                                                  sd_setup_mugal_list,
                                                  number_obs_list,
                                                  dhf_sensor_m_list)),
                                         columns=self._SETUP_DF_COLUMNS)
            self.set_reference_time(ref_delta_t_dt)  # Save reference time for `delta_t_h`, i.e. for the survey.

    def create_setup_obs_list(self):
        """Create list of all observations that contribute to the calculation of setup data in this Survey object.

        The list is created as pabdas dataframe. It holds information whether an observation in the `obs_df` is
        active or inactive (`keep_obs` flag).
        """
        self.setup_obs_list_df = self.obs_df.loc[:, self._SETUP_OBS_LIST_DF_COLUMNS].copy(deep=True)

    @property
    def number_of_setups(self):
        """Returns the number of setups in the survey"""
        if self.obs_df is not None:
            return len(self.obs_df['setup_id'].unique())
        else:
            return None

    @property
    def number_of_observations(self):
        """Returns the total number of observations in the survey"""
        if self.obs_df is not None:
            return len(self.obs_df)
        else:
            return None

    @property
    def number_of_active_observations(self):
        """Returns the number of active observations in the survey"""
        if self.obs_df is not None:
            return len(self.obs_df[self.obs_df['keep_obs']])
        else:
            return None

    @property
    def start_time(self):
        """Returns datetime of the first observation in the survey."""
        if self.obs_df is not None:
            return self.obs_df['obs_epoch'].min()
        else:
            return None

    @property
    def start_time_hhmmss_str(self):
        """Returns the time of the first observation in hh:mm:ss format as string."""
        if self.obs_df is not None:
            return self.obs_df['obs_epoch'].min().strftime('%H:%M:%S')
        else:
            return ''

    @property
    def end_time(self):
        """Returns datetime of the last observation in the survey."""
        if self.obs_df is not None:
            return self.obs_df['obs_epoch'].max()
        else:
            return None

    @property
    def end_time_hhmmss_str(self):
        """Returns the time of the last observation in hh:mm:ss format as string."""
        if self.obs_df is not None:
            return self.obs_df['obs_epoch'].max().strftime('%H:%M:%S')
        else:
            return ''

    @property
    def duration_hhmmss_str(self):
        """Returns the timespan (Timedalta) between first and last observation in the survey in hh:mm:ss format."""
        if self.obs_df is not None:
            duration = self.obs_df['obs_epoch'].max() - self.obs_df['obs_epoch'].min()
            return format_seconds_to_hhmmss(duration.seconds)
        else:
            return ''

    @property
    def is_active(self):
        """`True` indicates that the survey contains at least one active observation."""
        if self.obs_df is not None:
            return self.obs_df['keep_obs'].any()
        else:
            return False

    @property
    def observed_stations(self):
        """Returns a list of all stations observed in this survey."""
        if self.obs_df is not None:
            return self.obs_df['station_name'].unique().tolist()
        else:
            return None

    @property
    def observed_stations_str(self):
        """Returns a string with a list of all stations observed in this survey."""
        if self.obs_df is not None:
            return f', '.join(self.observed_stations)
        else:
            return None


