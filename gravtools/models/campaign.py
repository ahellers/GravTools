"""Modelling of relative gravity campaigns.

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

import datetime as dt
import pytz
import numpy as np
import pickle
import os
import sys

import pandas as pd
from gravtools.models.lsm import LSM
from gravtools.models.lsm_diff import LSMDiff
from gravtools.models.vg_lsm import VGLSM
from gravtools.models.lsm_nondiff import LSMNonDiff
from gravtools.models.mlr_bev_legacy import BEVLegacyProcessing
from gravtools.models.survey import Survey
from gravtools.models.station import Station
from gravtools.tides.correction_time_series import CorrectionTimeSeries
from gravtools.settings import ADDITIVE_CONST_ABS_GRTAVITY, GRAVIMETER_TYPES_KZG_LOOKUPTABLE, \
    GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE, GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m, EXPORT_OBS_LIST_COLUMNS, \
    MAX_SD_FOR_EXPORT_TO_NSB_FILE, WRITE_COMMENT_TO_NSB, PICKLE_PROTOCOL_VERSION
from gravtools import __version__ as GRAVTOOLS_VERSION


class Campaign:
    """Gravity Campaign dataset.

    A gravity campaign datasets consists of:

    - Campaign name
    - One or more gravity surveys that belong together and
      - Each survey was observed with one gravimeter on a single day
    - Station data (datum and non-datum stations)
    - Reductions and corrections
      - All observations (from surveys) are corrected and reduced in the same way
    - Time series data for the correction of gravity observations (optional)

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
    ref_delta_t_dt : datetime object
        Reference epoch for relative times within the campaign, e.g. for the determination of the drift polynomials.
        This reference time is equalt to the first (active) observation in the campaign considering all surveys.
    correction_time_series : :py:obj:`CorrectionTimeSeries`
        Contains time series data for correcting gravity observations.

    gravtools_version : str
        Version of the gravtools software that was used to create the dataset.
    """

    def __init__(self,
                 campaign_name,
                 output_directory,
                 surveys=None,  # Always use non-mutable default arguments!
                 stations=None,  # Always use non-mutable default arguments!
                 lsm_runs=None,  # Always use non-mutable default arguments!
                 ref_delta_t_dt=None  # Reference time for drift determination
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

        # Check ref_delta_t_dt:
        if ref_delta_t_dt is not None:
            if not isinstance(ref_delta_t_dt, dt.datetime):
                raise TypeError('`ref_delta_t_dt` needs to be a datetime object.')
        self.ref_delta_t_dt = ref_delta_t_dt

        # Version of gravtools:
        self.gravtools_version = GRAVTOOLS_VERSION

        # Correction time series object:
        self.correction_time_series = CorrectionTimeSeries()

    def add_empty_correction_time_series(self):
        """Adds an empty `CorrectionTimeSeries` object to the campaign.

        Notes
        -----
        This is required, e.g. if a campaign object os loaded into GravTools from a previous GRavTools version without
        support of time series corrections.
        """
        self.correction_time_series = CorrectionTimeSeries()

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
        """
        # Check if a survey with the dame name ("survey_add.name") already exists in this campaign:
        # - Raise warning:
        if survey_add.name in self.surveys.keys():
            raise RuntimeError(f'The campaign already contains a survey named {survey_add.name}. Survey names need to '
                               f'be unique within a campaign!')
        else:
            # Add survey:
            self.surveys[survey_add.name] = survey_add
            if verbose:
                print(f"Survey {survey_add.name} added to the campaign.")

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
        except Exception:
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
    def survey_names(self):
        """Returns a list with the names of all surveys in the campaign."""
        return list(self.surveys.keys())

    @property
    def number_of_surveys(self) -> int:
        """int : Returns the number of surveys in this campaign."""
        return len(self.surveys)

    @property
    def number_of_stations(self) -> int:
        """int : Returns the number of stations in this campaign."""
        return self.stations.get_number_of_stations

    def reduce_observations_in_all_surveys(self,
                                           target_ref_height=None,
                                           target_tide_corr=None,
                                           target_atm_pres_corr=None,
                                           atm_pres_admittance=None,
                                           tide_corr_timeseries_interpol_method='',
                                           verbose=False):
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
        target_atm_pres_corr : str (default = `None`)
            Specifying the atmospheric press ure correction type to be applied to all surveys. Valid types to be listed
            in :py:obj:`gravtools.settings.ATM_PRES_CORRECTION_TYPES`. Default is `None` indicating that the respective
            corrections of the input data are not changed.
        atm_pres_admittance : float, optional (default = `None`)
            Admittance factor for the determination of pressure corrections based on the difference between measured and
            normal air pressure. If `target_atm_pres_corr` is not None (i.e. atmospheric pressure corrections will be
            calculated), `atm_pres_admittance` has to be provided too (float). Otherwise, an error is raised.
        tide_corr_timeseries_interpol_method : str, optional (default='')
            Interpolation method used to calculate tidal corrections from time series data. If tidal corrections are
            obtained from other sources or models, this attribute is irrelevant and has to be empty!
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.
        """
        if verbose:
            print(f'## Reduce all observation in this campaign:')
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f'Survey {survey_name}:')
            if target_ref_height is not None:
                if verbose:
                    print(f' - Get vertical gradients')
                survey.obs_df_populate_vg_from_stations(self.stations, verbose=verbose)
            survey.reduce_observations(
                target_ref_height=target_ref_height,
                target_tide_corr=target_tide_corr,
                target_atm_pres_corr=target_atm_pres_corr,
                atm_pres_admittance=atm_pres_admittance,
                tide_corr_timeseries_interpol_method=tide_corr_timeseries_interpol_method,
                correction_time_series=self.correction_time_series,
                verbose=verbose)

    def add_stations_from_oesgn_table_file(self, oesgn_filename, is_datum=False, verbose=False):
        """Add station from an OESGN table file.

        Parameters
        ----------
        oesgn_filename : string, specifying the path and filename of the OESGN file
            Stations in the specified OESGN table file are added to the campaign.
        is_datum : bool, optional (default = False)
            `True` indicates that all loaded OESGN stations are initially selected as datum stations (is_datum=True)
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        self.stations.add_stations_from_oesgn_table(filename=oesgn_filename, is_datum=is_datum, verbose=verbose)

    def add_stations_from_csv_file(self, csv_filename, verbose=False):
        """Add station from a CSV file.

        Parameters
        ----------
        csv_filename : string, specifying the path and filename of the station csv file
            Stations in this csv file are added to the campaign.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        self.stations.add_stations_from_csv_file(filename=csv_filename, verbose=verbose)

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

    def calculate_setup_data(self,
                             obs_type='reduced',
                             active_obs_only_for_ref_epoch=True,
                             method='variance_weighted_mean',
                             method_sd='sd_from_obs_file',
                             default_sd_mugal=100.0,
                             verbose=False):
        """Calculate accumulated pseudo observations for each active setup in all active surveys.

        Notes
        -----
        Two relative reference epochs are calculated for each setup: (a) w.r.t. the first (active) observation in the
        whole campaign and (b) w.r.t. the first (active) observation in each survey. Both reference time do not differ
        for the first survey on a campaign. Whether active observations only are considered is defined by the input
        parameter `active_obs_only_for_ref_epoch`.

        Parameters
        ----------
        obs_type : str, 'observed' or 'reduced' (default)
            Defines whether the observed (as loaded from an observation file) or the reduced observations from
            `self.obs_df` are used to determine the weighted mean values per setup.
        active_obs_only_for_ref_epoch: bool, optional (default=True)
            `True` implies that the relative reference epochs are determined by considering active observations only.
        method : str, optional (default=`variance_weighted_mean`)
            Select method for the calculation of setup data. `variance_weighted_mean` implies that setup observations
            (observed gravity, standard deviations and reference time) are calculated by variance weighted mean of the
            individual observations. `individual_obs` implies that the original observations are used as setup data
            without any aggregation.
        method_sd : str, optional (default='sd_from_obs_file')
            Method for the determination of standard deviations (SD) of setup observations. `sd_from_obs_file` implies that
            SD are taken from the observation file. `sd_default_per_obs` and `sd_default_per_setup` imply that the
            given default SD is used, where the default SD is applied the individual observations in the first case and
            to setups in the second case. If applied to observations, the number of observations per setup still plays a
            role for weighting the setup observations in the adjustment.
        default_sd_mugal : float, optional (default=100.0)
            Default standard deviation [µGal] that is used to determine the SD of setup observations when `method_sd` is
            `sd_default_per_obs` or `sd_default_per_setup`
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Get reference epoch:
        self.ref_delta_t_dt = self.get_epoch_of_first_observation(active_obs_only_for_ref_epoch)

        # Loop over all surveys in the campaign:
        if verbose:
            print(f'Calculate setup data:')
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f' - Survey: {survey_name}')
            if survey.is_active and survey.keep_survey:
                survey.calculate_setup_data(obs_type=obs_type,
                                            ref_delta_t_campaign_dt=self.ref_delta_t_dt,
                                            active_obs_only_for_ref_epoch=active_obs_only_for_ref_epoch,
                                            method=method,
                                            method_sd=method_sd,
                                            default_sd_mugal=default_sd_mugal,
                                            verbose=verbose)
            else:
                survey.reset_setup_data(verbose)  # Remove setup data from previous calculations

    def get_epoch_of_first_observation(self, active_obs_only_for_ref_epoch=True):
        """Returns the epoch of the first (active) observation in this campaign.

        Parameters
        ----------
        active_obs_only_for_ref_epoch: bool, optional (default=True)
            `True` implies that the reference epoch is determined by considering active observations only.

        Returns
        -------
        datetime object
        """
        first_obs_epoch_dt = None
        flag_first_survey_in_campaign = True

        for survey_name, survey in self.surveys.items():

            # Set filter to select active observations only:
            if active_obs_only_for_ref_epoch:
                filter_tmp = survey.obs_df['keep_obs'] == True  # Select active observations only
            else:
                filter_tmp = [True] * len(survey.obs_df)  # Select all observations

            if len(survey.obs_df.loc[filter_tmp, 'obs_epoch']) > 0:
                if flag_first_survey_in_campaign:
                    flag_first_survey_in_campaign = False
                    first_obs_epoch_dt = survey.obs_df.loc[filter_tmp, 'obs_epoch'].min()
                else:
                    if survey.obs_df.loc[filter_tmp, 'obs_epoch'].min() < first_obs_epoch_dt:
                        first_obs_epoch_dt = survey.obs_df.loc[filter_tmp, 'obs_epoch'].min()
        return first_obs_epoch_dt

    def initialize_and_add_lsm_run(self, lsm_method, comment='', write_log=True):
        """Initialize and add an least-squares adjustment run (object) to the campaign.

        Parameters
        ----------
        lsm_method : str
            Defines the adjustment method. Has to b listed in :py:obj:`gravtools.settings.ADJUSTMENT_METHODS`.
        comment : str, optional (default = '')
            Optional comment on the adjustment run.
        write_log : bool, optional (default=True)
            Flag that indicates whether log string should be written or not.
        """
        # Initialize LSM object:
        if lsm_method == 'LSM_diff':
            lsm_run = LSMDiff.from_campaign(self, comment, write_log)
        elif lsm_method == 'LSM_non_diff':
            lsm_run = LSMNonDiff.from_campaign(self, comment, write_log)
        elif lsm_method == 'MLR_BEV':
            lsm_run = BEVLegacyProcessing.from_campaign(self, comment, write_log)
        elif lsm_method == 'VG_LSM_nondiff':
            lsm_run = VGLSM.from_campaign(self, comment, write_log)
        else:
            raise AssertionError(f'Unknown LSM method: {lsm_method}')
        # Add LSM object to campaign:
        self.lsm_runs.append(lsm_run)

    @property
    def lsm_run_times(self):
        """Returns a list of lsm-run times/dates that can be used to identify individual runs.

        Returns
        -------
        list : List of string stating the epochs of lsm adjustment runs that can be used to identify individual runs.
        """
        lsm_run_times = []
        for lsm_run in self.lsm_runs:
            lsm_run_times.append(lsm_run.time_str)
        return lsm_run_times

    def delete_lsm_run(self, idx):
        """Delete the LSM run with the specified index in the list.

        Parameters
        ----------
        idx : int
            Index of the LSM object in the list.
        """
        if idx != -1:
            del self.lsm_runs[idx]

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

    def write_nsb_file(self, filename: str, lsm_run_index, vertical_offset_mode: str = 'first',
                       exclude_datum_stations=False, formal_error_type='se', verbose=False):
        """Write the results of an LSM run to an nsb file (input for NSDB database).

        Notes
        -----
        A nsb file can only be written, if station results are available which is not the case e.g. for the estimation
        of vertical gravity gradients!

        Parameters
        ----------
        filename : str
            Name and path of the output nsb file (e.g. /home/johnny/example.nsb)
        lsm_run_index : int
            Index of the lsm run in `campaign.lsm_runs` of which the results are exported to the nsb file.
        vertical_offset_mode : str, optional (default='first')
            Defines how the vertical offsets between instrument top and ground (dhb) and reference marker (dhf),
            respectively, are determined in the case of multiple measurements (setups) on the same point. In the nsb
            file only one dhf/dhb pair per station is allowed. Two options: (1) 'first' indicates that dhb and dhf are
            taken from the first setup at a station. (2) 'mean' indicates that mean values over all setups are taken.
        exclude_datum_stations : boolean, optional (default=False)
            `True` indicates that datum stations are excluded from the nsb file.
        formal_error_type : str, optional (default=`se`, alternative: `sd`)
            Select which formal errors (of station gravities) are exported to the nsd file. The user may choose
            between post-fit standard deviations (`sd`) and standard errors ('se'). The standard errors consider the
            noise floor defined in the estimation settings.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Init.:
        nsb_string = ''
        formal_error_type_options = ('sd', 'se')

        if formal_error_type not in formal_error_type_options:
            raise AssertionError(f'Invalid formal error type! Valid: "sd" and "se".')

        # Get and prepare data:
        # - lsm_run
        lsm_run = self.lsm_runs[lsm_run_index]
        results_stat_df = lsm_run.get_results_stat_df

        # Check, if station results are available (e.g. nor the case for VG estimation):
        if results_stat_df is not None:

            # Check if the required columns are available:
            if 'g_est_mugal' in results_stat_df.columns and 'sd_g_est_mugal' in results_stat_df.columns:

                # Check if the data is suitable for export to the nsb file:
                if results_stat_df.loc[results_stat_df['sd_g_est_mugal'] > MAX_SD_FOR_EXPORT_TO_NSB_FILE,
                                       'sd_g_est_mugal'].any():
                    raise AssertionError(f"The SD of at least one station's estimated gravity is larger than {MAX_SD_FOR_EXPORT_TO_NSB_FILE} µGal! ")

                # Loop over stations in results dataframe:
                for index, row in results_stat_df.iterrows():

                    # Skip datum stations:
                    if exclude_datum_stations:
                        if row['is_datum']:
                            continue

                    station_name = row['station_name']
                    observed_in_surveys = []
                    dhb_list_m = []
                    dhf_list_m = []

                    # Get surveys at which the station was observed:
                    for survey_name, setup_data in lsm_run.setups.items():
                        setup_df = setup_data['setup_df']
                        if len(setup_df.loc[setup_df['station_name'] == station_name]) > 0:  # was observed in this setup!
                            observed_in_surveys.append(survey_name)
                            obs_df = self.surveys[survey_name].obs_df
                            setup_ids = obs_df.loc[obs_df['station_name'] == station_name, 'setup_id'].unique()
                            setup_ids = setup_df.loc[setup_df['station_name'] == station_name, 'setup_id'].to_list()
                            # Get list of dhb and dhf:
                            for setup_id in setup_ids:
                                dhb_list_m.append(obs_df.loc[obs_df['setup_id'] == setup_id, 'dhb_m'].values[0])
                                dhf_list_m.append(obs_df.loc[obs_df['setup_id'] == setup_id, 'dhf_m'].values[0])

                    if vertical_offset_mode == 'first':
                        dhb_m = dhb_list_m[0]
                        dhf_m = dhf_list_m[0]
                    elif vertical_offset_mode == 'mean':
                        dhb_m = np.mean(dhb_list_m)
                        dhf_m = np.mean(dhf_list_m)

                    # Get gravimeter S/N and gravimeter type of first survey in the list:

                    if len(observed_in_surveys) > 1:
                        if verbose:
                            print(f'WARNING: station {station_name} was observed in {len(observed_in_surveys)} surveys! Hence, '
                                  f'the gravimeter serial number/type and the observation date may be ambiguous in the nsb file!')
                            print(
                                f' - {station_name} was observed the following surveys: {", ".join(observed_in_surveys)}')
                        # Get the latest survey in which the station was observed:
                        latest_survey_name = ''
                        latest_survey_date = dt.date(1900,1,1)
                        for survey in observed_in_surveys:
                            if self.surveys[survey].date > latest_survey_date:
                                latest_survey_date = self.surveys[survey].date
                                latest_survey_name = survey
                    else:
                        latest_survey_name = observed_in_surveys[0]
                    gravimeter_type = self.surveys[latest_survey_name].gravimeter_type
                    gravimeter_serial_number = self.surveys[latest_survey_name].gravimeter_serial_number
                    date_str = self.surveys[latest_survey_name].date.strftime('%Y%m%d')

                    # Comment string:
                    # - Max. 5 characters!
                    if WRITE_COMMENT_TO_NSB == 'cg5_serial_number':
                        comment_str = str(gravimeter_serial_number)
                    elif WRITE_COMMENT_TO_NSB == 'inst_id':
                        comment_str = GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE[gravimeter_serial_number]
                    elif WRITE_COMMENT_TO_NSB == 'gravtools_version':
                        comment_str = 'GT'+''.join(GRAVTOOLS_VERSION.split('.'))
                    else:
                        raise AssertionError(f'Invalid choice for the nsb file comment: {WRITE_COMMENT_TO_NSB}!')

                    if formal_error_type == 'se':
                        formal_error = row['se_g_est_mugal']
                    elif formal_error_type == 'sd':
                        formal_error = row['sd_g_est_mugal']

                    nsb_string += '{:10s} {:8s}  {:9.0f} {:3.0f} {:1s}{:>5s} {:4.0f} {:4.0f}\n'.format(
                        station_name,
                        date_str,
                        row['g_est_mugal'] + ADDITIVE_CONST_ABS_GRTAVITY,
                        formal_error,
                        GRAVIMETER_TYPES_KZG_LOOKUPTABLE[gravimeter_type],
                        comment_str,
                        (dhb_m + GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[gravimeter_type]) * 100,
                        (dhf_m + GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[gravimeter_type]) * 100,
                    )

                # Write file:
                with open(filename, 'w') as out_file:
                    out_file.write(nsb_string)

            else:  # Required columns are not available
                if verbose:
                    print(f'The nsb file cannot be written as the required station data is not available.')

        else:  # No station results available
            if verbose:
                print(f'The nsb file cannot be written as the required station data is not available.')
            
    def write_log_file(self, filename: str, lsm_run_index, verbose=False):
        """Write log file of a selected LSM run.
        
        Parameters
        ----------
        filename : str
            Name and path of the output nsb file (e.g. /home/johnny/example.nsb)
        lsm_run_index : int
            Index of the lsm run in `campaign.lsm_runs` of which the results are exported to the nsb file.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Init.:
        log_string = ''

        # Get and prepare data:
        lsm_run = self.lsm_runs[lsm_run_index]
        log_string = lsm_run.get_log_string

        # Append additional information to log file string:
        time_now_str = dt.datetime.now(tz=pytz.UTC).strftime('%Y-%m-%d, %H:%M:%S %Z')
        append_str = ''
        append_str += f'------------------------------------------\n'
        append_str += f'LSM run comment: {lsm_run.comment}\n'
        append_str += f'Log file created: {time_now_str}\n'
        append_str += f'Log file written with GravTools version: {GRAVTOOLS_VERSION}\n'
        out_string = log_string + '\n' + append_str
        
        # Write file:
        if verbose:
            print(f'Write log file to {filename}.')
        with open(filename, 'w') as out_file:
            out_file.write(out_string)

    def save_to_pickle(self, filename=None, verbose=True):
        """Save the campaign object to a pickle file at the given path.

        Parameters
        ----------
        filename : str, optional (default=`None`)
            Path and name of the pickle file, e.g. /home/user1/data/camp1.pkl. `None` indicates that the campaign object
            is saved to the default output directory past (`campaign.output_directory`). In this case the file is named
            `<campaign_name>.pkl`.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        Returns
        -------
        str : Name and path of the saved file.
        """
        if filename is None:
            filename = os.path.join(self.output_directory, f'{self.campaign_name}.pkl')
        # Open file:
        if verbose:
            print(f'Export campaign data to {filename}.')
        with open(filename, 'wb') as outfile:
            if PICKLE_PROTOCOL_VERSION == '999':
                pickle.dump(self, outfile, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                pickle.dump(self, outfile, protocol=PICKLE_PROTOCOL_VERSION)
        return filename

    @classmethod
    def from_pkl(cls, filename: str, verbose: bool = True):
        """Loads a campaign object from a pickle file.

        Parameters
        ----------
        filename : str
            Name and path of the pickle file containing the campaign data.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        Returns
        -------
        `Campaign` object
        """
        campaign = pd.read_pickle(filename)
        if verbose:
            print(
                f'Loaded campaign "{campaign.campaign_name}" with {campaign.number_of_stations} station(s) and {campaign.number_of_surveys} survey(s).')
        campaign.check_survey_data(verbose=verbose)
        return campaign

    def check_survey_data(self, verbose: bool = True):
        """Check survey data in campaign, correct missing elements or raise an error.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        Returns
        -------
        `Campaign` object
        """
        for survey_name, survey in self.surveys.items():
            # Check observation dataframe:
            is_valid, error_msg = survey.check_obs_df(verbose=verbose)
            if not is_valid:
                if verbose:
                    print('Error: ' + error_msg)
                raise RuntimeError(error_msg)
            # Check survey object attributes:
            survey.init_missing_attributes(verbose=verbose)

    def set_output_directory(self, output_directory):
        """Change the campaign's output directory.

        Parameters
        ----------
        output_directory : str
            New output directory. The specified directory has to exist on the computer/os!
        """
        # Check output directory:
        if not isinstance(output_directory, str):
            raise TypeError('The argument "output_directory" needs to be a string.')
        else:
            if not os.path.isdir(output_directory):
                raise AssertionError(f'The directory "{output_directory}" does not exist!')
        self.output_directory = output_directory

    def write_obs_list_csv(self, filename_csv: str, export_type: str = 'all_obs', verbose: bool = False):
        """Export a list of all observations in the current campaign.

        Parameters
        ----------
        filename_csv : str
            Name and path of the output CSV file.
        export_type : str, optional (default = 'all_obs')
            Defines which observations are exported. There are three options: (1) all observations ('all_obs'), (2) only
            active observations ('active_only') or (3) only inactive observations ('inactive_only').
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Prepare dataframe with all observations of all surveys
        export_survey_df_list = []
        for survey_name, survey in self.surveys.items():
            tmp_obs_df = survey.obs_df.copy(deep=True)
            tmp_obs_df.insert(0, 'survey_name', survey_name)
            export_survey_df_list.append(tmp_obs_df)
        export_obs_df = pd.concat(export_survey_df_list, ignore_index=True, sort=False)
        export_obs_df.sort_values('obs_epoch', inplace=True)

        # Filter data:
        if export_type != 'all_obs':
            if export_type == 'active_only':
                tmp_filter = export_obs_df['keep_obs']
            elif export_type == 'inactive_only':
                tmp_filter = ~export_obs_df['keep_obs']
            export_obs_df = export_obs_df.loc[tmp_filter, :]

        # Export to CSV file:
        if verbose:
            print(f'Write observation list to: {filename_csv}')
        export_obs_df.to_csv(filename_csv, index=False, columns=EXPORT_OBS_LIST_COLUMNS)

    def write_obs_list_of_lsm_run_csv(self, filename_csv: str, lsm_run_index: int,
                                      export_type: str = 'all_obs', verbose: bool = False):
        """Export a list of all observations that were used for calcualting the setup data of a lsm_run.

        Parameters
        ----------
        filename_csv : str
            Name and path of the output CSV file.
        lsm_run_index : int
            Index of the lsm_run for which the observation list should be exported. If `None`, the current observation
            selection (in `Survey.obs_df`) is exported.
        export_type : str, optional (default = 'all_obs')
            Defines which observations are exported. There are three options: (1) all observations ('all_obs'), (2) only
            active observations ('active_only') or (3) only inactive observations ('inactive_only').
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Get LSM run
        try:
            lsm_run = self.lsm_runs[lsm_run_index]
        except:
            raise AssertionError(f'Invalid LSM run index: {lsm_run_index}')

        # Prepare dataframe with all observations of all surveys
        export_survey_df_list = []
        for survey_name, setup in lsm_run.setups.items():
            setup_obs_list_df = setup['setup_obs_list_df'].copy(deep=True)
            setup_obs_list_df.insert(0, 'survey_name', survey_name)
            export_survey_df_list.append(setup_obs_list_df)
        export_df = pd.concat(export_survey_df_list, ignore_index=True, sort=False)
        export_df.sort_values('obs_epoch', inplace=True)

        # Merge with survey data, if available:
        obs_df_list = []
        for survey_name, survey in self.surveys.items():
            obs_df = self.surveys[survey_name].obs_df.copy(deep=True)
            obs_df_list.append(obs_df)
        obs_df_all = pd.concat(obs_df_list, ignore_index=True, sort=False)
        obs_df_all.sort_values('obs_epoch', inplace=True)

        export_df = pd.merge(export_df, obs_df_all, how='left', left_on=['station_name', 'obs_epoch'],
                                     right_on=['station_name', 'obs_epoch']).rename(
            columns={'keep_obs_x': 'keep_obs', 'station_name_x': 'station_name', 'obs_epoch_x': 'obs_epoch'})

        drop_col_list = list(set(export_df.columns) - set(EXPORT_OBS_LIST_COLUMNS))
        export_df.drop(columns=drop_col_list, inplace=True)

        # Filter data:
        if export_type != 'all_obs':
            if export_type == 'active_only':
                tmp_filter = export_df['keep_obs']
            elif export_type == 'inactive_only':
                tmp_filter = ~export_df['keep_obs']
            export_df = export_df.loc[tmp_filter, :]

        # Export to CSV file:
        if verbose:
            print(f'Write observation list to: {filename_csv}')
        export_df.to_csv(filename_csv, index=False, columns=EXPORT_OBS_LIST_COLUMNS)

    def flag_observations_based_on_obs_list_csv_file(self, obs_list_filename: str, update_type: str = 'all_obs',
                                                     verbose: bool = False):
        """ Flag observations in campaign based on an observation list file.

        Parameters
        ----------
        obs_list_filename : str
            Name and path of the input CSV file containing the observation list.
        update_type : str, optional (default = 'all')
            Defines which observations in the input list are used to update the `keep_obs` status of the matched
            observations in the campaign:. There are 3 options: (1) `all _obs` indicates that all matched observations
            are updated, (2) `inactive_only` indicates that only inactive observations (in the list) are updated and (3)
            `active_only` indicates that only active observations in the list are updated.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        """
        if verbose:
            print(f'Flag observations based on observation list in: {obs_list_filename}')
        flag_log_str = ''

        # Read csv file:
        obs_list_df = pd.read_csv(obs_list_filename)
        if len(obs_list_df) > 0:
            # Check availability of needed columns:
            # EXPORT_OBS_LIST_COLUMNS
            invalid_cols = list(set(obs_list_df.columns) - set(EXPORT_OBS_LIST_COLUMNS))
            if len(invalid_cols) > 0:
                raise AssertionError(f'Invalid columns in the observation list csv file: {", ".join(invalid_cols)}')

            # Get filter for observations to be updated according to the "update_type":
            if update_type == 'all_obs':
                pass
            elif update_type == 'inactive_only':
                obs_list_df = obs_list_df[~obs_list_df['keep_obs']]
            elif update_type == 'active_only':
                obs_list_df = obs_list_df[obs_list_df['keep_obs']]
            else:
                raise AssertionError(f'Invalid input argument for "update_type": {update_type}')

            # Apply flagging:
            surveys = obs_list_df['survey_name'].unique().tolist()
            for survey_name in surveys:
                count_changed = 0
                count_matched = 0
                flag_log_str += f'Survey: {survey_name}\n'
                if survey_name not in self.surveys:
                    flag_log_str += '  - Does not exist in this campaign.\n'
                else:
                    num_obs_in_survey = len(obs_list_df[obs_list_df['survey_name'] == survey_name])
                    for index, row in obs_list_df[obs_list_df['survey_name'] == survey_name].iterrows():
                        obs_epoch = row['obs_epoch']
                        if sys.version_info >= (3, 7):
                            epoch_dt = dt.datetime.fromisoformat(obs_epoch)
                        else:
                            if (len(obs_epoch) == 25) & (obs_epoch[-3] == ':'):  # E.g.: '2018-10-23 05:52:01+00:00'
                                epoch_dt = dt.datetime.strptime(obs_epoch[:-3] + obs_epoch[-2:], '%Y-%m-%d %H:%M:%S%z')
                            else:
                                raise AssertionError('Unknown tim format!')
                        filter_tmp = (self.surveys[survey_name].obs_df['obs_epoch'] == epoch_dt) & (
                                    self.surveys[survey_name].obs_df['station_name'] == row['station_name'])
                        num_matched_rows = len(filter_tmp[filter_tmp])
                        if num_matched_rows > 1:
                            raise AssertionError(f'In survey {survey_name} the are {num_matched_rows} observations at the '
                                                 f'same time and station!')
                        if num_matched_rows == 1:  # OK!
                            count_matched += 1
                            if self.surveys[survey_name].obs_df.loc[filter_tmp, 'keep_obs'].bool() != row['keep_obs']:
                                self.surveys[survey_name].obs_df.loc[filter_tmp, 'keep_obs'] = row['keep_obs']
                                count_changed += 1
                    flag_log_str += f'  - Matched observations: {count_matched} of {num_obs_in_survey}\n'
                    flag_log_str += f'  - Changed "keep_obs" flags: {count_changed}\n'
                    # Activate/deactivate survey according to the keep_obs flag of their observations:

        else:
            flag_log_str = 'Empty observation list file!'
        if verbose:
            print(flag_log_str)
        return flag_log_str

    def change_campaign_name(self, name: str):
        """Change the name of the campaign.

        Parameters
        ----------
        name : str
            New campaign name. Non-empty string.
        """
        if isinstance(name, str):
            if name:  # Non-empty string
                if ' ' in name:
                    raise AssertionError(f'Blanks in the campaign name are not allowed!')
                else:
                    self.campaign_name = name
            else:
                raise AssertionError('The campaign name is emtpy!')
        else:
            raise AssertionError('The campaign name is not a string!')

    def __str__(self):
        return f'Campaign "{self.campaign_name}" with {self.number_of_surveys} surveys ' \
               f'and {self.stations.get_number_of_stations} stations.'


if __name__ == '__main__':
    """Main function, primarily for debugging and testing."""
    filename = '/pyProjects/gravtools/out/test.pkl'
    camp = Campaign.from_pkl(filename)
