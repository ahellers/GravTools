"""Classes for least-squares adjustment of gravimeter-surveys.

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

import numpy as np
from scipy import stats
import datetime as dt
import pytz
import copy
import pandas as pd

from gravtools.models.misc import time_it, get_nonunique_items

# optional imports:
try:
    import geopandas
except ImportError:
    _has_geopandas = False
else:
    _has_geopandas = True

from gravtools import settings


class LSM:
    """Base class for LSM classes.

    Attributes
    ----------
    lsm_method : str
        Defines the adjustment method. Has to b listed in :py:obj:`gravtools.settings.ADJUSTMENT_METHODS`.
    stat_df : :py:obj:`gravtools.Station.stat_df`
        The station dataframe contains all relevant station data.
    setups : dict of dicts
        The setups dictionary contains all observation data used for the adjustment. The keys of the dictionary
        are the survey names (str). The items are again keys with the following items:

        - ref_epoch_delta_t_h : datetime object
            Reference epoch for the relative reference times in the column `delta_t_h` in the `setup_df` dataframe.
            The reference epoch is determined as the epoch of the first (active) observation in this survey.
        -  ref_epoch_delta_t_campaign_h : datetime object
            Reference epoch for the relative reference times in the column `delta_t_campaign_h` in the `setup_df`
            dataframe. The reference epoch is determined as the epoch of the first (active) observation in the campaign.
        - setup_df : Pandas DataFrame
            Pandas dataframes containing the observation data (see :py:obj:`gravtools.Survey.setup_df`).
        - tide_correction_type : str, optional (default='')
            Type of the tidal corrections applied on the reduced observations (column `g_red_mugal` in `obs_df`) that
            were used to calculate the setup data store in this dict. Valid entries have to be listed in
            :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`.
        - reference_height_type : str, optional (default='')
            Reference height type of the reduced observations (column `g_red_mugal` in `obs_df`) that
            were used to calculate the setup data store in this dict. Valid entries have to be listed in
            :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`.

    comment : str, optional (default = '')
        Optional comment on the adjustment run.
    init_time : datetime object
        Representing the date and time the LSM object has been initialized.
    write_log : bool, optional (default=True)
        Flag that indicates whether log string should be written or not.
    log_str : str
        String to log the status of the adjustment.
    setup_obs_df : Pandas DataFrame
        Pandas Dataframes for logging (differential or absolute) setup observations, the related metadata, estimation
        results and statistics. The columns of the dataframe may differ between adjustment methods.
    number_of_iterations : int (default=0)
        Indicates the number of iterations if iterative adjustment was applied. `0` indicated that a non-iterative
        adjustment was applied.
    drift_ref_epoch_type : string ('survey' or 'campaign'), optional (default='survey')
        Defines whether the reference epoch t0 for the estimation of the drift polynomials for each survey in the
        campaign is the reference epoch of the first (active) observation in each survey (option: 'survey') or the
        first (active) observation in the whole campaign (option: 'campaign').
        """

    def __init__(self, lsm_method, stat_df, setups, comment='', write_log=True):
        """
        Parameters
        ----------
        lsm_method : str
            Defines the adjustment method. Has to b listed in :py:obj:`gravtools.settings.ADJUSTMENT_METHODS`.
        stat_df : :py:obj:`gravtools.Station.stat_df`
            The station dataframe contains all relevant station data.
        setups : dict of dicts
            The setups dictionary contains all observation data used for the adjustment. The keys of the dictionary
            are the survey names (str). The items are again keys with the following items:

            - ref_epoch_delta_t_h : datetime object
                Reference epoch for the relative reference times in the column `delta_t_h` in the `setup_df` dataframe.
                The reference epoch is determined as the epoch of the first (active) observation in this survey.
            -  ref_epoch_delta_t_campaign_h : datetime object
                Reference epoch for the relative reference times in the column `delta_t_campaign_h` in the `setup_df`
                dataframe. The reference epoch is determined as the epoch of the first (active) observation in the campaign.
            - setup_df : Pandas DataFrame
                Pandas dataframes containing the observation data (see :py:obj:`gravtools.Survey.setup_df`).
            - tide_correction_type : str, optional (default='')
                Type of the tidal corrections applied on the reduced observations (column `g_red_mugal` in `obs_df`) that
                were used to calculate the setup data store in this dict. Valid entries have to be listed in
                :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`.
        - red_reference_height_type : str, optional (default='')
                Reference height type of the reduced observations (column `g_red_mugal` in `obs_df`) that
                were used to calculate the setup data store in this dict. Valid entries have to be listed in
                :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`.

        comment : str, optional (default = '')
            Optional comment on the adjustment run.
        write_log : bool, optional (default=True)
            Flag that indicates whether log string should be written or not.

        Notes
        -----
        In order to prevent altering the original input data (in case of assignment via reference), deep copies of the
        input data are stored in objects of this class.
        """
        # Initial input checks:
        if lsm_method in settings.ADJUSTMENT_METHODS.keys():
            self.lsm_method = lsm_method
        else:
            raise ValueError(f'"{lsm_method}" is not a valid adjustment method identifier! All valid methods are '
                             f'defined in gravtools.settings.ADJUSTMENT_METHODS.')
        if isinstance(comment, str):
            self.comment = comment
        else:
            raise ValueError(f'"{comment}" need to be a string.')
        if isinstance(write_log, bool):
            self.write_log = write_log
        else:
            raise ValueError(f'"{write_log}" needs to be a boolean.')
        self.init_time = dt.datetime.now(tz=pytz.utc)
        self.log_str = ''

        # Create deep copies of the input data items:
        self.stat_df = stat_df.copy(deep=True)
        self.setups = copy.deepcopy(setups)

        # Initialize attributes:
        self.setup_obs_df = None  # Dataframe that contain the observation (setup) related results
        self.observed_stations = None  # Unique list of observed stations; defines the station IDs for matrices
        self.stat_obs_df = None  # Station dataframe that contains estimation results of observed stations
        self.drift_pol_df = None  # DataFrame that contains the estimated parameters of the drift polynomials and their statistics for each survey
        self.vg_pol_df = None  # DataFrame that contains the estimated parameters of the vertical gravity gradient polynomial

        # Estimation settings:
        self.drift_polynomial_degree = None
        self.sig0_a_priori = None  # A priori standard deviation of unit weight
        self.scaling_factor_datum_observations = None
        self.confidence_level_chi_test = None
        self.confidence_level_tau_test = None
        self.drift_ref_epoch_type = ''  # 'survey' or 'campaign'
        self.vg_polynomial_degree = None
        self.vg_polynomial_ref_height_offset_m = 0.0

        # General statistics:
        self.number_of_stations = None
        self.number_of_datum_stations = None
        self.number_of_estimates = None
        self.degree_of_freedom = None

        # A posteriori statistics:
        self.s02_a_posteriori = None  # A posteriori variance of unit weight

        # Matrices:
        self.Cxx = None  # Co-variance matrix of estimated parameters
        self.Rxx = None  # Correlation matrix of estimates parameters
        self.x_estimate_names = []  # List of items in the Cxx matrix (index 0 to n)
        self.mat_A = None  # Design matrix
        self.mat_x = None  # Estimated parameters

        # Iterative adjustment:
        self.number_of_iterations = 0  # `0` indicates no iterations.

        # Statistical tests:
        self.number_of_outliers = None
        self.global_model_test_status = ''  # str

    @classmethod
    def from_campaign(cls, campaign, comment='', write_log=True):
        """Constructor that generates and populates the LSM object (child!) from a Campaign class object.

        Notes
        -----
        Put all checks and data preparations relevant for all LSM methods into this method!

        Parameters
        ----------
        campaign : :py:obj:`gravtools.models.survey.Campaign`
            The campaign object needs to provide setup data for all active surveys and the related station data.
        comment : str, optional (default = '')
            Arbitrary comment on the LSM run.
        write_log : bool, optional (default=True)
            Flag that indicates whether log string should be written or not.

        Returns
        -------
        LSM Object of the child class
            Contains all information required for adjusting the campaign.
        """
        # Check if all required data is available in the campaign object:

        # Comment:
        if not isinstance(comment, str):
            raise TypeError(f'"comment" needs to be a string!')

        # Station data:
        if campaign.stations is None:
            raise AssertionError(f'The campaign "{campaign.campaign_name}" does not contain any station data!')
        else:
            if not hasattr(campaign.stations, 'stat_df'):
                raise AssertionError(f'The campaign "{campaign.campaign_name}" has not station dataframe!')
            else:
                if len(campaign.stations.stat_df) == 0:
                    raise AssertionError(f'The campaign "{campaign.campaign_name}" has an empty station dataframe!')

        # Survey data:
        if campaign.surveys is None:
            raise AssertionError(f'The campaign "{campaign.campaign_name}" does not contain any survey data!')
        else:
            if len(campaign.surveys) == 0:
                raise AssertionError(f'The campaign "{campaign.campaign_name}" contains no survey data!')

        # Create setups dict:
        # Loop over surveys in campaign:
        setups = {}
        for survey_name, survey in campaign.surveys.items():
            if survey.keep_survey and survey.is_active:
                if survey.setup_df is None:
                    raise AssertionError(f'Setup data is missing for survey "{survey_name}"')
                else:
                    if len(survey.setup_df) < 2:
                        raise AssertionError(f'Survey "{survey_name}" has less than two setup observations! '
                                             f'A minimum of two setups is required in order to define '
                                             f'differential observations!')
                    else:
                        setup_data_dict = {'ref_epoch_delta_t_h': survey.ref_delta_t_dt,
                                           'ref_epoch_delta_t_campaign_h': campaign.ref_delta_t_dt,
                                           'setup_df': survey.setup_df,
                                           'tide_correction_type': survey.setup_tide_correction_type,
                                           'reference_height_type': survey.setup_reference_height_type,
                                           'atm_pres_correction_type': survey.setup_atm_pres_correction_type,
                                           'scale_correction_type': survey.setup_scale_correction_type,
                                           'setup_calc_method': survey.setup_calc_method,
                                           'setup_obs_list_df': survey.setup_obs_list_df}
                        setups[survey_name] = setup_data_dict

        # Check if setup data is available:
        if len(setups) == 0:
            raise AssertionError(f'Setup data of campaign "{campaign.campaign_name}" does not contain observations!')

        # Initialize and return LSM object:
        return cls(campaign.stations.stat_df, setups, comment=comment, write_log=write_log)

    def adjust_autoscale_s0(self,
                            iteration_approach='Multiplicative',
                            s02_target=1,
                            s02_target_delta=0.1,
                            max_number_iterations=10,
                            add_const_to_sd_of_observations_step_size_mugal=5.0,
                            max_total_additive_const_to_sd_mugal=20.0,
                            multiplicative_factor_step_size_percent=10.0,  # [%]
                            max_multiplicative_factor_to_sd_percent=200.0,  # [%]
                            min_multiplicative_factor_to_sd_percent=50.0,  # [%]
                            drift_pol_degree=1,
                            sig0_mugal=1,
                            scaling_factor_datum_observations=1.0,
                            add_const_to_sd_of_observations_mugal=0.0,
                            scaling_factor_for_sd_of_observations=1.0,
                            confidence_level_chi_test=0.95,
                            confidence_level_tau_test=0.95,
                            drift_ref_epoch_type='survey',
                            noise_floor_mugal=0.0,
                            verbose=False,
                            ):  # Confidence level):
        """Run the adjustment iteratively in order to adjust s0 to the target value by adapting the SD of observations.

        Notes
        -----
        Iterative adjustment works for differential and non-differential LSM adjustment, not for VG estimation!

        Parameters
        ---------
        iteration_approach : str, optional (default='Multiplicative')
            Defines the iteration approach (Multiplicative od Additive)
        s02_target : float, optional (default=1.0)
            Target a posteriori s0² (variance of unit weight) for the iteration.
        s02_target_delta : float, optional (default=0.1)
            Permissible deviation of the target a posteriori s0². As soon as the a posteriori s0 of the current lsm run
            lies withing the threshold of `s02_target` +- `s02_target_delta` the iteration was successful and stops.
        max_number_iterations : int, optional (default=10)
            Maximum allowed number of iterations. If the a posteriori s0 does not lie within the defined threshold
            (`s02_target` +- `s02_target_delta`) after `max_number_iterations` iterations an assertion error is raised.
        add_const_to_sd_of_observations_step_size_mugal : float, optional (default=5.0)
            Initial iteration step size for the additive constant that is added to the standard deviation of all
            observations. This parameter is only considered when using the `additive` iteration approach.
        max_total_additive_const_to_sd_mugal : float, optional (default=20.0)
            If the additive factor for the SD of observations that is determined iteratively in order to reach the
            defined target s0 is larger than `max_total_additive_const_to_sd_mugal`, the iteration procedure
            failed and an assertion error is raised. This parameter is only considered when using the `additive`
            iteration approach.
        multiplicative_factor_step_size_percent : float, optional (default=10.0)
            Initial iteration step size for the multiplicative factor that is used to scale all setup
            observations. This parameter is only considered when using the `multiplicative` iteration approach.
        max_multiplicative_factor_to_sd_percent : float, optional (default=200.0)
            Maximum scaling factor when using the `multiplicative` iteration approach for scaling the SD of setup
            observations. Minimum = 100%.
        min_multiplicative_factor_to_sd_percent : float, optional (default=50.0)
            Minimum scaling factor when using the `multiplicative` iteration approach for scaling the SD of setup
            observations. Minimum = 1.0%, maximum = 100.0%
        drift_pol_degree : int, optional (default=1)
            Degree of estimated drift polynomial.
        sig0_mugal : int, optional (default=1)
            A priori standard deviation of unit weight of observations [µGal] for the stochastic model of the
            least-squares adjustment.
        scaling_factor_datum_observations : float, optional (default=1.0)
            Factor for scaling the standard deviation (SD) of g of datum stations. The scaled SD is is used for
            weighting the direct pseudo observations of g at the datum stations that are introduced as datum
            constraints.
        add_const_to_sd_of_observations_mugal : float, optional (default=0.0)
            The defined additive constant is added to the standard deviation (SD) of setup
            observations in order to scale the SD and the resulting weights to realistic values. In µGal. The scaling
            factor `scaling_factor_for_sd_of_observations` is applied before adding this constant!
        scaling_factor_for_sd_of_observations : float, optional (default=1.0)
            Scaling factor for the standard deviation of the setup observations. `add_const_to_sd_of_observations_mugal`
            is applied after applying the scaling factor!
        confidence_level_chi_test : float, optional (default=0.95)
            Confidence level for the goodness-of-fit test.
        confidence_level_tau_test : float, optional (default=0.95)
            Confidence level for the tau test.
        drift_ref_epoch_type : string ('survey' or 'campaign'), optional (default='survey')
            Defines whether the reference epoch t0 for the estimation of the drift polynomials for each survey in the
            campaign is the reference epoch of the first (active) observation in each survey (option: 'survey') or the
            first (active) observation in the whole campaign (option: 'campaign').
        noise_floor_mugal : float, optional (default=0.0)
            The standard error SE of the estimated gravity values at stations is calculated by SE = sqrt(SD**2 + NF**2),
            where SD is the estimated standard deviation of the station's gravity and NF is the noise floor value
            specified here.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line, e.g. for debugging and testing
        """

        # Init.:
        flag_s0_within_threshold = False
        # scale_factor = scaling_factor_for_sd_of_observations
        add_const = add_const_to_sd_of_observations_mugal  # initial
        add_const_step_size = add_const_to_sd_of_observations_step_size_mugal  # [mugal]
        mult_factor = scaling_factor_for_sd_of_observations  # initial
        mult_factor_step_size = multiplicative_factor_step_size_percent / 100  # [%] => factor
        last_iteration_step_action = ''  # 'increase' or 'decrease'
        iteration_log_str = ''
        complete_log_str = ''
        add_const_total_mugal = 0
        mult_factor_total = 1.0

        # Run iteration:
        for i_iteration in range(1, max_number_iterations + 1):

            # Only apply multiplicative facot or additive constant in the first iteration, depending on the iteration
            # approach:
            if i_iteration > 1:
                if iteration_approach == 'Additive':
                    mult_factor = 1.0
                elif iteration_approach == 'Multiplicative':
                    add_const = 0.0

            add_const_total_mugal += add_const  # Log total additive constant
            mult_factor_total *= mult_factor  # Log total multiplicative factor
            self.log_str = ''  # Reset log string of the previous run

            self.adjust(drift_pol_degree=drift_pol_degree,
                        sig0_mugal=sig0_mugal,
                        scaling_factor_datum_observations=scaling_factor_datum_observations,
                        add_const_to_sd_of_observations_mugal=add_const,  # Adjusted iteratively
                        scaling_factor_for_sd_of_observations=mult_factor,  # Adjusted iteratively
                        confidence_level_chi_test=confidence_level_chi_test,
                        confidence_level_tau_test=confidence_level_tau_test,
                        drift_ref_epoch_type=drift_ref_epoch_type,
                        noise_floor_mugal=noise_floor_mugal,
                        verbose=False
                        )

            iteration_log_str_tmp = f'########## Iteration {i_iteration} #########################\n'
            # if iteration_approach == 'Additive':
            iteration_log_str_tmp += f'Total additive const. to SD: {add_const_total_mugal:5.3f} µGal\n'
            iteration_log_str_tmp += f'Current additive const. to SD: {add_const:5.3f} µGal\n'
            # elif iteration_approach == 'Multiplicative':
            iteration_log_str_tmp += f'Total mult. factor for SD: {mult_factor_total:5.3f}\n'
            iteration_log_str_tmp += f'Current mult. factor for SD: {mult_factor:5.3f}\n'
            iteration_log_str_tmp += f's0² a posteriori: {self.s02_a_posteriori:5.3f} µGal²\n'
            iteration_log_str_tmp += f'\n'
            if verbose:
                print(iteration_log_str_tmp)
                print(self.log_str)
            if self.write_log:
                complete_log_str = complete_log_str + iteration_log_str_tmp + self.log_str + '\n\n'
            iteration_log_str += iteration_log_str_tmp

            if (self.s02_a_posteriori < (s02_target + s02_target_delta)) and (
                    self.s02_a_posteriori > (s02_target - s02_target_delta)):
                flag_s0_within_threshold = True
            elif self.s02_a_posteriori > (s02_target + s02_target_delta):  # Too large => Increase SD of obs.
                if last_iteration_step_action == 'decrease':
                    add_const_step_size = add_const_step_size / 2
                    mult_factor_step_size = mult_factor_step_size / 2
                last_iteration_step_action = 'increase'
                add_const = add_const_step_size
                mult_factor = 1.0 + mult_factor_step_size
            elif self.s02_a_posteriori < (s02_target - s02_target_delta):  # Too small => Decrease SD of obs.
                if last_iteration_step_action == 'increase':
                    add_const_step_size = add_const_step_size / 2
                    mult_factor_step_size = mult_factor_step_size / 2
                last_iteration_step_action = 'decrease'
                add_const = -add_const_step_size
                mult_factor = 1.0 - mult_factor_step_size

            if flag_s0_within_threshold:  # Exit loop and stop iteration
                break

        # Write status message to log according to iteration result:
        if iteration_approach == 'Additive':
            if flag_s0_within_threshold and (add_const_total_mugal <= max_total_additive_const_to_sd_mugal):
                iteration_log_str_tmp = f' => Iteration successful!\n'
                iteration_log_str_tmp += f' => s0² a posteriori of {self.s02_a_posteriori:1.3f} within [{s02_target - s02_target_delta:1.3f}, {s02_target + s02_target_delta:1.3f}]\n'
                iteration_log_str_tmp += f' => Total additive constant to SD of observations ({add_const_total_mugal:1.3f}) ' + \
                                         f'is smaller than the user defined threshold ' + \
                                         f'of {max_total_additive_const_to_sd_mugal:1.3f} µGal.\n'
            else:
                iteration_log_str_tmp = f' => ERROR: Iteration failed!\n'
                if not flag_s0_within_threshold:
                    iteration_log_str_tmp += f' => s0² a posteriori of {self.s02_a_posteriori:1.3f} not within ' + \
                                             f'[{s02_target - s02_target_delta:1.3f}, {s02_target + s02_target_delta:1.3f}]\n'
                if (add_const_total_mugal > max_total_additive_const_to_sd_mugal):
                    iteration_log_str_tmp += f' => Total additive constant to SD  ({add_const_total_mugal:1.3f}) ' + \
                                             f'exceeds the the user defined threshold ' + \
                                             f'of {max_total_additive_const_to_sd_mugal:1.3f} µGal.\n'
            iteration_log_str_tmp += f'\n'

        elif iteration_approach == 'Multiplicative':
            if flag_s0_within_threshold and (mult_factor_total <= (max_multiplicative_factor_to_sd_percent / 100)) and (
                    mult_factor_total >= (min_multiplicative_factor_to_sd_percent / 100)):
                iteration_log_str_tmp = f' => Iteration successful!\n'
                iteration_log_str_tmp += f' => s0² a posteriori of {self.s02_a_posteriori:1.3f} within [{s02_target - s02_target_delta:1.3f}, {s02_target + s02_target_delta:1.3f}]\n '
                iteration_log_str_tmp += f' => Total multiplicative factor for SD of observations ({mult_factor_total * 100:1.3f}%) ' + \
                                         f'is between the user defined threshold ' + \
                                         f'of {min_multiplicative_factor_to_sd_percent:1.3f}% and {max_multiplicative_factor_to_sd_percent:1.3f}%.\n'
            else:
                iteration_log_str_tmp = f' => ERROR: Iteration failed!\n'
                if not flag_s0_within_threshold:
                    iteration_log_str_tmp += f' => s0² a posteriori of {self.s02_a_posteriori:1.3f} not within ' + \
                                             f'[{s02_target - s02_target_delta:1.3f}, {s02_target + s02_target_delta:1.3f}]\n'
                if (mult_factor_total > (max_multiplicative_factor_to_sd_percent / 100)):
                    iteration_log_str_tmp += f' => Total multiplicative factor for SD  ({mult_factor_total * 100:1.3f}%) ' + \
                                             f'exceeds the the user defined threshold ' + \
                                             f'of {max_multiplicative_factor_to_sd_percent:1.3f}%.\n'
                if (mult_factor_total < (min_multiplicative_factor_to_sd_percent / 100)):
                    iteration_log_str_tmp += f' => Total multiplicative factor for SD  ({mult_factor_total * 100:1.3f}%) ' + \
                                             f'is smaller than the user defined threshold ' + \
                                             f'of {min_multiplicative_factor_to_sd_percent:1.3f}%.\n'
            iteration_log_str_tmp += f'\n'

        iteration_log_str += iteration_log_str_tmp

        if verbose:
            print(iteration_log_str_tmp)

        # Append iteration log to the log string of the last iteration:
        if self.write_log:
            self.log_str = complete_log_str + \
                           '\n########## Iteration log #########################\n\n' + \
                           f' - Iteration approach: {iteration_approach}\n\n' + \
                           iteration_log_str

        if not flag_s0_within_threshold:
            raise AssertionError(f'Iteration Error: s0² a posteriori of {self.s02_a_posteriori:1.3f} not within '
                                 f'[{s02_target - s02_target_delta:1.3f}, {s02_target + s02_target_delta:1.3f}] after '
                                 f'{i_iteration} iterations!')

        if iteration_approach == 'Additive':
            if abs(add_const_total_mugal) > max_total_additive_const_to_sd_mugal:
                raise AssertionError(f'Total additive constant to SD  ({add_const_total_mugal:1.3f}) '
                                     f'exceeds the the user defined threshold '
                                     f'of {max_total_additive_const_to_sd_mugal:1.3f} µGal.\n')
        elif iteration_approach == 'Multiplicative':
            if mult_factor_total > (max_multiplicative_factor_to_sd_percent / 100):
                raise AssertionError(f'Total multiplicative factor for SD  ({mult_factor_total * 100:1.3f}%) '
                                     f'exceeds the the user defined threshold '
                                     f'of {max_multiplicative_factor_to_sd_percent:1.3f}%.\n')
            if mult_factor_total < (min_multiplicative_factor_to_sd_percent / 100):
                raise AssertionError(f'Total multiplicative factor for SD  ({mult_factor_total * 100:1.3f}%) '
                                     f'is smaller than the the user defined threshold '
                                     f'of {min_multiplicative_factor_to_sd_percent:1.3f}%.\n')

        self.number_of_iterations = i_iteration

    def export_stat_results_shapefile(self, filename, epsg_code, verbose=True):
        """Export station-related results from the `stat_obs_df` dataframe to a shapefile.

        Notes
        -----
        This method relies on the optional module `geopandas` (optional dependency).

        Parameters
        ---------
        filename : str
            Name and path of the output shapefile.
        epsg_code: int
            EPSG code of the station coordinates' CRS.
        verbose : bool, optional (default=False)
            `True` implies that status messages are printed to the command line.
        """
        if not _has_geopandas:
            raise ImportError(f'Optional dependency geopandas not available, but needed for writing shapefiles!')
        if verbose:
            print(f'Save station results of lsm run "{self.comment}" to: {filename}')

        stat_obs_df = self.stat_obs_df.copy(deep=True)
        stat_obs_df['comment'] = self.comment
        stat_obs_gdf = geopandas.GeoDataFrame(stat_obs_df, geometry=geopandas.points_from_xy(stat_obs_df['long_deg'],
                                                                                             stat_obs_df['lat_deg']),
                                              crs=epsg_code)
        stat_obs_gdf.to_file(filename)

    @property
    def time_str(self):
        """Return time of lsm adjustment as formatted string."""
        return self.init_time.strftime("%Y-%m-%d, %H:%M:%S")

    @property
    def get_log_string(self):
        """Returns the log string."""
        return self.log_str

    @property
    # @time_it
    def get_correlation_matrix(self):
        """Calculates und returns the correlation matrix based on the Co-Variance matrix of the LSM run.

        Returns
        -------
        numpy array : Correlation matrix (n, n)

        Notes
        -----
        The co-variance matrix `self.Cxx` needs to exist.
        """
        if self.Cxx is None:
            # raise AssertionError('Correlation matrix cannot be calculated due to missing co-variance matrix.')
            return None
        mat_Rxx = np.full((self.Cxx.shape[0], self.Cxx.shape[1]), np.nan)
        for i_row in range(0, mat_Rxx.shape[0]):
            cxx_row_sqrt = np.sqrt(self.Cxx[i_row, i_row])
            for i_col in range(0, i_row + 1):
                mat_Rxx[i_row, i_col] = self.Cxx[i_row, i_col] / (
                        cxx_row_sqrt * np.sqrt(self.Cxx[i_col, i_col]))
                mat_Rxx[i_col, i_row] = mat_Rxx[i_row, i_col]
        return mat_Rxx

    @property
    def get_results_obs_df(self):
        """Getter for the observation-related results."""
        return self.setup_obs_df

    @property
    def get_results_drift_df(self):
        """Getter for the drift-related results."""
        return self.drift_pol_df

    @property
    def get_results_stat_df(self):
        """Getter for the station-related results."""
        return self.stat_obs_df

    @property
    def get_results_vg_df(self):
        """Getter for the results of vertical gravity gradient estimation."""
        try:
            return self.vg_pol_df
        except AttributeError:  # If "self" has no attribute "vg_pol_df" (not initialized)
            return None

    @property
    def check_for_unique_setup_id(self):
        """Indicates whether setup IDs should be unique based on the setup data calculation method.

        Returns
        -------
        bool : `True`, if the setup IDs should be unique, otherwise `False`.

        Notes
        -----
        If observations were not aggregated within at least one setup, checking for unique setup IDs makes no sense.
        """
        for survey_name, setup_data in self.setups.items():
            if setup_data['setup_calc_method'] == 'individual_obs':
                return False
        return True

    def check_unique_setups(self, setup_ids: list):
        """Raise an error, is the setup IDs in the input list are not unique.

        Notes
        -----
        Omit the check, if `check_for_unique_setup_id` returns `False`.
        """
        if self.check_for_unique_setup_id:
            if len(set(setup_ids)) != len(setup_ids):
                non_unique_items = get_nonunique_items(setup_ids)
                non_unique_items_str = [str(item) for item in non_unique_items]
                non_unique_items_str = ", ".join(non_unique_items_str)
                raise RuntimeError(
                    f'Setup IDs are not unique within the campaign! Setups with the following IDs are not unique: {non_unique_items_str}')

    @property
    def survey_info_string(self) -> str:
        """Returns information on all surveys."""
        survey_names_list = []
        num_setups_list = []
        num_obs_total_list = []
        num_obs_active_list = []
        tide_correction_type_list = []
        reference_height_type_list = []
        atm_pres_correction_type_list = []
        scale_correction_type_list = []
        setup_calc_method_list = []
        for survey_name, survey in self.setups.items():
            survey_names_list.append(survey_name)
            setup_df = survey['setup_df']
            num_setups_list.append(len(setup_df))
            setup_obs_list_df = survey['setup_obs_list_df']
            num_obs_total_list.append(len(setup_obs_list_df))
            num_obs_active_list.append(len(setup_obs_list_df.loc[setup_obs_list_df['keep_obs']]))
            tide_correction_type_list.append(survey['tide_correction_type'])
            reference_height_type_list.append(survey['reference_height_type'])
            atm_pres_correction_type_list.append(survey['atm_pres_correction_type'])
            scale_correction_type_list.append(survey['scale_correction_type'])
            setup_calc_method_list.append(survey['setup_calc_method'])
        tmp_df = pd.DataFrame(list(zip(survey_names_list,
                                       num_setups_list,
                                       num_obs_total_list,
                                       num_obs_active_list,
                                       tide_correction_type_list,
                                       reference_height_type_list,
                                       atm_pres_correction_type_list,
                                       scale_correction_type_list,
                                       setup_calc_method_list,
                                       )),
                              columns=['Survey', '#Setups', '#Obs', '#Active obs', 'Tide correction', 'Ref. height',
                                       'Atm. correction', 'Scale correction', 'Setup calc. method'])
        return tmp_df.to_string(index=False)

    # @property
    # def obs_level_corrections_info_string(self) -> str:
    #     """Returns information on the applied observation-level corrections."""
    #     info_str = ''
    #     return info_str


def bin_redundancy_components(mat_r):
    """Bin redundancy components for easier interpretatzion.

    See: Skriptum AG2 (Navratil, TU Wien), p. 70.
    """
    number_obs_r_equal_0 = 0  # No error control: Errors cannot be detected
    number_obs_r_0_to_03 = 0  # Bad error control
    number_obs_r_03_to_07 = 0  # Good error control
    number_obs_r_07_1 = 0  # Very good error control, but observation may be redundant
    number_obs_r_equal_1 = 0  # Observation is redundant

    results_dict = {}  # '': ''

    return number_obs_r_equal_0, number_obs_r_0_to_03


def tau_test(mat_w, dof, alpha, mat_r):
    """Tau-criterion test for outlier detection of a least-squares adjustment.

    This test considers that the true variance factor sigma_0² is unknown. It is based on normalized residuals
    (Studentized residuals) and the Tau distribution.

    If the test fails it indicates that ONE (!) observation is a gross error, i.e. ONE residual is an outlier!
    If the adjustment is affected by two or more gross errors this test is not applicable. In that case the
    following pragmatic approach is recommended (e.g. by Caspary, 1987): The observation with the largest test
    statistic is discarded. Then the adjustment is repeated with the remaining n-1 observations and the test is
    repeated based on the new results, etc.

    Notes
    -----
    Ref.: Pope (1976): The statistics of residuals and the detection of outliers, NOAA Technical Report NOS 65 NGS 1
    Ref.: Caspari (1987): Concepts of network and deformation analysis. Monograph by UNSW Sydney, pp. 76-77
    This method implements a two-sided test, as the critical value is calculated with (alpha/2) and the abs. value of
    the normalized residual w is tested.

    Parameters
    ----------
    mat_w: np.array of floats
        Vector with standardized post-fit residuals of all observations in the least-squares adjustment.
        They are computed as post-fit residuals divides by their standard deviations.
    dof: int
        Degree of freedom of the least-squares adjustment.
    alpha: float (0 to 1)
        Significance level.
    mat_r: np.array of floats
        Vector with redundancy components derived in the least-squares adjustment for each observation.
    """
    # Critical value (two-sided test: the abs. value of w is tested und the Tau-distribution is symmetric about 0!):
    tsd_crt = stats.t.ppf(1 - alpha / 2, dof - 1)  # t distribution
    tau_crt = (tsd_crt * np.sqrt(dof)) / np.sqrt(dof - 1 + tsd_crt ** 2)  # Critical value (Pope, 1976, Eq. (6))
    tau_test_result = []
    for idx, w in enumerate(mat_w):
        # Check, if redancy component is larger than threshold:
        if mat_r[idx] > settings.R_POPE_TEST_THRESHOLD:
            if np.abs(w) > tau_crt:  # Equ. (6-36) in Caspary (1987)
                tau_test_result.append('failed')
            else:
                tau_test_result.append('passed')
        else:
            tau_test_result.append(f'r too small')  # Pope test not applied due to small redundancy component!
    return tau_test_result, tau_crt


def create_hist(mat_v):
    """Create histogram."""
    residuals = np.ndarray.tolist(mat_v)
    residuals = [item for sublist in residuals for item in sublist]
    hist_residuals, bin_edges = np.histogram(residuals, bins=5)
    return hist_residuals, bin_edges

def global_model_test(cf, dof, a_posteriori_variance_of_unit_weight, a_priori_variance_of_unit_weight):
    """Global model test based on the comparison of a posteriori and a priori variance of unit weight.

    This global model test asserts if the "model is correct and complete". If it fails it indicates that the
    observations contradict the mathematical adjustment model. However, this test cannot the validity of the model or
    the correctness of the observation!

    Notes
    -----
    This "global model test" is described by Caspari (1987): Concepts of network and deformation analysis. Monograph
    by UNSW Sydney, pp. 6-8 and pp.68-69.
    This method implements a two-sided test.

    Parameters
    ----------
    cf: float
        Confidence level for Chi²
    dof: int
        Degree of freedom.
    a_posteriori_variance_of_unit_weight: float
        A posteriori variance of unit weight (after adjustment).
    a_priori_variance_of_unit_weight: float
        A priori variance of unit weight (before adjustment).

    Returns
    -------
    """
    alpha = 1 - cf  # Significance level = Probability of committing a type 1 error (H0 wrongly dismissed)
    chi_crit_upper = stats.chi2.ppf(1 - alpha / 2, dof)  # critical value
    chi_crit_lower = stats.chi2.ppf(alpha / 2, dof)  # critical value
    # chi_critical_value = stats.chi2.ppf(alpha, dof)  # critical value
    test_value = dof * a_posteriori_variance_of_unit_weight / a_priori_variance_of_unit_weight  # Equ. (2-13) in Caspary (1987)
    if chi_crit_lower < test_value < chi_crit_upper:
        # if test_value <= test_value:  # Equ. (2-15) in Caspary (1987)
        chi_test_status = 'Passed'
    else:
        chi_test_status = 'Not passed'
    chi_crit = [chi_crit_lower, chi_crit_upper]
    return chi_crit, test_value, chi_test_status
