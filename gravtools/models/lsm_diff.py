"""Classes for least-squares adjustment of differential relative gravimeter observations.

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
import pandas as pd
import datetime as dt
import pytz
# from matplotlib import pyplot as plt

from gravtools import settings
from gravtools.models.lsm import LSM, create_hist, goodness_of_fit_test, tau_test
from gravtools.models import misc
from gravtools import __version__ as GRAVTOOLS_VERSION


class LSMDiff(LSM):
    """Least-squares adjustment of differential gravimeter observations with weighted constraints.

    Attributes
    ----------
    setup_obs_df : Pandas DataFrame
        Pandas Dataframes for logging (differential or absolute) setup observations, the related metadata, estimation
        results and statistics. The columns of the dataframe are defined in `self._SETUP_DIFF_COLUMNS`.
    observed_stations : list of str
        Unique list of names of observed stations, defining the order of stations (station IDs) for all matrices and
        vectors used in the adjustment scheme.
    """

    # Column names of self.setup_obs_df:
    # - keys: Column names of the pandas dataframe
    # - values: Short description for table headers, etc., in the GUI
    _SETUP_DIFF_COLUMNS_DICT = {
        'survey_name': 'Survey',
        'ref_epoch_dt': 'Epoch',
        'diff_obs_id': 'Obs. ID',
        'station_name_from': 'Station from',
        'station_name_to': 'Station to',
        'setup_id_from': 'Setup ID from',
        'setup_id_to': 'Setup ID to',
        'g_diff_mugal': 'g [??Gal]',
        'sd_g_diff_mugal': 'SD [??Gal]',
        'sd_g_diff_est_mugal': 'SD_est [??Gal]',
        'v_diff_mugal': 'Residuals [??Gal]',  # Post fit residuals
        'w_diff_mugal': 'Std. Residual []',
        'r_diff_obs': 'Redundancy []',
        'tau_test_result': 'Outlier Test',
    }
    _SETUP_DIFF_COLUMNS = list(_SETUP_DIFF_COLUMNS_DICT.keys())

    # Column names of self.drift_pol_df:
    # - keys: Column names of the pandas dataframe
    # - values: Short description for table headers, etc., in the GUI
    _DRIFT_POL_DF_COLUMNS_DICT = {
        'survey_name': 'Survey',
        'degree': 'Degree',
        'coefficient': 'Coefficient',
        'sd_coeff': 'SD',
        'coeff_unit': 'Unit',
        'ref_epoch_t0_dt': 't0',
    }
    _DRIFT_POL_DF_COLUMNS = list(_DRIFT_POL_DF_COLUMNS_DICT.keys())

    def __init__(self, stat_df, setups, comment='', write_log=True):
        """
        Parameters
        ----------
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

        comment : str, optional (default = '')
            Arbitrary comment on the LSM run.
        write_log : bool, optional (default=True)
            Flag that indicates whether log string should be written or not.
        """
        # Call constructor from abstract base class:
        lsm_method = 'LSM_diff'
        super().__init__(lsm_method, stat_df=stat_df, setups=setups, comment=comment, write_log=write_log)

    @classmethod
    def from_campaign(cls, campaign, comment='', write_log=True):
        """Constructor that generates and populates the LSM object from a Campaign class object.

        Notes
        -----
        From all active surveys in the campaign the setup data (= observations) are loaded, additionally to the station
        data (station dataframe).

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
        :py:obj:`.LSMDiff`
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
            if survey.keep_survey:
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
                                           'setup_df': survey.setup_df}
                        setups[survey_name] = setup_data_dict

        # Check if setup data is available:
        if len(setups) == 0:
            raise AssertionError(f'Setup data of campaign "{campaign.campaign_name}" does not contain observations!')

        # Initialize and return LSM object:
        return cls(campaign.stations.stat_df, setups, comment=comment, write_log=write_log)

    def adjust(self, drift_pol_degree=1,
               sig0_mugal=1,
               scaling_factor_datum_observations=1.0,
               add_const_to_sd_of_observations_mugal=0.0,
               scaling_factor_for_sd_of_observations=1.0,
               confidence_level_chi_test=0.95,
               confidence_level_tau_test=0.95,
               drift_ref_epoch_type='survey',
               verbose=False
               ):  # Confidence level):
        """Run the adjustment based on differential observations.

        Parameters
        ----------
        drift_pol_degree : int, optional (default=1)
            Degree of estimated drift polynomial.
        sig0_mugal : int, optional (default=1)
            A priori standard deviation of unit weight of observations [??Gal] for the stochastic model of the
            least-squares adjustment.
        scaling_factor_datum_observations : float, optional (default=1.0)
            Factor for scaling the standard deviation (SD) of g of datum stations. The scaled SD is used for
            weighting the direct pseudo observations of g at the datum stations that are introduced as datum
            constraints.
        add_const_to_sd_of_observations_mugal : float, optional (default=0.0)
            The defined additive constant is added to the standard deviation (SD) of setup
            observations in order to scale the SD and the resulting weights to realistic values. In ??Gal. The scaling
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
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line, e.g. for debugging and testing

        Notes
        -----
        The adjustemnt method (weighted constraints adjustment based on differential observations) is described in
        Wijaya et al. (2019) - pyGABEUR-ITB: A free Software for Adjustment of Relative Gravimeter Data
        and in
        Hwang et al. (2002) - Adjustment of relative gravity measurements using weighted and datum-free constraints,
        Computers & Geosciences 28 (2002), pp. 1005-1015 (implemented in the software "gravnet").
        """

        # Prepare lists and indices:
        # - Observations and parameters:
        self.observed_stations = []
        survey_names = []
        number_of_observations = 0
        setup_ids = []
        for survey_name, setup_data in self.setups.items():
            setup_df = setup_data['setup_df']
            self.observed_stations = self.observed_stations + setup_df['station_name'].to_list()
            number_of_observations = number_of_observations + len(setup_df)
            survey_names.append(survey_name)
            setup_ids = setup_ids + setup_df['setup_id'].to_list()
        self.observed_stations = misc.unique_ordered_list(
            self.observed_stations)  # Unique list of stations => order of stations in matrices!
        number_of_stations = len(self.observed_stations)
        number_of_surveys = len(self.setups)
        number_of_parameters = number_of_stations + drift_pol_degree * number_of_surveys  # Total number of parameters to be estimated
        number_of_diff_obs = number_of_observations - number_of_surveys

        # Check, if setup IDs are unique:
        if len(set(setup_ids)) != len(setup_ids):
            raise AssertionError('Setup IDs are not unique within the campaign!')

        # - Datum points for weighted constraints:
        # Get dataframe with subset of observed stations only:
        filter_tmp = self.stat_df['station_name'].isin(self.observed_stations)
        self.stat_obs_df = self.stat_df.loc[filter_tmp].copy(deep=True)  # All observed stations
        filter_tmp = self.stat_obs_df['is_datum'] == True
        stat_df_obs_datum = self.stat_obs_df.loc[filter_tmp]
        datum_stations = stat_df_obs_datum['station_name'].to_list()
        number_of_datum_stations = len(datum_stations)
        if number_of_datum_stations < 1:
            raise AssertionError('None of the observed (and active) stations is a datum station '
                                 '(minimum one is required)!')
        if stat_df_obs_datum.loc[:, 'g_mugal'].isnull().any() or stat_df_obs_datum.loc[:, 'sd_g_mugal'].isnull().any():
            raise AssertionError('g and/or sd(g) is missing for at least one datum station! Both values are required'
                                 'for ALL datum stations.')

        if verbose or self.write_log:
            time_now_str = dt.datetime.now(tz=pytz.UTC).strftime('%Y-%m-%d, %H:%M:%S %Z')
            tmp_str = f'#### Adjustment log (differential LSM) ####\n'
            tmp_str += f'Processed with GravTools {GRAVTOOLS_VERSION} ({time_now_str})\n'
            tmp_str += f'\n'
            tmp_str += f'---- Input data and settings ----\n'
            tmp_str += f'Number of surveys: {number_of_surveys}\n'
            tmp_str += f'Number of stations: {number_of_stations}\n'
            tmp_str += f'Number of differential observations: {number_of_diff_obs}\n'
            tmp_str += f'Number of estimated parameters: {number_of_parameters}\n'
            tmp_str += f'Number of datum stations: {number_of_datum_stations}\n'
            tmp_str += f'Degree of freedom (w/o datum constraints): {number_of_diff_obs - number_of_parameters}\n'
            tmp_str += f'Degree of freedom (with datum constraints): {number_of_diff_obs - number_of_parameters + number_of_datum_stations}\n'
            tmp_str += f'\n'
            tmp_str += f'Degree of drift polynomial: {drift_pol_degree}\n'
            tmp_str += f'One reference epoch for each: {drift_ref_epoch_type}\n'
            tmp_str += f'A priori std. deviation of unit weight [??Gal]: {sig0_mugal}\n'
            tmp_str += f'Scaling factor for datum constraints: {scaling_factor_datum_observations}\n'
            tmp_str += f'Scaling factor for SD of setup observations: {scaling_factor_for_sd_of_observations}\n'
            tmp_str += f'Additive const. to SD of setup obs. [??Gal]: {add_const_to_sd_of_observations_mugal}\n'
            tmp_str += f'Confidence level Chi-test: {confidence_level_chi_test:4.2f}\n'
            tmp_str += f'Confidence level Tau-test: {confidence_level_tau_test:4.2f}\n'
            tmp_str += f'\n'
            tmp_str += f'---- Set up and populate matrices ----\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Initialize matrices:
        # => Initialize complete matrices first and then populate them. This is most efficient!
        # - Observation model:
        mat_A0 = np.zeros([number_of_diff_obs, number_of_parameters])  # Model-matrix
        mat_L0 = np.zeros((number_of_diff_obs, 1))
        mat_sig_ll0 = np.zeros(number_of_diff_obs)
        # - Constraints:
        mat_Ac = np.zeros([number_of_datum_stations, number_of_parameters])  # Model-matrix for constraints
        mat_Lc = np.zeros((number_of_datum_stations, 1))
        #  => Convert to diagonal matrix: P0 = np.diag(np.array([1,2,3,4])) = np.diag(mat_p0)
        mat_sig_llc = np.zeros(number_of_datum_stations)

        # Populate matrices:
        diff_obs_id = -1  # Index of differential observations in vectors mat_L0, rows of mat_A0 and mat_p0
        pd_drift_col_offset = number_of_stations - 1  # Column offset for drift parameters in A-matrix
        survey_count = -1  # Survey counter for indexing the drift parameters in the A-matrix
        g_diff_obs_mugal_list = []
        station_name_from_list = []
        station_name_to_list = []
        setup_id_from_list = []
        setup_id_to_list = []
        sd_g_diff_obs_mugal_list = []
        diff_obs_id_list = []
        survey_names_list = []
        ref_epoch_dt_list = []  # Reference epochs (datetime objects) of differential observations

        # # Scale and manipulate the SD of setup observations in order to adjust their weights in the adjustment:
        # # - Apply scaling factor to SD of setup observations:
        # setup_df['sd_g_mugal'] = setup_df['sd_g_mugal'] * scaling_factor_for_sd_of_observations
        # # - Add additive constant to standard deviation of setup observations in order to scale them to realistic values:
        # setup_df['sd_g_mugal'] = setup_df['sd_g_mugal'] + add_const_to_sd_of_observations_mugal

        for survey_name, setup_data in self.setups.items():
            setup_df = setup_data['setup_df']
            previous_row = None
            survey_count += 1

            # Scale and manipulate the SD of setup observations in order to adjust their weights in the adjustment:
            # - Apply scaling factor to SD of setup observations:
            setup_df['sd_g_mugal'] = setup_df['sd_g_mugal'] * scaling_factor_for_sd_of_observations
            # - Add additive constant to standard deviation of setup observations in order to scale them to realistic values:
            setup_df['sd_g_mugal'] = setup_df['sd_g_mugal'] + add_const_to_sd_of_observations_mugal

            for index, row in setup_df.iterrows():
                if previous_row is None:  # First time the loop is entered
                    previous_row = row
                else:
                    diff_obs_id += 1  # Increment diff. observation ID
                    g_diff_mugal = row['g_mugal'] - previous_row['g_mugal']
                    sd_g_diff_mugal = np.sqrt(row['sd_g_mugal'] ** 2 + previous_row['sd_g_mugal'] ** 2)
                    station_name_from = previous_row['station_name']
                    station_name_to = row['station_name']
                    setup_id_from = previous_row['setup_id']
                    setup_id_to = row['setup_id']
                    # epoch_from = previous_row['epoch_unix']  # [sec]
                    if drift_ref_epoch_type == 'survey':
                        epoch_from = previous_row['delta_t_h']  # [hours]
                    elif drift_ref_epoch_type == 'campaign':
                        epoch_from = previous_row['delta_t_campaign_h']  # [hours]
                    # epoch_to = row['epoch_unix']  # [sec]
                    if drift_ref_epoch_type == 'survey':
                        epoch_to = row['delta_t_h']  # [hours]
                    elif drift_ref_epoch_type == 'campaign':
                        epoch_to = row['delta_t_campaign_h']  # [hours]
                    ref_epoch_dt = previous_row['epoch_dt'] + (row['epoch_dt'] - previous_row['epoch_dt']) / 2  # mean

                    # Populate matrices and vectors:
                    mat_L0[(diff_obs_id, 0)] = g_diff_mugal
                    mat_sig_ll0[diff_obs_id] = sd_g_diff_mugal ** 2
                    # Partial derivative for g at stations:
                    if self.observed_stations.index(station_name_to) != self.observed_stations.index(station_name_from):
                        mat_A0[diff_obs_id, self.observed_stations.index(station_name_to)] = 1
                        mat_A0[diff_obs_id, self.observed_stations.index(station_name_from)] = -1
                    # Partial derivative for drift polynomial:
                    for pd_drift_id in range(drift_pol_degree):
                        mat_A0[diff_obs_id, pd_drift_col_offset + pd_drift_id + 1 + survey_count * drift_pol_degree] = \
                            epoch_to ** (pd_drift_id + 1) - epoch_from ** (pd_drift_id + 1)

                    # Log data in DataFrame:
                    g_diff_obs_mugal_list.append(g_diff_mugal)
                    sd_g_diff_obs_mugal_list.append(sd_g_diff_mugal)
                    station_name_from_list.append(station_name_from)
                    station_name_to_list.append(station_name_to)
                    diff_obs_id_list.append(diff_obs_id)
                    survey_names_list.append(survey_name)
                    setup_id_from_list.append(setup_id_from)
                    setup_id_to_list.append(setup_id_to)
                    ref_epoch_dt_list.append(ref_epoch_dt)

                    previous_row = row  # Store old row

        None_list_placeholder = [None] * len(survey_names_list)
        self.setup_obs_df = pd.DataFrame(list(zip(survey_names_list,
                                                  ref_epoch_dt_list,
                                                  diff_obs_id_list,
                                                  station_name_from_list,
                                                  station_name_to_list,
                                                  setup_id_from_list,
                                                  setup_id_to_list,
                                                  g_diff_obs_mugal_list,
                                                  sd_g_diff_obs_mugal_list,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  )),
                                         columns=self._SETUP_DIFF_COLUMNS)

        # - constraints:
        datum_station_id = -1
        for index, row in stat_df_obs_datum.iterrows():
            datum_station_id += 1
            station_name = row['station_name']
            station_id = self.observed_stations.index(station_name)
            mat_Ac[datum_station_id, station_id] = 1  # Partial derivative
            mat_Lc[(datum_station_id, 0)] = row['g_mugal']  # g for datum definition
            sd_mugal_for_weighting = row['sd_g_mugal'] / scaling_factor_datum_observations
            mat_sig_llc[datum_station_id] = sd_mugal_for_weighting ** 2

        # Set up all required matrices:
        mat_A = np.vstack((mat_A0, mat_Ac))  # Eq. (16)
        mat_L = np.vstack((mat_L0, mat_Lc))  # Eq. (16)
        mat_sig_ll = np.diag(np.hstack((mat_sig_ll0, mat_sig_llc)))
        mat_Qll = mat_sig_ll / (sig0_mugal ** 2)
        mat_P = np.linalg.inv(mat_Qll)

        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'---- Solve equation system: Results and statistics ----\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Solve equation system:
        mat_N = mat_A.T @ mat_P @ mat_A  # Normal equation matrix
        mat_Qxx = np.linalg.inv(mat_N)  # Co-factor matrix of estimates
        mat_x = mat_Qxx @ (mat_A.T @ mat_P @ mat_L)  # Estimates
        mat_v = (mat_A @ mat_x) - mat_L  # Post-fit residuals
        mat_Qldld = mat_A @ mat_Qxx @ mat_A.T  # A posteriori Co-factor matrix of adjusted observations
        # mat_Qll = np.linalg.inv(mat_P)
        mat_Qvv = mat_Qll - mat_Qldld  # Co-factor matrix of post-fit residuals

        # Test: "Gewichtsreziprokenprobe nach Ansermet" (see Skriptum AG1, p. 136, Eq. (6.86))
        u = np.sum(np.diag((mat_P @ mat_Qldld)))  # number of unknown parameters (estimates)
        tmp_diff = np.abs(number_of_parameters - u)
        if np.abs(number_of_parameters - u) > settings.ANSERMET_DIFF_TRESHOLD:
            raise AssertionError(f'"Gewichtsreziprokenprobe nach Ansermet" failed! Difference = {tmp_diff}')
        else:
            if verbose or self.write_log:
                tmp_str = f'# Gewichtsreziprokenprobe nach Ansermet (difference = {tmp_diff}) => Passed!\n'
                if verbose:
                    print(tmp_str)
                if self.write_log:
                    self.log_str += tmp_str

        # Condition of normal equation matrix:
        if verbose or self.write_log:
            cond = np.linalg.cond(mat_N)
            tmp_str = f'Condition of normal equation matix N = {cond:1.3f}\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # A posteriori variance of unit weight s02:
        dof = mat_A.shape[0] - mat_A.shape[1]  # degree of freedom
        par_r = mat_v.T @ mat_P @ mat_v  # = v^T * P * v
        if dof == 0:
            # s02_a_posteriori_mugal2 = par_r[0][0]
            raise AssertionError('Degree of freedom has to be larger than 0!')
        else:
            s02_a_posteriori_mugal2 = par_r[0][0] / dof  # Eq. (20)

        s0_mugal = np.sqrt(s02_a_posteriori_mugal2)  # A posteriori std. deviation of unit weight
        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'A posteriori variance (sd) of unit weight: ' \
                       f'{s02_a_posteriori_mugal2:5.3f} ({s0_mugal:5.3f})\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # ### Statistics and tests ###

        # Convert co-factor matrices to covariance matrices (variances in diagonal vector)
        mat_Cvv = s02_a_posteriori_mugal2 * mat_Qvv  # A posteriori Covariance matrix of post-fit residuals
        mat_Cxx = s02_a_posteriori_mugal2 * mat_Qxx  # A posteriori Covariance matrix of estimated paramaters
        mat_Cldld = s02_a_posteriori_mugal2 * mat_Qldld  # A posteriori Covariance matrix of adjusted observations
        # diag_Qxx = np.diag(mat_Qxx)

        # Calculate standard deviations:
        mat_sd_xx = np.sqrt(np.diag(mat_Cxx))  # A posteriori SD of estimates
        mat_sd_ldld = np.sqrt(np.diag(mat_Cldld))  # A posteriori SD of adjusted observations
        # mat_sd_vv = np.sqrt(np.diag(mat_Cvv))  # A posteriori SD of residuals
        mat_sd_vv = np.sqrt(
            np.diag(abs(mat_Cvv)))  # !!!! without "abs()" the sqrt operation fails because neg. values may occure!

        # creating histogram from residuals
        residual_hist, bin_edges = create_hist(mat_v)  # Calculate histogram

        # goodness-of-fit test
        chi_crit, chi_val, chi_test = goodness_of_fit_test(confidence_level_chi_test, dof,
                                                           s02_a_posteriori_mugal2, sig0_mugal ** 2)
        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'# Goodness-of-fit test results:\n'
            tmp_str += f'# Chi-val  Chi-crt-lower  Chi-crt-upper  Status\n'
            tmp_str += f'  {chi_val:7.3f}  {chi_crit[0]:13.3f}  {chi_crit[1]:13.3f}  {chi_test:s}\n'
            tmp_str += f'\n'
            tmp_str += f'# Histogram of the residuals:\n'
            tmp_str += f'# Lower-edge(??Gal)  Upper-edge(??Gal)  Frequency\n'
            for loop_1 in range(len(residual_hist)):
                tmp_str += f'  {bin_edges[loop_1]:16.4f}  {bin_edges[loop_1 + 1]:16.4f}  {residual_hist[loop_1]:9.0f}\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # outlier detection effectiveness (redundancy components)
        diag_Qvv = np.diag(mat_Qvv)  # TODO: Still needed?
        # mat_R = np.diag(mat_P) * diag_Qvv
        # mat_R is exactly the same as "mat_r"

        # Redundanzanteile (redundancy components):
        # - AG II, pp. 66-71
        mat_r = np.diag(mat_Qvv @ mat_P)

        # Standardisierte Versesserungen (AG II, p. 66)
        # - Normalverteilt mit Erwartwarungswert = 0 (wie Verbesserungen)
        # - Standardabweichung = 1 (standardisiert)
        mat_w = mat_v[:, 0] / mat_sd_vv

        # Tau test for outlier detection:
        alpha_tau = 1 - confidence_level_tau_test
        tau_test_result, tau_critical_value = tau_test(mat_w=mat_w, dof=dof, alpha=alpha_tau, mat_r=mat_r)
        number_of_outliers = tau_test_result.count("failed")

        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'# Tau-test results:\n'
            tmp_str += f'Critical value: {tau_critical_value:1.3f}\n'
            tmp_str += f' - Number of detected outliers: {number_of_outliers}\n'
            tmp_str += f' - Number low redundancy component: {tau_test_result.count("r too small")}\n'
            tmp_str += f'\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Blunder detection parameters calculation
        # sv_tau = 1 - confidence_level_tau_test
        # std_res, tau_val, tau_crt = tau_criterion_test(diag_Qvv, mat_r, mat_v, s02_a_posteriori_mugal2, dof, sv_tau)

        # #### Store results ####
        g_est_mugal = mat_x[0:number_of_stations, 0]
        sd_g_est_mugal = mat_sd_xx[0:number_of_stations]
        drift_pol_coeff = mat_x[number_of_stations:, 0]
        drift_pol_coeff_sd = mat_sd_xx[number_of_stations:]
        sd_diff_obs_mugal = mat_sd_ldld[:number_of_diff_obs]
        sd_pseudo_obs_mugal = mat_sd_ldld[number_of_diff_obs:]
        v_diff_obs_mugal = mat_v[:number_of_diff_obs, 0]
        v_pseudo_obs_mugal = mat_v[number_of_diff_obs:, 0]
        w_diff_mugal = mat_w[:number_of_diff_obs]
        w_pseudo_obs_mugal = mat_w[number_of_diff_obs:]
        r_diff_obs = mat_r[:number_of_diff_obs]
        r_pseudo_obs = mat_r[number_of_diff_obs:]
        # TODO: Add Tau criterion data here
        tau_test_result_diff_obs = tau_test_result[:number_of_diff_obs]
        tau_test_result_pseudo_obs = tau_test_result[number_of_diff_obs:]

        # Station related results:
        for idx, stat_name in enumerate(self.observed_stations):
            filter_tmp = self.stat_obs_df['station_name'] == stat_name
            self.stat_obs_df.loc[filter_tmp, 'g_est_mugal'] = g_est_mugal[idx]
            self.stat_obs_df.loc[filter_tmp, 'sd_g_est_mugal'] = sd_g_est_mugal[idx]
        # Calculate differences to estimates:
        self.stat_obs_df['diff_g_est_mugal'] = self.stat_obs_df['g_est_mugal'] - self.stat_obs_df['g_mugal']
        self.stat_obs_df['diff_sd_g_est_mugal'] = self.stat_obs_df['sd_g_est_mugal'] - self.stat_obs_df['sd_g_mugal']

        # Drift parameters:
        survey_name_list = []
        degree_list = []
        coefficient_list = []
        sd_coeff_list = []
        coeff_unit_list = []
        ref_epoch_t0_dt_list = []
        tmp_idx = 0
        x_estimate_drift_coeff_names = []
        for survey_name, setup_data in self.setups.items():
            for degree in range(drift_pol_degree):
                survey_name_list.append(survey_name)
                degree_list.append(degree + 1)  # starts with 1
                coefficient_list.append(drift_pol_coeff[tmp_idx])  # * (3600**(degree + 1))  # [??Gal/h]
                sd_coeff_list.append(drift_pol_coeff_sd[tmp_idx])
                coeff_unit_list.append(f'??Gal/h^{degree + 1}')
                if drift_ref_epoch_type == 'survey':
                    ref_epoch_t0_dt_list.append(setup_data['ref_epoch_delta_t_h'])
                elif drift_ref_epoch_type == 'campaign':
                    ref_epoch_t0_dt_list.append(setup_data['ref_epoch_delta_t_campaign_h'])
                else:
                    ref_epoch_t0_dt_list.append(None)  # Should not happen!
                tmp_idx += 1
                x_estimate_drift_coeff_names.append(f'{survey_name}-{degree + 1}')
        self.drift_pol_df = pd.DataFrame(list(zip(survey_name_list,
                                                  degree_list,
                                                  coefficient_list,
                                                  sd_coeff_list,
                                                  coeff_unit_list,
                                                  ref_epoch_t0_dt_list)),
                                         columns=self._DRIFT_POL_DF_COLUMNS)

        for idx, v_diff_mugal in enumerate(v_diff_obs_mugal):
            filter_tmp = self.setup_obs_df['diff_obs_id'] == idx
            self.setup_obs_df.loc[filter_tmp, 'v_diff_mugal'] = v_diff_mugal
            self.setup_obs_df.loc[filter_tmp, 'sd_g_diff_est_mugal'] = sd_diff_obs_mugal[idx]
            self.setup_obs_df.loc[filter_tmp, 'w_diff_mugal'] = w_diff_mugal[idx]  # standardized residuals
            self.setup_obs_df.loc[filter_tmp, 'r_diff_obs'] = r_diff_obs[
                idx]  # redundancy components
            self.setup_obs_df.loc[filter_tmp, 'tau_test_result'] = tau_test_result_diff_obs[idx]  # str

        # Print results to terminal:
        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f' - Station data:\n'
            tmp_str += self.stat_obs_df[['station_name', 'is_datum', 'g_mugal', 'g_est_mugal',
                                         'diff_g_est_mugal', 'sd_g_mugal',
                                         'sd_g_est_mugal']].to_string(index=False,
                                                                      float_format=lambda x: '{:.1f}'.format(x))
            tmp_str += f'\n\n'
            tmp_str += f' - Drift polynomial coefficients:\n'
            tmp_str += self.drift_pol_df.to_string(index=False, float_format=lambda x: '{:.6f}'.format(x))
            tmp_str += f'\n\n'
            tmp_str += f' - Differential observations:\n'
            for survey_name in survey_names:
                filter_tmp = self.setup_obs_df['survey_name'] == survey_name
                tmp_str += f'   - Survey: {survey_name}\n'
                tmp_str += self.setup_obs_df.loc[filter_tmp, ['station_name_from', 'station_name_to', 'g_diff_mugal',
                                                              'sd_g_diff_mugal', 'sd_g_diff_est_mugal',
                                                              'v_diff_mugal', 'r_diff_obs']].to_string(index=False,
                                                                                                       float_format=lambda
                                                                                                           x: '{:.2f}'.format(
                                                                                                           x))
                tmp_str += f'\n\n'
            tmp_str += f' - Pseudo observations at datum stations (constraints):\n'
            tmp_str += f'Station name  sd [??Gal]   v [??Gal]   w [??Gal]   r [0-1]    Tau test result\n'
            for idx, station_name in enumerate(datum_stations):
                tmp_str += f'{station_name:10}   {sd_pseudo_obs_mugal[idx]:8.3}     {v_pseudo_obs_mugal[idx]:+8.3}     {w_pseudo_obs_mugal[idx]:+5.3}     {r_pseudo_obs[idx]:+5.3}  {tau_test_result_pseudo_obs[idx]}\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Save data/infos to object for later use:
        self.drift_polynomial_degree = drift_pol_degree
        self.sig0_a_priori = sig0_mugal
        self.scaling_factor_datum_observations = scaling_factor_datum_observations
        self.confidence_level_chi_test = confidence_level_chi_test
        self.confidence_level_tau_test = confidence_level_chi_test
        self.number_of_stations = number_of_stations
        self.number_of_datum_stations = number_of_datum_stations
        self.number_of_estimates = number_of_parameters
        self.degree_of_freedom = dof
        self.s02_a_posteriori = s02_a_posteriori_mugal2
        self.Cxx = mat_Cxx
        self.x_estimate_names = self.observed_stations + x_estimate_drift_coeff_names
        self.goodness_of_fit_test_status = chi_test
        self.number_of_outliers = number_of_outliers
        self.drift_ref_epoch_type = drift_ref_epoch_type

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


if __name__ == '__main__':
    test2 = LSMDiff('a', 'b')  # This test will cause an error due to invalid data types!
    pass
