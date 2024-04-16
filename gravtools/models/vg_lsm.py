"""Classes for the estimation of vertical gravity gradient by least-squares adjustment of differential relative
gravimeter observations.

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
from gravtools.models.lsm import LSM, create_hist, global_model_test, tau_test
from gravtools.models import misc
from gravtools import __version__ as GRAVTOOLS_VERSION


class VGLSM(LSM):
    """VG estimation by least-squares adjustment of non-differential gravimeter observations.

    Notes
    -----
    The estimation of vertical gradients is restricted to th observations of a single survey and one station!


    Attributes
    ----------
    setup_obs_df : Pandas DataFrame
        Pandas Dataframes for logging (differential or absolute) setup observations, the related metadata, estimation
        results and statistics. The columns of the dataframe are defined in `self._SETUP_OBS_COLUMNS`.
    observed_stations : list of str
        Unique list of names of observed stations, defining the order of stations (station IDs) for all matrices and
        vectors used in the adjustment scheme.
    """

    # Column names of self.setup_obs_df:
    # - keys: Column names of the pandas dataframe
    # - values: Short description for table headers, etc., in the GUI
    _SETUP_OBS_COLUMNS_DICT = {
        'survey_name': 'Survey',
        'ref_epoch_dt': 'Epoch',
        'obs_id': 'Obs. ID',
        'station_name': 'Station',
        'setup_id': 'Setup ID',
        'g_obs_mugal': 'g [µGal]',
        'sd_g_obs_mugal': 'SD [µGal]',
        'dhf_sensor_m': 'Sensor height [m]',
        'sd_g_obs_est_mugal': 'SD_est [µGal]',
        'v_obs_est_mugal': 'Residuals [µGal]',  # Post fit residuals
        'sd_v_obs_est_mugal': 'SD_v [µGal]',  # SD of post-fit residuals
        'w_obs_est_mugal': 'Std. Residual []',
        'r_obs_est': 'Redundancy []',
        'tau_test_result': 'Outlier Test',
    }
    _SETUP_OBS_COLUMNS = list(_SETUP_OBS_COLUMNS_DICT.keys())

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

    # Column names of self.vg_pol_df:
    # - keys: Column names of the pandas dataframe
    # - values: Short description for table headers, etc., in the GUI
    _VG_POL_DF_COLUMNS_DICT = {
        'degree': 'Degree',
        'coefficient': 'Coefficient',
        'sd_coeff': 'SD',
        'coeff_unit': 'Unit',
    }
    _VG_POL_DF_COLUMNS = list(_VG_POL_DF_COLUMNS_DICT.keys())

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
        lsm_method = 'VG_LSM_nondiff'
        super().__init__(lsm_method, stat_df=stat_df, setups=setups, comment=comment, write_log=write_log)

    @classmethod
    def from_campaign(cls, campaign, comment='', write_log=True):
        """Constructor that generates and populates the LSM object from a Campaign class object.

        Notes
        -----
        Put checks that are dependent on the LSM method into here! All common checks and preparations are implemented
        in the `from_campaign` method of the parent class `LSM`.

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
        :py:obj:`.VGLSM`
            Contains all information required for adjusting the campaign.
        """
        # Station data:
        # - Just one station (at different height levels) allowed => Check coordinates!
        if campaign.stations is None:
            raise AssertionError(f'The campaign "{campaign.campaign_name}" does not contain any station data!')
        else:
            if not hasattr(campaign.stations, 'stat_df'):
                raise AssertionError(f'The campaign "{campaign.campaign_name}" has not station dataframe!')
            else:
                if len(campaign.stations.stat_df) == 0:
                    raise AssertionError(f'The campaign "{campaign.campaign_name}" has an empty station dataframe!')
                else:  # Just one station allowed => Check the coordinates (lon, lat)!
                    stat_df_observed = campaign.stations.stat_df.loc[campaign.stations.stat_df.is_observed]
                    if (len(stat_df_observed['lat_deg'].unique()) > 1) or (
                            len(stat_df_observed['long_deg'].unique()) > 1):
                        raise AssertionError(f'This campaign contains observed stations with differing latitudes and/or'
                                             f' longitudes. This is an indicator for observations of more than one '
                                             f'stations which is not allowed.')

        # Survey data:
        # - Just one survey allowed!
        if campaign.surveys is None:
            raise AssertionError(f'The campaign "{campaign.campaign_name}" does not contain any survey data!')
        else:
            if len(campaign.surveys) > 1:
                raise AssertionError(
                    f'The campaign "{campaign.campaign_name}" contains more than one survey! The VG estimation is '
                    f'restricted to just one survey.')

        # Check if the reduced setup observations refer to the sensor height of the instrument:
        # Loop over surveys in campaign (just one...):
        for survey_name, survey in campaign.surveys.items():
            if survey.keep_survey and survey.is_active:
                if survey.red_reference_height_type != 'sensor_height':
                    raise AssertionError(f'In survey {survey_name} the reduced gravity values do not refer to the '
                                         f'sensor height! Change the reference height to "Sensor" and recalculate the '
                                         f'setup observation data.')

        return super().from_campaign(campaign, comment=comment, write_log=True)

    def adjust(self, drift_pol_degree=1,
               vg_polynomial_degree=1,
               vg_polynomial_ref_height_offset_m=0.0,
               sig0_mugal=1,
               confidence_level_chi_test=0.95,
               confidence_level_tau_test=0.95,
               verbose=False
               ):
        """Run the VG adjustment based on non-differential observations.

        Parameters
        ----------
        drift_pol_degree : int, optional (default=1)
            Degree of estimated drift polynomial.
        vg_polynomial_degree : int, optional (default=1)
            Degree of the estimated polynomial of the vertical gravity gradient. Valid values are 1, 2 or 3.
        vg_polynomial_ref_height_offset_m : float, optional (default=0.0)
            Vertical offset [m] between the reference height of the control point and the zero-level of the estimated
            VG polynomial. A positive offset describes an offset above the control point.
        sig0_mugal : int, optional (default=1)
            A priori standard deviation of unit weight of observations [µGal] for the stochastic model of the
            least-squares adjustment.
        confidence_level_chi_test : float, optional (default=0.95)
            Confidence level for the goodness-of-fit test.
        confidence_level_tau_test : float, optional (default=0.95)
            Confidence level for the tau test.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line, e.g. for debugging and testing

        Notes
        -----
        """

        # Initial checks:
        if (vg_polynomial_degree < 1) or (vg_polynomial_degree > 3):
            raise AssertionError(f'Invalid degree of the VG polynomial (degree = {vg_polynomial_degree})! Valid values '
                                 f'are 1 to 3.')

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
        number_of_surveys = len(self.setups)  # Has to be 1 anyway

        # ### Station data ###
        # Get dataframe with subset of observed stations only:
        filter_tmp = self.stat_df['station_name'].isin(self.observed_stations)
        self.stat_obs_df = self.stat_df.loc[filter_tmp].copy(deep=True)  # All observed stations
        # Drop columns that are not required at VG estimation:
        self.stat_obs_df.drop(columns=['g_mugal', 'sd_g_mugal', 'is_datum'], inplace=True)
        # Calculate statistical data for the sensor height (relevant information for VG estimation):
        tmp_df = setup_df.loc[:, ['station_name', 'dhf_sensor_m']].groupby('station_name').mean().rename(
            columns={"dhf_sensor_m": "dhf_sensor_mean_m"})
        self.stat_obs_df = self.stat_obs_df.merge(tmp_df, left_on='station_name', right_on='station_name', how='left')
        tmp_df = setup_df.loc[:, ['station_name', 'dhf_sensor_m']].groupby('station_name').std().rename(
            columns={"dhf_sensor_m": "dhf_sensor_std_m"})
        self.stat_obs_df = self.stat_obs_df.merge(tmp_df, left_on='station_name', right_on='station_name', how='left')
        tmp_df = setup_df.loc[:, ['station_name', 'dhf_sensor_m']].groupby('station_name').min().rename(
            columns={"dhf_sensor_m": "dhf_sensor_min_m"})
        self.stat_obs_df = self.stat_obs_df.merge(tmp_df, left_on='station_name', right_on='station_name', how='left')
        tmp_df = setup_df.loc[:, ['station_name', 'dhf_sensor_m']].groupby('station_name').max().rename(
            columns={"dhf_sensor_m": "dhf_sensor_max_m"})
        self.stat_obs_df = self.stat_obs_df.merge(tmp_df, left_on='station_name', right_on='station_name', how='left')

        if number_of_surveys > 1:
            raise AssertionError(
                f'Invalid number of surveys ({number_of_surveys})! Only one survey allowed for VG estimation.')
        # Total number of parameters to be estimated:
        # - Drift polynomial coeff.: Polynomial degree * number of surveys
        # - 1 constant instrumental bias per survey
        # - VG polynomial coeff.: Polynomial degree * number of surveys
        number_of_parameters = drift_pol_degree + vg_polynomial_degree + 1

        # Check, if setup IDs are unique:
        self.check_unique_setups(setup_ids)

        if verbose or self.write_log:
            time_now_str = dt.datetime.now(tz=pytz.UTC).strftime('%Y-%m-%d, %H:%M:%S %Z')
            tmp_str = f'#### Adjustment log (VG LGM estimation based on non-differential observations) ####\n'
            tmp_str += f'Processed with GravTools {GRAVTOOLS_VERSION} ({time_now_str})\n'
            tmp_str += f'Comment: {self.comment}\n'
            tmp_str += f'\n'
            tmp_str += f'---- Input data and settings ----\n'
            tmp_str += f'Method: {settings.ADJUSTMENT_METHODS[self.lsm_method]}\n'
            tmp_str += f'Number of stations (i.e. height levels): {number_of_stations}\n'
            tmp_str += f'Number of observations: {number_of_observations}\n'
            tmp_str += f'Number of estimated parameters: {number_of_parameters}\n'
            tmp_str += f'Degree of freedom: {number_of_observations - number_of_parameters}\n'
            tmp_str += f'\n'
            tmp_str += f'Degree of drift polynomial: {drift_pol_degree}\n'
            tmp_str += f'Degree of VG polynomial: {vg_polynomial_degree}\n'
            tmp_str += f'VG polynomial height offset [m]: {vg_polynomial_ref_height_offset_m:4.3f}\n'
            tmp_str += f'A priori std. deviation of unit weight [µGal]: {sig0_mugal}\n'
            tmp_str += f'Confidence level Chi-test: {confidence_level_chi_test:4.2f}\n'
            tmp_str += f'Confidence level Tau-test: {confidence_level_tau_test:4.2f}\n'
            tmp_str += f'\n'
            tmp_str += f'---- Survey infos ----\n'
            tmp_str += self.survey_info_string
            tmp_str += f'\n'
            tmp_str += f'\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Initialize matrices:
        if verbose or self.write_log:
            tmp_str = f'---- Set up and populate matrices ----\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str
        # => Initialize complete matrices first and then populate them. This is most efficient!
        # - Observation model:
        mat_A = np.zeros([number_of_observations, number_of_parameters])  # Design matrix
        mat_L = np.zeros((number_of_observations, 1))  # Observations
        mat_sig_ll = np.zeros(number_of_observations)  # Variances of obs.

        # Populate matrices:
        obs_id = -1  # Index of differential observations in vectors mat_L, rows of mat_A and mat_p0
        survey_count = -1  # Survey counter for indexing the drift parameters in the A-matrix
        g_obs_mugal_list = []
        station_name_list = []
        setup_id_list = []
        sd_g_obs_mugal_list = []
        obs_id_list = []
        survey_names_list = []
        ref_epoch_dt_list = []  # Reference epochs (datetime objects) of observations
        dhf_sensor_m_list = []  # Height above the control point for VG estimation

        for survey_name, setup_data in self.setups.items():
            setup_df = setup_data['setup_df']
            survey_count += 1

            # Add offset to the height of the sensor above the control point:
            # setup_df['dhf_sensor_m'] = setup_df['dhf_sensor_m'] + vg_polynomial_ref_height_offset_m
            if (setup_df['sd_g_mugal'] < 0).any():
                raise AssertionError(
                    f'ERROR: SD of observations ("sd_g_mugal") in survey {survey_name} <= 0 are not allowed!')

            for index, row in setup_df.iterrows():
                obs_id += 1  # Increment diff. observation ID
                g_obs_mugal = row['g_mugal']
                sd_g_obs_mugal = row['sd_g_mugal']
                station_name = row['station_name']
                setup_id = row['setup_id']
                delta_t_h = row['delta_t_h']  # [hours]
                ref_epoch_dt = row['epoch_dt']
                dhf_sensor_m = row['dhf_sensor_m'] + vg_polynomial_ref_height_offset_m

                # Populate matrices and vectors:
                mat_L[(obs_id, 0)] = g_obs_mugal
                mat_sig_ll[obs_id] = sd_g_obs_mugal ** 2
                # Partial derivative for drift polynomial including constant instrumental bias (pol. degree = 0):
                for pd_drift_id in range(drift_pol_degree + 1):
                    mat_A[obs_id, pd_drift_id] = \
                        delta_t_h ** pd_drift_id
                # Partial derivative for VG polynomial:
                for pd_vg_id in range(vg_polynomial_degree):
                    mat_A[obs_id, drift_pol_degree + 1 + pd_vg_id] = \
                        dhf_sensor_m ** (pd_vg_id + 1)

                # Log data in DataFrame:
                g_obs_mugal_list.append(g_obs_mugal)
                sd_g_obs_mugal_list.append(sd_g_obs_mugal)
                station_name_list.append(station_name)
                obs_id_list.append(obs_id)
                survey_names_list.append(survey_name)
                setup_id_list.append(setup_id)
                ref_epoch_dt_list.append(ref_epoch_dt)
                dhf_sensor_m_list.append(dhf_sensor_m)

        None_list_placeholder = [None] * len(survey_names_list)
        self.setup_obs_df = pd.DataFrame(list(zip(survey_names_list,
                                                  ref_epoch_dt_list,
                                                  obs_id_list,
                                                  station_name_list,
                                                  setup_id_list,
                                                  g_obs_mugal_list,
                                                  sd_g_obs_mugal_list,
                                                  dhf_sensor_m_list,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  None_list_placeholder,
                                                  )),
                                         columns=self._SETUP_OBS_COLUMNS)

        # Set up all required matrices:
        mat_sig_ll = np.diag(mat_sig_ll)
        mat_Qll = mat_sig_ll / (sig0_mugal ** 2)
        mat_P = np.linalg.inv(mat_Qll)

        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'---- Results and statistics ----\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Solve equation system:
        mat_N = mat_A.T @ mat_P @ mat_A  # Normal equation matrix
        mat_Qxx = np.linalg.inv(mat_N)  # Co-factor matrix of estimates
        mat_x = mat_Qxx @ (mat_A.T @ mat_P @ mat_L)  # Estimates
        mat_v = (mat_A @ mat_x) - mat_L  # Post-fit residuals
        mat_v = misc.numpy_array_set_zero(mat_v, atol=1e-4)  # ToDo: Check atol!
        mat_Qldld = mat_A @ mat_Qxx @ mat_A.T  # A posteriori Co-factor matrix of adjusted observations
        mat_Qldld = misc.numpy_array_set_zero(mat_Qldld)
        # mat_Qll = np.linalg.inv(mat_P)
        mat_Qvv = mat_Qll - mat_Qldld  # Co-factor matrix of post-fit residuals
        mat_Qvv = misc.numpy_array_set_zero(mat_Qvv)

        # Test: "Gewichtsreziprokenprobe nach Ansermet" (see Skriptum AG1, p. 136, Eq. (6.86))
        u = np.sum(np.diag((mat_P @ mat_Qldld)))  # number of unknown parameters (estimates)
        tmp_diff = np.abs(number_of_parameters - u)
        if np.abs(number_of_parameters - u) > settings.ANSERMET_DIFF_THRESHOLD:
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
        # mat_sd_vv = np.sqrt(np.diag(abs(mat_Cvv)))  # !!!! without "abs()" the sqrt operation fails because neg. values may occur!

        # A posteriori SD of residuals:
        # - Check just in case, whether all diagonal elements of the Qvv matrix are positive!
        if (np.diag(mat_Qvv) < 0).any():
            mat_sd_vv = np.sqrt(np.diag(abs(mat_Cvv)))
            if verbose or self.write_log:
                tmp_str = f' - Warning: At least one diagonal element of the Qvv matrix is negative!\n'
                if verbose:
                    print(tmp_str)
                if self.write_log:
                    self.log_str += tmp_str
        else:
            mat_sd_vv = np.sqrt(np.diag(mat_Cvv))

        # creating histogram from residuals
        residual_hist, bin_edges = create_hist(mat_v)  # Calculate histogram

        # goodness-of-fit test
        chi_crit, chi_val, chi_test = global_model_test(confidence_level_chi_test, dof,
                                                        s02_a_posteriori_mugal2, sig0_mugal ** 2)
        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'# Goodness-of-fit test results:\n'
            tmp_str += f'# Chi-val  Chi-crt-lower  Chi-crt-upper  Status\n'
            tmp_str += f'  {chi_val:7.3f}  {chi_crit[0]:13.3f}  {chi_crit[1]:13.3f}  {chi_test:s}\n'
            tmp_str += f'\n'
            tmp_str += f'# Histogram of the residuals:\n'
            tmp_str += f'# Lower-edge(µGal)  Upper-edge(µGal)  Frequency\n'
            for loop_1 in range(len(residual_hist)):
                tmp_str += f'  {bin_edges[loop_1]:16.4f}  {bin_edges[loop_1 + 1]:16.4f}  {residual_hist[loop_1]:9.0f}\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Redundanzanteile (redundancy components):
        # - Measure for the outlier detection effectiveness
        # - AG II, pp. 66-71
        mat_r = np.diag(mat_Qvv @ mat_P)

        # Standardized residuals used as test statistics for outlier detection:
        # - AG II, p. 66
        # - Taking care of zero-elements in the mat_sd_vv vector to prevent division by zero errors!
        mat_w = np.zeros(len(mat_v))
        tmp_filter = ~np.isclose(mat_v[:, 0], 0.0)
        mat_w[tmp_filter] = mat_v[tmp_filter, 0] / mat_sd_vv[tmp_filter]

        # Tau test for outlier detection:
        alpha_tau = 1 - confidence_level_tau_test
        tau_test_result, tau_critical_value = tau_test(mat_w=mat_w, dof=dof, alpha=alpha_tau, mat_r=mat_r)
        number_of_outliers = tau_test_result.count("failed")

        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'# Tau-test results:\n'
            tmp_str += f'Critical value (for testing w): {tau_critical_value:1.3f}\n'
            tmp_str += f' - Number of detected outliers: {number_of_outliers}\n'
            tmp_str += f' - Number low redundancy component: {tau_test_result.count("r too small")}\n'
            tmp_str += f'\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # #### Store results ####
        drift_pol_coeff = mat_x[:drift_pol_degree + 1, 0]
        drift_pol_coeff_sd = mat_sd_xx[:drift_pol_degree + 1]
        vg_pol_coeff = mat_x[drift_pol_degree + 1:, 0]
        vg_pol_coeff_sd = mat_sd_xx[drift_pol_degree + 1:]
        sd_g_obs_est_mugal = mat_sd_ldld
        v_obs_est_mugal = mat_v
        sd_v_obs_est_mugal = mat_sd_vv
        w_obs_est_mugal = mat_w
        r_obs_est = mat_r
        tau_test_result_obs = tau_test_result

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
            for degree in range(drift_pol_degree + 1):
                survey_name_list.append(survey_name)
                degree_list.append(degree)  # starts with 1
                coefficient_list.append(drift_pol_coeff[tmp_idx])  # * (3600**(degree + 1))  # [µGal/h]
                sd_coeff_list.append(drift_pol_coeff_sd[tmp_idx])
                ref_epoch_t0_dt_list.append(setup_data['ref_epoch_delta_t_h'])
                if degree == 0:
                    coeff_unit_list.append(f'µGal')
                else:
                    coeff_unit_list.append(f'µGal/h^{degree}')
                tmp_idx += 1
                x_estimate_drift_coeff_names.append(f'drift-{degree}')
        self.drift_pol_df = pd.DataFrame(list(zip(survey_name_list,
                                                  degree_list,
                                                  coefficient_list,
                                                  sd_coeff_list,
                                                  coeff_unit_list,
                                                  ref_epoch_t0_dt_list)),
                                         columns=self._DRIFT_POL_DF_COLUMNS)

        # VG parameters:
        degree_list = []
        coefficient_list = []
        sd_coeff_list = []
        coeff_unit_list = []
        tmp_idx = 0
        x_estimate_vg_coeff_names = []
        for degree in range(1, vg_polynomial_degree + 1):
            degree_list.append(degree)  # starts with 1
            coefficient_list.append(vg_pol_coeff[tmp_idx])  # * (3600**(degree + 1))  # [µGal/h]
            sd_coeff_list.append(vg_pol_coeff_sd[tmp_idx])
            if degree == 0:
                coeff_unit_list.append(f'µGal')
            else:
                coeff_unit_list.append(f'µGal/m^{degree}')
            tmp_idx += 1
            x_estimate_vg_coeff_names.append(f'vg-{degree}')
        self.vg_pol_df = pd.DataFrame(list(zip(degree_list,
                                               coefficient_list,
                                               sd_coeff_list,
                                               coeff_unit_list)),
                                      columns=self._VG_POL_DF_COLUMNS)

        # Observation-related results:
        for idx, v_obs_mugal in enumerate(v_obs_est_mugal):
            filter_tmp = self.setup_obs_df['obs_id'] == idx
            self.setup_obs_df.loc[filter_tmp, 'v_obs_est_mugal'] = v_obs_mugal
            self.setup_obs_df.loc[filter_tmp, 'sd_v_obs_est_mugal'] = sd_v_obs_est_mugal[idx]
            self.setup_obs_df.loc[filter_tmp, 'sd_g_obs_est_mugal'] = sd_g_obs_est_mugal[idx]
            self.setup_obs_df.loc[filter_tmp, 'w_obs_est_mugal'] = w_obs_est_mugal[idx]  # standardized residuals
            self.setup_obs_df.loc[filter_tmp, 'r_obs_est'] = r_obs_est[idx]  # redundancy components
            self.setup_obs_df.loc[filter_tmp, 'tau_test_result'] = tau_test_result_obs[idx]  # str

        # Print results to terminal and to log string:
        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f' - Station data:\n'
            tmp_str += self.stat_obs_df[['station_name', 'dhf_sensor_mean_m', 'dhf_sensor_std_m',
                                         'dhf_sensor_min_m', 'dhf_sensor_max_m']].to_string(index=False,
                                                                      float_format=lambda x: '{:.4f}'.format(x))
            tmp_str += f'\n\n'
            tmp_str += f' - Drift polynomial coefficients:\n'
            tmp_str += self.drift_pol_df.to_string(index=False, float_format=lambda x: '{:.6f}'.format(x))
            tmp_str += f'\n\n'
            tmp_str += f' - VG polynomial coefficients:\n'
            tmp_str += self.vg_pol_df.to_string(index=False, float_format=lambda x: '{:.6f}'.format(x))
            tmp_str += f'\n\n'
            tmp_str += f' - Observations:\n'
            for survey_name in survey_names:
                filter_tmp = self.setup_obs_df['survey_name'] == survey_name
                tmp_str += f'   - Survey: {survey_name}\n'
            tmp_str += self.setup_obs_df.loc[filter_tmp, ['station_name', 'g_obs_mugal',
                                                          'sd_g_obs_mugal', 'sd_g_obs_est_mugal',
                                                          'v_obs_est_mugal']].to_string(index=False,
                                                                                        float_format=lambda
                                                                                            x: '{:.1f}'.format(x))
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Save data/infos to object for later use:
        self.drift_polynomial_degree = drift_pol_degree
        self.sig0_a_priori = sig0_mugal
        self.confidence_level_chi_test = confidence_level_chi_test
        self.confidence_level_tau_test = confidence_level_chi_test
        self.number_of_stations = number_of_stations
        self.number_of_estimates = number_of_parameters
        self.degree_of_freedom = dof
        self.s02_a_posteriori = s02_a_posteriori_mugal2
        self.Cxx = mat_Cxx
        self.x_estimate_names = x_estimate_drift_coeff_names + x_estimate_vg_coeff_names
        self.global_model_test_status = chi_test
        self.number_of_outliers = number_of_outliers
        self.drift_ref_epoch_type = 'survey'  # Not relevant anyway, because only ONE survey allowed in the campaign!
        self.vg_polynomial_ref_height_offset_m = vg_polynomial_ref_height_offset_m
        self.vg_polynomial_degree = vg_polynomial_degree
        self.mat_A = mat_A
        self.mat_x = mat_x
        # The following attributes are None as initialized:
        # self.scaling_factor_datum_observations
        # self.number_of_datum_stations
