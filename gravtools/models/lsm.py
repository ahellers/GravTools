"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
Contains classes for least-squares adjustment of gravimeter-surveys.
"""

import pandas as pd
import numpy as np
from scipy import stats
import datetime as dt
import os
import copy

from gravtools.settings import SURVEY_DATA_SOURCE_TYPES, STATION_DATA_SOURCE_TYPES, GRAVIMETER_ID_BEV, \
    TIDE_CORRECTION_TYPES, DEFAULT_GRAVIMETER_ID_CG5_SURVEY, REFERENCE_HEIGHT_TYPE, NAME_OBS_FILE_BEV, \
    PATH_OBS_FILE_BEV, BEV_GRAVIMETER_TIDE_CORR_LOOKUP, GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m
from gravtools.const import VG_DEFAULT
from gravtools.models.exceptions import FileTypeError
from gravtools.CG5_utils.cg5_survey import CG5Survey


class LSM_diff:
    """Least-squares adjustment of differential gravimeter observations with weighted constraints.

    Attributes
    ----------
    stat_df : :py:obj:`gravtools.Station.stat_df`, optional (default=None)
            The station dataframe contains all relevant station data.
    setups : dict, optional (default=None)
        The setups dictionary contains all observation data used for the adjustment. The keys of the dictionary
        are the survey names (str) and the items are pandas dataframes containing the observation data (see
        :py:obj:`gravtool.Survey.setup_df`)
    setups_diff_df : Pandas DataFrame
        Pandas Dataframes for logging differential observations, the related metadata, estimation results and
        statistics. The columns of the dataframe are defined in `self._SETUP_DIFF_COLUMNS`.
    observed_stations : list of str
        Unique list of names of observed stations, defining the order of stations (station IDs) for all matrices and
        vectors used in the adjustment scheme.
    """

    _SETUP_DIFF_COLUMNS = (
        'survey_name',
        'diff_obs_id',
        'station_name_from',
        'station_name_to',
        'setup_id_from',
        'setup_id_to',
        'g_diff_mugal',
        'sd_g_diff_mugal',
        'sd_g_diff_est_mugal',
        'v_diff_mugal',
    )

    _DRIFT_POL_DF_COLUMNS = (
        'survey_name',
        'degree',
        'coefficient',
        'sd_coeff',
    )

    def __init__(self, stat_df=None, setups=None):
        """
        Parameters
        ----------
        stat_df : :py:obj:`gravtools.Station.stat_df`, optional (default=None)
            The station dataframe contains all relevant station data.
        setups : dict, optional (default=None)
            The setups dictionary contains all observation data used for the adjustment. The keys of the dictionary
            are the survey names (str) and the items are pandas dataframes containing the observation data (see
            :py:obj:`gravtool.Survey.setup_df`)

        Notes
        -----
        In order to prevent altering the original input data (in case of assignment via reference), deep copies of the
        input data are stored in objects of this class.

        """
        # Create deep copies of the input data items:
        self.stat_df = stat_df.copy(deep=True)
        self.setups = copy.deepcopy(setups)
        self.setups_diff_df = None  # Dataframe that containe the observation (setup) related results
        self.observed_stations = None  # Unique list of observed stations; defines the station IDs for matrices
        self.stat_obs_df = None  # Station dataframe that contains estimation results of observed stations
        self.drift_pol_df = None  # DataFrame that contains the estimated parameters of the drift polynomials an
        # their statistics for each survey

    @classmethod
    def from_campaign(cls, campaign):
        """Constructor that generates and populates the LSM object from a Campaign class object.

        Notes
        -----
        From all active surveys in the campaing the setup data (= observations) are loaded, additionally to the station
        data (station dataframe).

        Parameters
        ----------
        campaign : :py:obj:`gravtools.models.survey.Campaign`
            The campaign object needs to provide setup data for all active surveys and the related station data.

        Returns
        -------
        :py:obj:`.LSM_diff`
            Contains all information required for adjusting the campaign.
        """
        # Check if all required data is available in the campaign object:
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
                        setups[survey_name] = survey.setup_df

        # Check if setup data is available:
        if len(setups) == 0:
            raise AssertionError(f'Setup data of campaign "{campaign.campaign_name}" does not contain observations!')

        # Initialize and return LSM object:
        return cls(campaign.stations.stat_df, setups)

    def adjust(self, drift_pol_degree=1,
               sig0_mugal=1,
               scaling_factor_datum_observations=1.0,
               confidence_level_chi_test=0.95,
               confidence_level_tau_test=0.95,
               verbose=False
               ):  # Confidence level):
        """Run the adjustment.

        Parameters
        ----------
        drift_pol_degree : int, optional (default=1)
            Degree of estimated drift polynomial.
        sig0_mugal : int, optional (default=1)
            A priori standard deviation of unit weight of observations [µGal] for the stochastic model of the
            least-squares adjustment.
        scaling_factor_datum_observations : float, optional (default=1.0)
            Factor for scaling the standard deviation (SD) of g of datum stations. The scaled SD is is used for
            weighting the direct pseudo observations of g at the datum stations that are introduced as datum
            constraints.
        confidence_level_chi_test : float, optional (default=0.95)
            Confidence level for the goodness-of-fit test. 
        confidence_level_tau_test : float, optional (default=0.95)
            Confidence level for the tau test.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line, e.g. for debugging and testing

        Notes
        -----
        The adjustemnt method (weighted constraints adjustment based on differential observations) is documented in
        Wijaya et al. (2019) - pyGABEUR-ITB: A free Software for Adjustment of Relative Gravimeter Data.
        """

        # Prepare lists and indices:
        # - Observations and parameters:
        self.observed_stations = []
        survey_names = []
        number_of_observations = 0
        setup_ids = []
        for survey_name, setup_df in self.setups.items():
            self.observed_stations = self.observed_stations + setup_df['station_name'].to_list()
            number_of_observations = number_of_observations + len(setup_df)
            survey_names.append(survey_name)
            setup_ids = setup_ids + setup_df['setup_id'].to_list()
        self.observed_stations = list(
            set(self.observed_stations))  # Unique list of stations => order of stations in matrices!
        number_of_stations = len(self.observed_stations)
        number_of_surveys = len(self.setups)
        number_of_parameters = number_of_stations + drift_pol_degree * number_of_surveys  # Total number of parameters to be estimated
        number_of_diff_obs = number_of_observations - number_of_surveys

        if verbose:
            print(f'')
            print(f'#### Adjustment ####')
            print(f'---- Collect Information ----')
            print(f'Number of surveys: {number_of_surveys}')
            print(f'Number of stations: {number_of_stations}')
            print(f'Number of differential observations: {number_of_diff_obs}')
            print(f'Number of estimated parameters: {number_of_parameters}')

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

        if verbose:
            print(f'Number of datum stations: {number_of_datum_stations}')
            print(f'-------------------------------')
            print(f'')
            print(f'---- Set up and populate matrices ----')

        # Initialize matrices:
        # => Initialize complete matrices first and then populate them. This is most efficient!
        # - Observation model:
        mat_A0 = np.zeros([number_of_diff_obs, number_of_parameters])  # Model-matrix
        mat_L0 = np.zeros((number_of_diff_obs, 1))
        mat_p0 = np.zeros(number_of_diff_obs)  # Vector of weights for differential observations.
        # - Constraints:
        mat_Ac = np.zeros([number_of_datum_stations, number_of_parameters])  # Model-matrix for constraints
        mat_Lc = np.zeros((number_of_datum_stations, 1))
        mat_pc = np.zeros(number_of_datum_stations)
        #  => Convert to diagonal matrix: P0 = np.diag(np.array([1,2,3,4])) = np.diag(mat_p0)

        # Populate matrices:
        diff_obs_id = -1  # Index of differential observations in vectors mat_L0, rows of mat_A0 and mat_p0
        pd_drift_col_offset = number_of_stations - 1  # Column offset for drift parameters in A-matrix
        g_diff_obs_mugal_list = []
        station_name_from_list = []
        station_name_to_list = []
        setup_id_from_list = []
        setup_id_to_list = []
        sd_g_diff_obs_mugal_list = []
        diff_obs_id_list = []
        survey_names_list = []

        for survey_name, setup_df in self.setups.items():
            previous_row = None
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
                    epoch_from = previous_row['epoch_unix']  # [sec]
                    epoch_to = row['epoch_unix']  # [sec]

                    # TODO: Evtl. die startecpoche t0 anders setzen! Auf nummerische Stabilität achten!
                    # Jetzt in UNIX time [sec]

                    # Populate matrices and vectors:
                    mat_L0[(diff_obs_id, 0)] = g_diff_mugal
                    mat_p0[diff_obs_id] = sig0_mugal ** 2 / sd_g_diff_mugal ** 2
                    # Partial derivative for g at stations:
                    mat_A0[diff_obs_id, self.observed_stations.index(station_name_to)] = 1
                    mat_A0[diff_obs_id, self.observed_stations.index(station_name_from)] = -1
                    # Partial derivative for drift polynomial:
                    for pd_drift_id in range(drift_pol_degree):
                        mat_A0[diff_obs_id, pd_drift_col_offset + pd_drift_id + 1] = \
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

                    previous_row = row  # Store old row

        None_list_placeholder = [None] * len(survey_names_list)
        self.setups_diff_df = pd.DataFrame(list(zip(survey_names_list,
                                                    diff_obs_id_list,
                                                    station_name_from_list,
                                                    station_name_to_list,
                                                    setup_id_from_list,
                                                    setup_id_to_list,
                                                    g_diff_obs_mugal_list,
                                                    sd_g_diff_obs_mugal_list,
                                                    None_list_placeholder,
                                                    None_list_placeholder,
                                                    )),
                                           columns=self._SETUP_DIFF_COLUMNS)

        # - constraints:
        datum_station_id = -1
        for index, row in stat_df_obs_datum.iterrows():
            # print(f' - Station {row["station_name"]:10s} (row index {index:d})')
            datum_station_id += 1
            station_name = row['station_name']
            station_id = self.observed_stations.index(station_name)
            mat_Ac[datum_station_id, station_id] = 1  # Partial derivative
            mat_Lc[(datum_station_id, 0)] = row['g_mugal']  # g for datum definition
            sd_mugal_for_weighting = row['sd_g_mugal'] * scaling_factor_datum_observations
            mat_pc[datum_station_id] = sig0_mugal ** 2 / sd_mugal_for_weighting ** 2  # Weighting

        # Set up all required matrices:
        mat_A = np.vstack((mat_A0, mat_Ac))  # Eq. (16)
        mat_L = np.vstack((mat_L0, mat_Lc))  # Eq. (16)
        mat_P = np.diag(np.hstack((mat_p0, mat_pc)))  # Eq. (18)

        if verbose:
            print(f'')
            print(f'---- Solve equation system: Results and statistics ----')

        # Solve equation system:
        mat_N = mat_A.T @ mat_P @ mat_A
        mat_Qxx = np.linalg.inv(mat_N)  # Eq. (19)
        mat_x = mat_Qxx @ (mat_A.T @ mat_P @ mat_L)
        mat_v = (mat_A @ mat_x) - mat_L  # a posteriori residuals
        mat_Qldld = mat_A @ mat_Qxx @ mat_A.T
        mat_Qll = np.linalg.inv(mat_P)
        mat_Qvv = mat_Qll - mat_Qldld

        # Test: "Gewichtsreziprokenprobe nach Ansermet" (see Skriptum AG1, p. 136, Eq. (6.86))
        u = np.sum(np.diag((mat_P @ mat_Qldld)))
        tmp_diff = np.abs(number_of_parameters - u)
        if np.abs(number_of_parameters - u) > 1e-6:
            raise AssertionError(f'"Gewichtsreziprokenprobe nach Ansermet" failed! Difference = {tmp_diff}')
        else:
            if verbose:
                print(f'Gewichtsreziprokenprobe nach Ansermet: Difference = {tmp_diff} => Passed!')

        # A posteriori variance of unit weight s02:
        dof = mat_A.shape[0] - mat_A.shape[1]  # degree of freedom
        par_r = mat_v.T @ mat_P @ mat_v  # = v^T * P * v
        if dof == 0:
            s02_a_posteriori_mugal2 = par_r[0][0]  # ref_var
        else:
            s02_a_posteriori_mugal2 = par_r[0][0] / dof  # Eq. (20)
        if verbose:
            print(f'Degree of freedom: {dof}')
            print(f'A posteriori variance (sd) of unit weight '
                  f': {s02_a_posteriori_mugal2} ({np.sqrt(s02_a_posteriori_mugal2)})')

        # ### Statistics and tests ###

        # Convert co-factor matrices to covariance matrices (variances in diagonal vector)
        mat_Cvv = s02_a_posteriori_mugal2 * mat_Qvv  # Residuals of observations
        mat_Cxx = s02_a_posteriori_mugal2 * mat_Qxx  # Estimated parameters
        mat_Cldld = s02_a_posteriori_mugal2 * mat_Qldld  # Adjusted observations
        diag_Qxx = np.diag(mat_Qxx)

        # Calculate standard deviations:
        mat_sd_xx = np.sqrt(np.diag(mat_Cxx))
        mat_sd_ldld = np.sqrt(np.diag(mat_Cldld))

        # creating histogram from residuals
        residual_hist, bin_edges = create_hist(mat_v)  # Calculate histogram

        # goodness-of-fit test
        chi_crit, chi_val, chi_test = goodness_of_fit_test(confidence_level_chi_test, dof,
                                                           s02_a_posteriori_mugal2, sig0_mugal ** 2)
        if verbose:
            print(f'')
            print(f'# Goodness-of-fit test results:')
            print(f'# Chi-val  Chi-crt-lower  Chi-crt-upper  Status')
            print(f'  {chi_val:7.3f}  {chi_crit[0]:13.3f}  {chi_crit[1]:13.3f}  {chi_test:s}')
            print(f'')
            print(f'# Histogram of the residuals:')
            print(f'# Lower-edge(mGal)  Upper-edge(mGal)  Frequency')
            for loop_1 in range(len(residual_hist)):
                print(f'  {bin_edges[loop_1]:16.4f}  {bin_edges[loop_1 + 1]:16.4f}  {residual_hist[loop_1]:9.0f}')

        # outlier detection effectiveness (redundancy number)
        diag_Qvv = np.diag(mat_Qvv)
        mat_R = np.diag(mat_P) * diag_Qvv
        # blunder detection parameters calculation
        sv_tau = 1 - confidence_level_tau_test
        std_res, tau_val, tau_crt = tau_criterion_test(diag_Qvv, mat_R, mat_v, s02_a_posteriori_mugal2, dof, sv_tau)

        # #### Store results ####
        g_est_mugal = mat_x[0:number_of_stations, 0]
        sd_g_est_mugal = mat_sd_xx[0:number_of_stations]
        drift_pol_coeff = mat_x[number_of_stations:, 0]
        drift_pol_coeff_sd = mat_sd_xx[number_of_stations:]
        sd_diff_obs_mugal = mat_sd_ldld[:number_of_diff_obs]
        sd_pseudo_obs_mugal = mat_sd_ldld[number_of_diff_obs:]
        v_diff_obs_mugal = mat_v[:number_of_diff_obs, 0]
        v_pseudo_obs_mugal = mat_v[number_of_diff_obs:, 0]
        # std_res ...
        # tau_val ...
        # TODO: Add Tau criterion data here


        # Station related results:
        for idx, stat_name in enumerate(self.observed_stations):
            filter_tmp = self.stat_obs_df['station_name'] == stat_name
            self.stat_obs_df.loc[filter_tmp, 'g_est_mugal'] = g_est_mugal[idx]
            self.stat_obs_df.loc[filter_tmp, 'sd_g_est_mugal'] = sd_g_est_mugal[idx]
        # Calculate differences to estimates:
        self.stat_obs_df['diff_g_est_mugal'] = self.stat_obs_df['g_mugal'] - self.stat_obs_df['g_est_mugal']
        self.stat_obs_df['diff_sd_g_est_mugal'] = self.stat_obs_df['sd_g_mugal'] - self.stat_obs_df['sd_g_est_mugal']

        # Drift parameters:
        survey_name_list = []
        degree_list = []
        coefficient_list = []
        sd_coeff_list = []
        tmp_idx = 0
        for survey_name, setup_df in self.setups.items():
            for degree in range(drift_pol_degree):
                survey_name_list.append(survey_name)
                degree_list.append(degree + 1)  # starts with 1
                coefficient_list.append(drift_pol_coeff[tmp_idx])
                sd_coeff_list.append(drift_pol_coeff_sd[tmp_idx])
        self.drift_pol_df = pd.DataFrame(list(zip(survey_name_list,
                                                  degree_list,
                                                  coefficient_list,
                                                  sd_coeff_list)),
                                         columns=self._DRIFT_POL_DF_COLUMNS)

        # TODO: Add content (estimation results and statistics) to self.setups_diff_df here:
        # Observation related results:
        # self.setups_diff_df.columns
        # Index(['survey_name', 'diff_obs_id', 'station_name_from', 'station_name_to', 'setup_id_from', 'setup_id_to',
        #        'g_diff_mugal', 'sd_g_diff_mugal'], dtype='object')
        #         'sd_g_diff_est_mugal',
        #         'v_diff_mugal'
        for idx, v_diff_mugal in enumerate(v_diff_obs_mugal):
            filter_tmp = self.setups_diff_df['diff_obs_id'] == idx
            self.setups_diff_df.loc[filter_tmp, 'v_diff_mugal'] = v_diff_mugal
            self.setups_diff_df.loc[filter_tmp, 'sd_g_diff_est_mugal'] = sd_diff_obs_mugal[idx]

        # Print results to terminal:
        if verbose:
            print(f'')
            print(f'Station data:')
            pd.options.display.float_format = '{:.3f}'.format
            pd.set_option('display.max_columns', None)
            pd.set_option('display.width', 1000)
            # print here...
            print(self.stat_obs_df[['station_name', 'is_datum',
                                    'g_mugal', 'g_est_mugal', 'diff_g_est_mugal', 'sd_g_mugal', 'sd_g_est_mugal']])

            pd.options.display.float_format = '{:.6f}'.format
            print(f'')
            print(f'Drift polynomial coefficients:')
            print(self.drift_pol_df)

            pd.options.display.float_format = '{:.3f}'.format
            print(f'')
            print(f'Differential observations:')
            for survey_name in survey_names:
                print(f' - Survey: {survey_name}')
                filter_tmp = self.setups_diff_df['survey_name'] == survey_name
                print(self.setups_diff_df.loc[filter_tmp, ['station_name_from', 'station_name_to', 'g_diff_mugal', 'sd_g_diff_mugal', 'sd_g_diff_est_mugal', 'v_diff_mugal']])
                # TODO: Print relevant content of self.setups_diff_df here!

            pd.reset_option('max_columns')

            print(f'')
            print(f'Pseudo observations at datum stations (constraints):')
            # Datum constraints (pseudo observations at datum stations):
            print(f'Station name  sd [µGal]   v [µGal]')
            for idx, station_name in enumerate(datum_stations):
                print(f'{station_name:10}   {sd_pseudo_obs_mugal[idx]:8.3}     {v_pseudo_obs_mugal[idx]:+8.3}')


def tau_criterion_test(diag_Qvv, mat_R, mat_V, ref_var, dof, sv):
    std_res = []
    tau_val = []
    tsd_crt = stats.t.ppf(1 - sv / 2, dof - 1)  # calc. t-distribution crit. val.
    tau_crt = (tsd_crt * np.sqrt(dof)) / np.sqrt(dof - 1 + tsd_crt ** 2)  # calc. tau-criterion crit. val.
    for n in range(len(diag_Qvv)):
        if diag_Qvv[n] == 0:
            std_res.append(0)
        else:
            std_res.append(mat_V[n][0] / np.sqrt(abs(diag_Qvv[n])))  # standardized residual
        if mat_R[n] > 1e-06:
            tau_val.append(abs(std_res[n]) / np.sqrt(ref_var))  # value tested against tau-criterion crit. val.
            # tau_val = abs(std_res)
        else:
            tau_val.append(0)  # if redundancy number is too small blunder won't be detected so why bother?
    return std_res, tau_val, tau_crt


def create_hist(mat_v):
    """Create histogram."""
    residuals = np.ndarray.tolist(mat_v)
    residuals = [item for sublist in residuals for item in sublist]
    hist_residuals, bin_edges = np.histogram(residuals, bins=5)
    return hist_residuals, bin_edges


def goodness_of_fit_test(cf, dof, ref_var, apriori_var):
    """Statistical testing using chi square."""
    # apriori_var = 1
    sv = 1 - cf
    chi_crit_upper = stats.chi2.ppf(1 - sv / 2, dof)  # critical value
    chi_crit_lower = stats.chi2.ppf(sv / 2, dof)  # critical value
    chi_val = dof * ref_var / apriori_var  # tested value
    if chi_crit_lower < chi_val < chi_crit_upper:
        chi_test = 'Passed'
    else:
        chi_test = 'Not passed'
    chi_crit = [chi_crit_lower, chi_crit_upper]
    return chi_crit, chi_val, chi_test
