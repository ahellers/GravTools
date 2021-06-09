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
import pytz
import copy
import matplotlib.pyplot as plt
# from abc import ABC
from gravtools import settings

from gravtools.settings import SURVEY_DATA_SOURCE_TYPES, STATION_DATA_SOURCE_TYPES, GRAVIMETER_ID_BEV, \
    TIDE_CORRECTION_TYPES, DEFAULT_GRAVIMETER_ID_CG5_SURVEY, REFERENCE_HEIGHT_TYPE, NAME_OBS_FILE_BEV, \
    PATH_OBS_FILE_BEV, BEV_GRAVIMETER_TIDE_CORR_LOOKUP, GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m
from gravtools.const import VG_DEFAULT
from gravtools.models.exceptions import FileTypeError
from gravtools.CG5_utils.cg5_survey import CG5Survey

# Using abstract base classes, see e.g.:
#  - https://www.python-course.eu/python3_abstract_classes.php
#  - https://stackoverflow.com/questions/5133262/python-abstract-base-class-init-initializion-or-validation


class LSM:
    """Base class for LSM classes.

    Attributes
    ----------
    lsm_method : str
        Defines the adjustment method. Has to b listed in :py:obj:`gravtools.settings.ADJUSTMENT_METHODS`.
    stat_df : :py:obj:`gravtools.Station.stat_df`
        The station dataframe contains all relevant station data.
    setups : dict of pandas DataFrames
        The setups dictionary contains all observation data used for the adjustment. The keys of the dictionary
        are the survey names (str) and the items are pandas dataframes containing the observation data (see
        :py:obj:`gravtool.Survey.setup_df`)
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
        """
    def __init__(self, lsm_method, stat_df, setups, comment='', write_log=True):
        """
        Parameters
        ----------
        lsm_method : str
            Defines the adjustment method. Has to b listed in :py:obj:`gravtools.settings.ADJUSTMENT_METHODS`.
        stat_df : :py:obj:`gravtools.Station.stat_df`
            The station dataframe contains all relevant station data.
        setups : dict of pandas DataFrames
            The setups dictionary contains all observation data used for the adjustment. The keys of the dictionary
            are the survey names (str) and the items are pandas dataframes containing the observation data (see
            :py:obj:`gravtool.Survey.setup_df`)
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
        self.drift_pol_df = None  # DataFrame that contains the estimated parameters of the drift polynomials an
        # their statistics for each survey

        # Estimation settings:
        self.drift_polynomial_degree = None
        self.sig0_a_priori = None  # A priori standard deviation of unit weight
        self.scaling_factor_datum_observations = None
        self.confidence_level_chi_test = None
        self.confidence_level_tau_test = None

        # General statistics:
        self.number_of_stations = None
        self.number_of_datum_stations = None
        self.number_of_estimates = None
        self.degree_of_freedom = None

        # A posteriori statistics:
        self.s02_a_posteriori = None  # A posteriori variance of unit weight

        # Matrices:
        self.Cxx = None # Co-variance matrix of estimated parameters

    @property
    def time_str(self):
        """Return time of lsm adjustment as formatted string."""
        return self.init_time.strftime("%Y-%m-%d, %H:%M:%S")


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
    _SETUP_DIFF_COLUMNS = (
        'survey_name',
        'ref_epoch_dt',
        'diff_obs_id',
        'station_name_from',
        'station_name_to',
        'setup_id_from',
        'setup_id_to',
        'g_diff_mugal',
        'sd_g_diff_mugal',
        'sd_g_diff_est_mugal',
        'v_diff_mugal',
        'w_diff_mugal',
        'r_diff_obs',
        'tau_test_result',
    )

    # Column names of self.drift_pol_df:
    _DRIFT_POL_DF_COLUMNS = (
        'survey_name',
        'degree',
        'coefficient',
        'sd_coeff',
        'coeff_unit',
    )

    def __init__(self, stat_df, setups, comment='', write_log=True):
        """
        Parameters
        ----------
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
                        setups[survey_name] = survey.setup_df

        # Check if setup data is available:
        if len(setups) == 0:
            raise AssertionError(f'Setup data of campaign "{campaign.campaign_name}" does not contain observations!')

        # Initialize and return LSM object:
        return cls(campaign.stations.stat_df, setups, comment=comment, write_log=write_log)

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
            tmp_str = f'#### Adjustment log ####\n'
            tmp_str += f'\n'
            tmp_str += f'---- Input data and settings ----\n'
            tmp_str += f'Number of surveys: {number_of_surveys}\n'
            tmp_str += f'Number of stations: {number_of_stations}\n'
            tmp_str += f'Number of differential observations: {number_of_diff_obs}\n'
            tmp_str += f'Number of estimated parameters: {number_of_parameters}\n'
            tmp_str += f'Number of datum stations: {number_of_datum_stations}\n'
            tmp_str += f'Degree of freedom: {number_of_diff_obs - number_of_parameters}\n'
            tmp_str += f'\n'
            tmp_str += f'Degree of drift polynomial: {drift_pol_degree}\n'
            tmp_str += f'A priori std. deviation of unit weight: {sig0_mugal}\n'
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
        g_diff_obs_mugal_list = []
        station_name_from_list = []
        station_name_to_list = []
        setup_id_from_list = []
        setup_id_to_list = []
        sd_g_diff_obs_mugal_list = []
        diff_obs_id_list = []
        survey_names_list = []
        ref_epoch_dt_list = []  # Reference epochs (datetime objects) of differential observations

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
                    # epoch_from = previous_row['epoch_unix']  # [sec]
                    epoch_from = previous_row['delta_t_h']  # [hours]
                    # epoch_to = row['epoch_unix']  # [sec]
                    epoch_to = row['delta_t_h']  # [hours]
                    ref_epoch_dt = previous_row['epoch_dt'] + (row['epoch_dt'] - previous_row['epoch_dt']) / 2  # mean

                    # Populate matrices and vectors:
                    mat_L0[(diff_obs_id, 0)] = g_diff_mugal
                    mat_sig_ll0[diff_obs_id]  = sd_g_diff_mugal ** 2
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
            # print(f' - Station {row["station_name"]:10s} (row index {index:d})')
            datum_station_id += 1
            station_name = row['station_name']
            station_id = self.observed_stations.index(station_name)
            mat_Ac[datum_station_id, station_id] = 1  # Partial derivative
            mat_Lc[(datum_station_id, 0)] = row['g_mugal']  # g for datum definition
            sd_mugal_for_weighting = row['sd_g_mugal'] * scaling_factor_datum_observations
            mat_sig_llc[datum_station_id] = sd_mugal_for_weighting ** 2

        # Set up all required matrices:
        mat_A = np.vstack((mat_A0, mat_Ac))  # Eq. (16)
        mat_L = np.vstack((mat_L0, mat_Lc))  # Eq. (16)
        mat_sig_ll = np.diag(np.hstack((mat_sig_ll0, mat_sig_llc)))
        mat_Qll = mat_sig_ll / (sig0_mugal**2)
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

        # A posteriori variance of unit weight s02:
        dof = mat_A.shape[0] - mat_A.shape[1]  # degree of freedom
        par_r = mat_v.T @ mat_P @ mat_v  # = v^T * P * v
        if dof == 0:
            s02_a_posteriori_mugal2 = par_r[0][0]
            # dof = 0 should not be the case here!
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
        mat_sd_vv = np.sqrt(np.diag(mat_Cvv))  # A posteriori SD of residuals

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
            tmp_str += f'# Lower-edge(mGal)  Upper-edge(mGal)  Frequency\n'
            for loop_1 in range(len(residual_hist)):
                tmp_str += f'  {bin_edges[loop_1]:16.4f}  {bin_edges[loop_1 + 1]:16.4f}  {residual_hist[loop_1]:9.0f}\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # outlier detection effectiveness (redundancy components)
        diag_Qvv = np.diag(mat_Qvv)
        # mat_R = np.diag(mat_P) * diag_Qvv
        # mat_R is exactly the same as "mat_r"

        # Redundanzanteile (redundancy components):
        # - AG II, pp. 66-71
        mat_r = np.diag(mat_Qvv @ mat_P)

        # Standardisierte Versesserungen (AG II, p. 66)
        # - Normalverteilt mit Erwartwarungswert = 0 (wie Verbesserungen)
        # - Standardabweichung = 1 (standardisiert)
        mat_w = mat_v[:,0] / mat_sd_vv

        # Tau test for outlier detection:
        alpha_tau = 1 - confidence_level_tau_test
        tau_test_result, tau_critical_value = tau_test(mat_w=mat_w, dof=dof, alpha=alpha_tau, mat_r=mat_r)

        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'# Tau-test results:\n'
            tmp_str += f'Critical value: {tau_critical_value:1.3f}\n'
            tmp_str += f' - Number of detected outliers: {tau_test_result.count("failed")}\n'
            tmp_str += f' - Number low redundancy component: {tau_test_result.count("r too small")}\n'
            tmp_str += f'\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # blunder detection parameters calculation
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
        tmp_idx = 0
        for survey_name, setup_df in self.setups.items():
            for degree in range(drift_pol_degree):
                survey_name_list.append(survey_name)
                degree_list.append(degree + 1)  # starts with 1
                coefficient_list.append(drift_pol_coeff[tmp_idx])  # * (3600**(degree + 1))  # [µGal/h]
                sd_coeff_list.append(drift_pol_coeff_sd[tmp_idx])
                coeff_unit_list.append(f'µGal/h^{degree + 1}')
                tmp_idx += 1
        self.drift_pol_df = pd.DataFrame(list(zip(survey_name_list,
                                                  degree_list,
                                                  coefficient_list,
                                                  sd_coeff_list,
                                                  coeff_unit_list)),
                                         columns=self._DRIFT_POL_DF_COLUMNS)

        # TODO: Add content (estimation results and statistics) to self.setup_obs_df here:
        # Observation related results:
        # self.setup_obs_df.columns
        # Index(['survey_name', 'diff_obs_id', 'station_name_from', 'station_name_to', 'setup_id_from', 'setup_id_to',
        #        'g_diff_mugal', 'sd_g_diff_mugal'], dtype='object')
        #         'sd_g_diff_est_mugal',
        #         'v_diff_mugal'
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
                                                            'v_diff_mugal']].to_string(index=False,
                                                                                       float_format=lambda x: '{:.1f}'.format(x))
            tmp_str += f'\n\n'
            tmp_str += f' - Pseudo observations at datum stations (constraints):\n'
            tmp_str += f'Station name  sd [µGal]   v [µGal]   w [µGal]   r [0-1]    Tau test result\n'
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

    def create_drift_plot_matplotlib(self):
        """Create a drift plot with matplotlib."""

        # Prep data:
        setup_df = self.setups['20200701a'].copy(deep=True)
        stat_obs_df_short = self.stat_obs_df.loc[:, ['station_name', 'g_est_mugal', 'sd_g_est_mugal']]
        setup_df = pd.merge(setup_df, stat_obs_df_short, on="station_name")
        setup_df['g_plot_mugal'] = setup_df['g_mugal'] - setup_df['g_est_mugal']
        setup_df.sort_values(by='delta_t_h', inplace=True)

        # Evaluate drift polynomial:
        coeff_list = self.drift_pol_df['coefficient'].to_list()
        coeff_list.reverse()
        coeff_list.append(0)
        t_min_h = 0
        t_max_h = setup_df['delta_t_h'].max()
        dt_h = np.linspace(t_min_h, t_max_h, 100)
        drift_polynomial_mugal = np.polyval(coeff_list, dt_h)

        # !!! Due to the differential observations, the constant bias (N0) of the gravity reading cannot be estimated!
        # In order to draw the drift polynomial function w.r.t. the gravity meter observations (for the sake of visual
        # assessment of the drift function), the const. bias N0 is approximated, see below.
        offset_mugal = setup_df['g_plot_mugal'].mean()
        yy_mugal = drift_polynomial_mugal - drift_polynomial_mugal.mean() + offset_mugal

        # plot
        fig, ax = plt.subplots()
        for station_name in setup_df['station_name'].unique():
            label_str = station_name
            delta_t_h = setup_df.loc[setup_df['station_name'] == station_name, 'delta_t_h']
            g_plot_mugal = setup_df.loc[setup_df['station_name'] == station_name, 'g_plot_mugal']
            ax.plot(delta_t_h, g_plot_mugal, 'o', label=label_str)

        ax.plot(dt_h, yy_mugal, 'k--', label='drift function')
        # - Legend and labels:
        plt.legend(loc='best')
        ax.grid()
        plt.title(f'Drift Polynomial (unknown vertical offset!)')
        plt.xlabel('time [h]')
        plt.ylabel('gravity reading [µGal]')

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


def bin_redundacy_components(mat_r):
    """Bin redundancy components for easier interpretatzion.

    See: Skriptum AG2 (Navratil, TU Wien), p. 70.
    """
    number_obs_r_equal_0 = 0  # No error control: Errors cannot be detected
    number_obs_r_0_to_03 = 0  # Bad error control
    number_obs_r_03_to_07 = 0  # Good error control
    number_obs_r_07_1 = 0  # Very good error control, but observation may be redundant
    number_obs_r_equal_1 = 0  # Observatiuon is redundant

    results_dict = {}  # '': ''

    return number_obs_r_equal_0, number_obs_r_0_to_03


def tau_test(mat_w, dof, alpha, mat_r):
    """Tau-criterion test for outlier detection of a least-squares adjustment.

    See: Pope (1976): The statistics of residuals and the detection of outliers.

    Parameters
    ----------
    mat_w: np.array of floats
        Vector with standardized post-fit residuals of all observations in the least-squares adjustment.
        They are computed as post-fit residuls divides by their standard deviations.
    dof: int
        Degree of freedom of the least-squares adjustment.
    alpha: float (0 to 1)
        Significance level.
    mat_r: np.array of floats
        Vector with redundancy components derived in the least-squares adjustment for each observation.
    """
    # Critical value:
    tsd_crt = stats.t.ppf(1 - alpha / 2, dof - 1)  # calc. t-distribution crit. val.
    tau_crt = (tsd_crt * np.sqrt(dof)) / np.sqrt(dof - 1 + tsd_crt ** 2)  # Critical value (Pope, 1976, Eq. (6))
    tau_test_result = []
    for idx, w in enumerate(mat_w):
        # Check, if redancy component is larger than threshold:
        if mat_r[idx] > settings.R_POPE_TEST_TRESHOLD:
            if np.abs(w) > tau_crt:
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

def goodness_of_fit_test(cf, dof, a_posteriori_variance_of_unit_weight, a_priori_variance_of_unit_weight):
    """Statistical testing using chi square.

    Notes
    -----
    This "global model test" is dewscribed by Caspari (1987): Concepts of network and deformation analysis.
    pp. 6-8 and pp.68-69.

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
    # a_priori_variance_of_unit_weight = 1
    alpha = 1 - cf  # Significance level = Probability of commiting a type 1 error (H0 wrongly dismissed)
    chi_crit_upper = stats.chi2.ppf(1 - alpha / 2, dof)  # critical value
    chi_crit_lower = stats.chi2.ppf(alpha / 2, dof)  # critical value
    chi_val = dof * a_posteriori_variance_of_unit_weight / a_priori_variance_of_unit_weight  # tested value
    #TODO: Why is there an upper AND a lower critical value? In the literatur only an upper critical value ist defined!
    if chi_crit_lower < chi_val < chi_crit_upper:
        chi_test_status = 'Passed'
    else:
        chi_test_status = 'Not passed'
    chi_crit = [chi_crit_lower, chi_crit_upper]
    return chi_crit, chi_val, chi_test_status


if __name__ == '__main__':
    test2 = LSMDiff('a', 'b')  # This test will cause an error due to invalid data types!
    pass

# TODO: Save relevant estimation settings and matrices/vectors in LSM object for later analysis and documentation!
# TODO: Global model test: Why is there an upper and lower critical value? => In literarture only upper!
# TODO: Add information on tests (global model AND Tau-test) to log string!
# TODO: Check the documentation/docstrings!

# TODO: Add pseudo observation results to results df and label them as constraints!
# - Problem: contriaints do not have a reference epoch.
#   - Lösung:
#     (1) Alle Resultate in results_df
#     (2) Constraint-Beob mit col "is_constraint" flaggen
#     (3) {get_model_data_df_for_plotting} verwenden, um plotable Daten zu bekommen!
#     (4) contraint Beob. in observations results table farblich kennzeichnen!

# TODO: Plot co-variance matrix for estimates

# TODO: Properly implement drift plot option (see comments above at method: create_drift_plot_matplotlib)
# => A slier in the gui would be nice to adjust the unknown offset!

# TODO: Treat redundancy components roperly and add determination of inner and outer reliability (AG2, pp. 70-72)
# r: In obs results table die einzelnen obs nach der Kategorisierung (in settings definiert) auf p. 70 einteilen!
