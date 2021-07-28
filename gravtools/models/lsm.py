"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
Contains classes for least-squares adjustment of gravimeter-surveys.
"""

import numpy as np
from scipy import stats
import datetime as dt
import pytz
import copy
# from abc import ABC
from gravtools import settings


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

    @property
    def get_log_string(self):
        """Returns the log string."""
        return self.log_str


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


# TODO: Save relevant estimation settings and matrices/vectors in LSM object for later analysis and documentation!
# TODO: Global model test: Why is there an upper and lower critical value? => In literarture only upper!
# TODO: Add information on tests (global model AND Tau-test) to log string!
# TODO: Check the documentation/docstrings!

# TODO: Plot co-variance matrix for estimates

# TODO: Treat redundancy components roperly and add determination of inner and outer reliability (AG2, pp. 70-72)
# r: In obs results table die einzelnen obs nach der Kategorisierung (in settings definiert) auf p. 70 einteilen!
