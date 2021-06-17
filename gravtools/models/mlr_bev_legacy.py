"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
Contains classes for the legacy processign scheme of gravimeter observations at BEV using multiple linear regression for
drift adjustment. This scheme also involves the determination of absulte gravity values based on datum stations.
"""

from gravtools.models.lsm import LSM
import pandas as pd
import numpy as np
from scipy import stats
import datetime as dt
import pytz
import copy
from sklearn.linear_model import LinearRegression
from gravtools import settings


def prep_polyval_coef(pol_coef_dict):
    """Preparation of the polynomial coefficients for np.polyval().

    Parameters
    ----------
    pol_coef_dict : dict
        Dictionary with polynomial coefficients from mlr.

    Returns
    -------
    list : list of polynomial coeficients suitable for np.polyval()
    """
    poly_coef = list()
    for degree, value in pol_coef_dict.items():
        poly_coef.append(value)
    poly_coef.reverse()
    poly_coef.append(0)
    return poly_coef


class BEVLegacyProcessing(LSM):
    """Lagacy processing scheme for relative gravimeter observations at BEV.

    First, the instrumental drift is estimated by multiple linear regression. Afterwards, the drift-corrected
    gravimeter readings per station are referred to absolute gravity values at the datum stations.

    Attributes
    ----------

    """

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
        lsm_method = 'MLR_BEV'
        super().__init__(lsm_method, stat_df=stat_df, setups=setups, comment=comment, write_log=write_log)

    @classmethod
    def from_campaign(cls, campaign, comment='', write_log=True):
        """Constructor that generates and populates the LSM-inherited object from a Campaign class object.

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
               verbose=False
               ):  # Confidence level):
        """Run the adjustment.

        Parameters
        ----------
        drift_pol_degree : int, optional (default=1)
            Degree of estimated drift polynomial.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line, e.g. for debugging and testing

        Notes
        -----
        Only one drift polynomial is adjusted in the input data! Hence, it might be problematic to use observations
        of multiple surveys that were taken over a longer time period.
        """

        # ##### 1.) Initial preparations #####

        # - Observations and parameters:
        all_setups_df = None
        observed_stations = []
        survey_names = []
        number_of_observations = 0
        setup_ids = []
        for survey_name, setup_df in self.setups.items():
            observed_stations = observed_stations + setup_df['station_name'].to_list()
            number_of_observations = number_of_observations + len(setup_df)
            survey_names.append(survey_name)
            setup_ids = setup_ids + setup_df['setup_id'].to_list()
            if all_setups_df is None:
                all_setups_df = setup_df.copy(deep=True)
                all_setups_df['survey_name'] = survey_name
            else:
                tmp_setup_df = setup_df.copy(deep=True)
                tmp_setup_df['survey_name'] = survey_name
                all_setups_df = pd.concat([all_setups_df, tmp_setup_df]).reset_index()  # concat dataframes
        observed_stations = list(
            set(observed_stations))  # Unique list of stations => order of stations in matrices!
        number_of_stations = len(observed_stations)
        number_of_surveys = len(self.setups)
        # number_of_parameters = number_of_stations + drift_pol_degree * number_of_surveys  # Total number of parameters to be estimated
        # number_of_diff_obs = number_of_observations - number_of_surveys

        # Check, if setup IDs are unique:
        if len(set(setup_ids)) != len(setup_ids):
            raise AssertionError('Setup IDs are not unique within the campaign!')

        # Write log:
        if verbose or self.write_log:
            tmp_str = f'#### Adjustment log ####\n'
            tmp_str += f'\n'
            tmp_str += f'---- Input data and settings ----\n'
            tmp_str += f'Number of surveys: {number_of_surveys}\n'
            tmp_str += f'Number of stations: {number_of_stations}\n'
            tmp_str += f'Number of observations: {number_of_observations}\n'
            # tmp_str += f'Number of estimated parameters: {number_of_parameters}\n'  # obsolet
            # tmp_str += f'Number of datum stations: {number_of_datum_stations}\n'  # ????????????????????ß
            # tmp_str += f'Degree of freedom: {number_of_diff_obs - number_of_parameters}\n'  # obsolet
            tmp_str += f'\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # Get dataframe with subset of observed stations only:
        filter_tmp = self.stat_df['station_name'].isin(observed_stations)
        self.stat_obs_df = self.stat_df.loc[filter_tmp].copy(deep=True)  # All observed stations

        # ##### 2.) Drift Adjustment #####
        # All obs data is here: all_setups_df['delta_t_h']
        # - delta_t [h] since first obs in session [also deactivated obs count!] => all_setups_df['delta_t_h']
        # Test dataset: "n20200701_1" and "2020-07-01_hermannskogel.txt" (use last obs of each setup! => then equal!)

        # Set up obs_df with categorical variables for each point:
        df_short = all_setups_df[['station_name', 'g_mugal', 'delta_t_h']].copy(deep=True)

        # Create obs_df with categorical dummy variables for each station name:
        prefix = 'is_pkt'
        df_dummy = pd.get_dummies(df_short, prefix=[prefix])
        categorical_variables = [prefix + '_' + str for str in observed_stations]

        # Set up target parameters according to polynomial degree:
        target_parameters = []
        for i in range(1, drift_pol_degree + 1):  # i = 1 to polynomial_degree
            col_name = 'dt_{}'.format(i)
            df_dummy[col_name] = df_dummy['delta_t_h'].pow(i)
            target_parameters.append(col_name)

        selection_list = target_parameters + categorical_variables

        # Calculate regression:
        mlr = LinearRegression(fit_intercept=False)  # No b0 coefficient estimated
        mlr.fit(df_dummy[selection_list], df_dummy['g_mugal'])

        pol_coef_dict = {i + 1: mlr.coef_[i] for i in range(0, drift_pol_degree)}  # Dict with polynomial coefficients

        # Store drift-corrected gravitymeter readings ("g_drift_est_mugal"):
        tmp_list = mlr.coef_[-number_of_stations:]  #
        for idx, station_name in enumerate(observed_stations):
            filter_tmp = self.stat_obs_df['station_name'] == station_name
            self.stat_obs_df.loc[filter_tmp, 'g_drift_est_mugal'] = tmp_list[idx]

        # Calculate drift correction for all observations:
        poly_coef = prep_polyval_coef(pol_coef_dict)

        all_setups_df['corr_drift_mugal'] = np.polyval(poly_coef, all_setups_df['delta_t_h'])

        # Calculate statistics:
        for station_name in observed_stations:
            # Estimate for this station:
            filter_tmp_stat_obs_df = self.stat_obs_df['station_name'] == station_name
            g_drift_est_mugal = self.stat_obs_df.loc[filter_tmp_stat_obs_df, 'g_drift_est_mugal'].values[0]
            # print(station_name, g_est_mugal)
            # print(filter_tmp_stat_obs_df)

            # Calculate difference between drift corrected reading and the estimated gravity reading at this point:
            filter_tmp_all_setups_df = all_setups_df['station_name'] == station_name
            all_setups_df.loc[filter_tmp_all_setups_df, 'abw_mugal'] = \
                all_setups_df.loc[filter_tmp_all_setups_df, 'g_mugal'] - \
                all_setups_df.loc[filter_tmp_all_setups_df, 'corr_drift_mugal'] - g_drift_est_mugal

            # Sigma of "g_drift_est_mugal":
            number_of_observations_of_station = filter_tmp_all_setups_df[filter_tmp_all_setups_df].count()
            if number_of_observations_of_station > 1:  # Get all drift-calibration points
                self.stat_obs_df['is_drift_point'] = True
                # Warum "+4"? Eher willkürlich in DRIFT2011 gesetzt, oder? 4 = 2*2.
                self.stat_obs_df.loc[filter_tmp_stat_obs_df, 'sig_g_drift_est_mugal'] = \
                    np.sqrt((sum(all_setups_df.loc[filter_tmp_all_setups_df, 'abw_mugal'] ** 2) /
                             (number_of_observations_of_station - 1)) + 4)
            else:
                self.stat_obs_df['is_drift_point'] = False

        pol_coef_sig_mugal = np.sqrt(
            sum(all_setups_df['abw_mugal'] ** 2) / (
                    number_of_observations - 1 - sum(self.stat_obs_df['is_drift_point'] == False)))

        # Sigma for all estimats with only one observation:
        # - Why + 25? Just an assumption (measuremenet accuracy is 5 µGal => 5**2 = 25)?
        self.stat_obs_df.loc[~self.stat_obs_df['is_drift_point'],
                             'sig_g_drift_est_mugal'] = np.sqrt(pol_coef_sig_mugal ** 2 + 25)

        # Write log:
        if verbose or self.write_log:
            tmp_str = f'#### Drift adjustment (MLR) ####\n'
            tmp_str += f'\n'
            tmp_str += f'Degree of drift polynomial: {drift_pol_degree}\n'
            tmp_str += f'\n'
            tmp_str += f'Polynomial coefficients (sigma = {pol_coef_sig_mugal:5.2f} µGal)\n'
            for degree, value in pol_coef_dict.items():
                tmp_str += f' - b{degree:} = {value:9.5f} (µGal/h)^{degree} ({value*24/1000:9.5f} (mGal/Tag)^{degree})\n'
            tmp_str += f'\n'
            tmp_str += f'Drift-corrected gravimeter reading per station:\n'
            for idx, row in self.stat_obs_df.iterrows():
                tmp_str += f' - {row["station_name"]:12s}: {row["g_drift_est_mugal"]:12.2f} µGal (sig = {row["sig_g_drift_est_mugal"]:5.2f} µGal)\n'
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # ##### 3.) Determination of absolute gravity values #####
        sig_reading = 2  # [µGal], Assumption of the reading error, independent of the gravimeter type, etc.

        # Get datum stations:
        filter_tmp = self.stat_obs_df['is_datum'] == True
        stat_df_datum = self.stat_obs_df.loc[filter_tmp, :].copy(deep=True)
        stat_df_datum['sd_g_mugal'] = stat_df_datum['sd_g_mugal'].astype(float)  # Convert to float
        stat_df_datum['g_mugal'] = stat_df_datum['g_mugal'].astype(float)  # Convert to float
        number_of_datum_stations = len(stat_df_datum)

        # Differences between drift-corrected readings and datum g-values:
        stat_df_datum['g0_mugal'] = stat_df_datum['g_mugal'] - stat_df_datum['g_drift_est_mugal']
        stat_df_datum['sig_g0_mugal'] = np.sqrt(stat_df_datum['sig_g_drift_est_mugal'] ** 2
                                                + stat_df_datum['sd_g_mugal'] ** 2
                                                + sig_reading ** 2)

        # Calculation of the weighted mean of all g0:
        stat_df_datum['go_p'] = 1 / (stat_df_datum['sig_g0_mugal'] ** 2)
        g0_mean_mugal = sum(stat_df_datum['go_p'] * stat_df_datum['g0_mugal']) / sum(
            stat_df_datum['go_p'])  # Weighted mean of g0 values (weights: go_p)
        if number_of_datum_stations == 1:  # Only one OEGSN stations measured
            sig_g0_mean_mugal = stat_df_datum['sd_g_mugal']
        else:
            sig_g0_mean_mugal = np.sqrt(
                ((sum(stat_df_datum['g0_mugal'] ** 2) - ((sum(stat_df_datum['g0_mugal'])) ** 2 / number_of_datum_stations)) /
                 (number_of_datum_stations * (number_of_datum_stations - 1))) + (1 / (sum(stat_df_datum['go_p']))))

        # Calculation of the absolute g values at all observed stations:
        self.stat_obs_df['g_est_mugal'] = self.stat_obs_df['g_drift_est_mugal'] + g0_mean_mugal
        self.stat_obs_df['sd_g_est_mugal'] = np.sqrt(sig_g0_mean_mugal ** 2 + self.stat_obs_df['sig_g_drift_est_mugal'] ** 2)
        self.stat_obs_df['diff_g_est_mugal'] = self.stat_obs_df['g_est_mugal'] - self.stat_obs_df['g_mugal']
        self.stat_obs_df['diff_sd_g_est_mugal'] = self.stat_obs_df['sd_g_est_mugal'] - self.stat_obs_df['sd_g_mugal']
        # self.stat_obs_df['g_est_full_mugal'] = self.stat_obs_df['g_est_mugal'] + 9.8e8  # Get absolute gravity value at station [µGal]

        # Add g0 and sig_go to "self.stat_obs_df":
        self.stat_obs_df['g0_mugal'] = None
        self.stat_obs_df['sig_g0_mugal'] = None
        self.stat_obs_df.update(stat_df_datum[['g0_mugal', 'sig_g0_mugal']])

        # Write log:
        if verbose or self.write_log:
            tmp_str = f'\n'
            tmp_str += f'#### Absolute gravity determination ####\n'
            tmp_str += f'\n'
            tmp_str += f'Number of datum stations: {number_of_datum_stations}\n'
            tmp_str += self.stat_obs_df[['station_name', 'g_est_mugal', 'sd_g_est_mugal',
                                         'diff_g_est_mugal']].to_string(
                index=False,float_format=lambda x: '{:.1f}'.format(x))
            if verbose:
                print(tmp_str)
            if self.write_log:
                self.log_str += tmp_str

        # ##### 4.) Save data/infos for later use #####
        self.observed_stations = observed_stations
        self.drift_polynomial_degree = drift_pol_degree
        # self.number_of_estimates = number_of_parameters
        self.number_of_stations = number_of_stations
        # self.degree_of_freedom = dof

        # Drift dataframe:
        survey_name_list = []
        degree_list = []
        coefficient_list = []
        sd_coeff_list = []
        coeff_unit_list = []
        for degree, pol_coeff in pol_coef_dict.items():
            survey_name_list.append(None)  # One drift polynomial estimated for all surveys in campaign
            degree_list.append(degree)
            coefficient_list.append(pol_coeff)
            sd_coeff_list.append(None)  # No standard deviation available
            coeff_unit_list.append(f'µGal/h^{degree}')
        self.drift_pol_df = pd.DataFrame(list(zip(survey_name_list,
                                                  degree_list,
                                                  coefficient_list,
                                                  sd_coeff_list,
                                                  coeff_unit_list)),
                                         columns=self._DRIFT_POL_DF_COLUMNS)

        # Rename columns to be compatible with the table view model:
        all_setups_df.rename(columns={'epoch_dt': 'ref_epoch_dt'}, inplace=True)
        self.setup_obs_df = all_setups_df

        # Not used here => =None as initialized!
        # self.sig02_a_priori = sig0_mugal
        # self.scaling_factor_datum_observations = scaling_factor_datum_observations
        # self.confidence_level_chi_test = confidence_level_chi_test
        # self.confidence_level_tau_test = confidence_level_chi_test
        # self.s02_a_posteriori = s02_a_posteriori_mugal2

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


# pd.options.display.float_format = '{:.2f}'.format
# pd.set_option('display.max_columns', 500)
# pd.set_option('display.width', 1000)

# TODO: Write output files
# TODO: Create drift plot
