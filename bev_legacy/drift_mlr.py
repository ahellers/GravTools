"""
Adjustment of a drift polynomial of dregreee 1 to 3 using multiple linear regression.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

# Imports from other gravtools modules:
from bev_legacy import utils


def calc_drift_corr_mlr(obs_df, stat_df, pol_degree):
    """
    Berechnung eines Polynoms (Grad n) zur Korrektur der Instrumenten-Drift mittels multipler linearer Regression und
    Korrektur der Gravimeter-Lesungen.
    Berechnung analog zur Fortran Routine 'DRIFT2011.for' von Meurers & Ruess.
    :param obs_df:
    :param stat_df:
    :param pol_degree:
    :return:
    """
    # ##### Multiple Linear Regression: #####

    # Set up obs_df with categorical variables for each point:
    df_short = obs_df.loc[:,
               obs_df.columns.intersection(['punktnummer', 'dt_h', 'g_red_mugal'])]  # Only kep the listed columns

    # Create obs_df with categorical dummy variables for each point name:
    prefix = 'is_pkt'
    df_dummy = pd.get_dummies(df_short, prefix=[prefix])
    categorical_variables = [prefix + '_' + str for str in stat_df.punktnummer]

    # Set up target parameters according to polynomial degree:
    target_parameters = []
    for i in range(1, pol_degree + 1):  # i = 1 to polynomial_degree
        col_name = 'dt_{}'.format(i)
        df_dummy[col_name] = df_dummy['dt_h'].pow(i)
        target_parameters.append(col_name)

    selection_list = target_parameters + categorical_variables

    # Calculate regression:
    mlr = LinearRegression(fit_intercept=False)  # No b0 coefficient estimated
    # mlr = LinearRegression()
    mlr.fit(df_dummy[selection_list], df_dummy['g_red_mugal'])

    pol_coef_dict = {i+1: mlr.coef_[i] for i in range(0, pol_degree)}  # Create dict with polynomial coefficients

    # Store corrected reading (stat_df.g_est_mugal):
    tmp_list = mlr.coef_[-stat_df.shape[0]:]  #
    for i, value in enumerate(tmp_list):
        # print(i, value)
        stat_df.iloc[i, 3] = value

    # Calculate drift correction for all observations:
    poly_coef = utils.prep_polyval_coef(pol_coef_dict)
    np.polyval(poly_coef, obs_df.dt_h)
    obs_df['corr_drift_mugal'] = np.polyval(poly_coef, obs_df.dt_h)

    # Calculate statistics:
    for stat in stat_df.punktnummer:

        # Estimate for this station:
        g_est_mugal = stat_df.loc[stat_df['punktnummer'] == stat].g_est_mugal.values[0]  # Schätzwert an der Station

        # Calculate difference between drift corrected reading and the estimates gravity reading at this point (abw):
        obs_df.loc[obs_df['punktnummer'] == stat, 'abw_mugal'] = \
            obs_df.loc[obs_df['punktnummer'] == stat, 'g_red_mugal'] - \
            obs_df.loc[obs_df['punktnummer'] == stat, 'corr_drift_mugal'] - g_est_mugal

        # Mind. 2 Beobachtungen (= Driftpunkt)?
        if stat_df.loc[stat_df['punktnummer'] == stat].is_drift_point.bool():
            num_of_obs = obs_df.loc[obs_df['punktnummer'] == stat].shape[0]
            # Warum "+4"? Eher willkürlich in DRIFT2011 gesetzt, oder? 4 = 2*2.
            stat_df.loc[stat_df['punktnummer'] == stat, 'sig_g_est_mugal'] = \
                np.sqrt((sum(obs_df.loc[obs_df['punktnummer'] == stat, 'abw_mugal'] ** 2) / (num_of_obs - 1)) + 4)

    pol_coef_sig_mugal = np.sqrt(
        sum(obs_df.abw_mugal ** 2) / (obs_df.shape[0] - 1 - sum(stat_df.is_drift_point == False)))

    # Sigma for all estimats with only one observation:
    stat_df.loc[~stat_df['is_drift_point'], 'sig_g_est_mugal'] = np.sqrt(pol_coef_sig_mugal ** 2 + 25)

    results_dict = dict()
    results_dict['pol_coef'] = pol_coef_dict
    results_dict['pol_coef_sig_mugal'] = pol_coef_sig_mugal
    results_dict['stat_df'] = stat_df
    results_dict['obs_df'] = obs_df  # Erweitert mit Verbessungen und Drift-Korekturen je Beobachtung

    return results_dict
