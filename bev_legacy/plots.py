"""
Plotting utilities for gravtools.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

import numpy as np
import matplotlib.pyplot as plt

# Imports from other gravtools modules:
from bev_legacy import utils


def create_drift_plot(obs_df, stat_df, poly_coef_dict, save_pdf=True, session_name='', path_save_file=''):
    """
    Erstellen eines Drift-Plots.
    :param path_save_file: filepath for saving teh plot
    :param session_name: Name of the observing session (e.g. obs-file name)
    :param obs_df: observation dataframe
    :param stat_df: station dataframe
    :param poly_coef_dict: dict with coefficients of drift polanomial
    :param save_pdf: flag, save plot as PDF?
    :return:
    """
    # Check input arguments
    if save_pdf and (len(session_name) == 0):
        print('ERROR: Name for output PDF file (drift plot) is not defined!')
        exit()

    # Evaluate drift polynomial:
    pol_degree = len(poly_coef_dict)
    poly_coef = utils.prep_polyval_coef(poly_coef_dict)
    dt_h = np.linspace(obs_df.dt_h.min(), obs_df.dt_h.max(), 100)
    yy_mugal = np.polyval(poly_coef, dt_h)

    g_plot_min = 0.0
    g_plot_max = 0.0

    fig, ax = plt.subplots()
    ax.plot(dt_h, yy_mugal, '--', label='Polynom (n={})'.format(pol_degree))
    for pkt_num in stat_df.punktnummer:
        # print(pkt_num)
        df_tmp = obs_df.loc[obs_df['punktnummer'] == pkt_num]
        g_est = stat_df.loc[stat_df['punktnummer'] == pkt_num].g_est_mugal.values[0]
        if stat_df.loc[stat_df['punktnummer'] == pkt_num].is_drift_point.bool():
            label_str = pkt_num + '*'
        else:
            label_str = pkt_num
        y_tmp = df_tmp.g_red_mugal - g_est
        ax.plot(df_tmp.dt_h, y_tmp, 'o', label=label_str)
        if y_tmp.max() > g_plot_max:
            g_plot_max = y_tmp.max()
        if y_tmp.min() < g_plot_min:
            g_plot_min = y_tmp.min()

    gangbreite_mugal = g_plot_max - g_plot_min

    # - Legend and labels:
    plt.legend(loc='best')
    # ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    ax.grid()
    plt.title('Drift-Polynom: {} (Gangbreite: {:4.1f} µGal)'.format(session_name, gangbreite_mugal))
    plt.xlabel('Zeit [h]')
    plt.ylabel('Lesung [µGal]')

    # Save as PDF file:
    if save_pdf:
        plt.savefig(path_save_file + session_name + '_drift.pdf')

    # plt.show(block=True)  # Keep plot open
