"""
Gravity network adjustment described by Reilly (1970).
Output functions of gravtools.
Author: Andreas Hellerschmied
"""

import sys
import numpy as np

# Imports from other gravtools modules:
from bev_legacy import settings
from bev_legacy import init


def adjust_reilly1970(obs_df, stat_df):
    """
    Adjustment of gravity meter observations according to Reilly 1970.
    Estimated parameters:

     - offset and linear drift per instrument
     - correction of the calibration constant of each instrument
     - gravity at observed stations

    :return:
    """

    # datum_stations_list = ['0-059-21']
    datum_stations_list = []
    stat_df['is_datum'] = False
    obs_df['is_datum'] = False
    if not datum_stations_list:  # is not empty
        # flag_test = False
        # if flag_test:
        # All ÖSGN stations are datum stations
        stat_df['is_datum'] = stat_df.is_oesgn
        obs_df.loc[~obs_df.g_oesgn_mugal.isna(), 'is_datum'] = True
    else:
        # Use datum stations from list:
        for stat in datum_stations_list:
            stat_df.loc[stat_df['punktnummer'] == stat, 'is_datum'] = True
            obs_df.loc[obs_df['punktnummer'] == stat, 'is_datum'] = True
            # Check, if all datum stationsade ÖSGN stations:
            if not stat_df[stat_df['punktnummer'] == stat].is_oesgn.bool():
                print('ERROR: the datum station {} is no ÖSGN station!'.format(stat))

    obs_df = obs_df.sort_values('is_datum')  # Sort: [obs of new stations; obs of datum stations]
    obs_df = obs_df.reset_index()  # enumeration as index instead of time-tag
    stat_df = stat_df.sort_values('is_datum')  # Sort: [obs of new stations; obs of datum stations]

    # Lookup-table for station numbers:
    # - e.g. indicating the position of station-related matrix entries
    stat_num_dict = stat_df['punktnummer'].to_dict()
    stat_num_dict = dict((v, k) for k, v in stat_num_dict.items())  # Exchange keys and values

    # Lookup-table for gravimeters:
    gravimeter_type_dict = {}
    gravimeter_type = obs_df.gravimeter_typ.unique()
    for num, grav_type in enumerate(gravimeter_type):
        gravimeter_type_dict[grav_type] = num

    # Get parameters:
    num_gravimeter = obs_df.gravimeter_typ.unique().shape[0]
    num_est_para_per_instr = 2  # number of estimated parameters per gravimeter
    num_est_instr_para = num_est_para_per_instr * num_gravimeter  # total number of estimated instrument parameters
    num_obs = obs_df.shape[0]
    num_stations = stat_df.shape[0]
    num_datum_stations = stat_df[stat_df['is_datum'] == True].shape[0]  # all ÖSGN stations
    num_new_stations = num_stations - num_datum_stations

    # Preallocate matices and vectors:
    P = np.zeros([num_obs, num_stations])
    Q = np.zeros([num_datum_stations, num_stations])
    T = np.zeros([num_obs, num_est_instr_para])
    Z = np.zeros([num_datum_stations, num_est_instr_para])
    y = np.empty(num_obs)
    y = y.reshape(-1, 1)
    y[:] = np.nan

    # Set up vectors:
    h = stat_df[stat_df.is_datum].g_oesgn_mugal.to_numpy()
    h = h.reshape(-1, 1)

    # Set up individual matrices:
    for i_obs, values in obs_df.iterrows():
        if not values['is_datum']:  # Observed station is not a datum station
            # print(values['punktnummer'])
            idx_stat = stat_num_dict[values['punktnummer']]  # Get index of observed station
            P[i_obs, idx_stat] = 1
            y[i_obs] = values['g_red_mugal']
        else:  # Observed station is a datum station
            y[i_obs] = values['g_red_mugal'] - values['g_oesgn_mugal']
        for grav_type, idx_start in gravimeter_type_dict.items():
            T[i_obs, idx_start * num_est_para_per_instr] = 1  # a_k
            T[i_obs, idx_start * num_est_para_per_instr + 1] = values['dt_h']  # b_k
            if num_est_para_per_instr == 3:
                T[i_obs, idx_start * num_est_para_per_instr + 2] = values['g_mugal'] - values['tide_mugal']  # beta_k

    i_col = 0
    for i_stat, values in stat_df.iterrows():
        if values['is_datum']:
            Q[i_col, stat_num_dict[values['punktnummer']]] = 1
            i_col = i_col + 1

    # combine matrics and build equation system:
    # - Formulation by Reilly 1970:
    # N11 = P.transpose().dot(P) + Q.transpose().dot(Q)
    # N12 = -P.transpose().dot(T)
    # N21 = -T.transpose().dot(P)
    # N22 = T.transpose().dot(T)
    # N1 = np.concatenate((N11, N12), axis=1)
    # N2 = np.concatenate((N21, N22), axis=1)
    # N = np.concatenate((N1, N2), axis=0)
    #
    # B11 = P.transpose().dot(y) + Q.transpose().dot(h)
    # B21 = -T.transpose().dot(y)
    # B = np.concatenate((B11, B21), axis=0)
    #
    # # Solve equation system:
    # N_inv = np.linalg.inv(N)
    # x = N_inv.dot(B)

    # - Standard formulation (same results as Reilly 1970)
    A = np.concatenate((np.concatenate((P, -T), axis=1), np.concatenate((Q, Z), axis=1)), axis=0)
    l = np.concatenate((y, h), axis=0)
    N = A.transpose().dot(A)

    print(' - Rank of N: {}'.format(np.linalg.matrix_rank(N)))
    print(' - Rank of A: {}'.format(np.linalg.matrix_rank(A)))
    print(' - Size of N: {} x {}'.format(N.shape[0], N.shape[0]))
    print(' - Rangdefizit: {}'.format(N.shape[0] - np.linalg.matrix_rank(N)))

    N_inv = np.linalg.inv(N)
    x = N_inv.dot(A.transpose().dot(l))

    # Retrive and print results:
    #
    print('Geschätzte Schwere an Stationen:')
    for stat, num in stat_num_dict.items():
        print(' - {}: {:10.3f} µGal'.format(stat, x[num].item()))
    print('Sonstige Parameter:')
    for i in range(num_stations, len(x)):
        print(' - a{}: {:15.3f}'.format(i - num_stations + 1, x[i].item()))


def main(path_oesgn_table, name_oesgn_table, path_obs_file, name_obs_file, out_path):

    obs_df, stat_df, df_oesgn, obs_info_dict = init.read_and_prep_data(path_oesgn_table,
                                                                       name_oesgn_table,
                                                                       path_obs_file,
                                                                       name_obs_file)

    # Adjustment
    adjust_reilly1970(obs_df, stat_df)

    print('...finished!')


# Run as standalone program:
if __name__ == "__main__":

    if len(sys.argv) == 1:
        print('Eingangsparameter aus options.py bezogen (Default-Parameter).')
        path_oesgn_table = settings.PATH_OESGN_TABLE
        name_oesgn_table = settings.NAME_OESGN_TABLE
        path_obs_file = settings.PATH_OBS_FILE_BEV
        name_obs_file = settings.NAME_OBS_FILE_BEV
        out_path = settings.OUT_PATH

    elif len(sys.argv) == 2:  # Name of obervation file as input argument
        name_obs_file = sys.argv[1]
        # Default-parameters from options-file:
        path_oesgn_table = settings.PATH_OESGN_TABLE
        name_oesgn_table = settings.NAME_OESGN_TABLE
        path_obs_file = settings.PATH_OBS_FILE_BEV
        out_path = settings.OUT_PATH

    else:
        print('Error: Invalid number of input arguments!')
        exit()

    # Run init module:
    main(path_oesgn_table,
         name_oesgn_table,
         path_obs_file,
         name_obs_file,
         out_path)

else:
    # not run as standalone program, but as module
    pass