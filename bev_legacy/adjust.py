"""
<description>
Author: Andreas Hellerschmied
"""

import sys
import numpy as np
from scipy.linalg import block_diag

# Imports from other gravtools modules:
from gravtools import settings
from bev_legacy import init


def gauss_markoff(obs_df, stat_df, pol_degree=1):
    """
    Adjustment of gravity meter observations using a Gauss-Markoff model.
    Linear problem => Without linearization - no a priori g values required at non-datum stations.
    Estimated parameters:

     - offset and linear drift per instrument
     - gravity at observed stations

    :return:
    """

    # Options:
    # datum_stations_list = ['0-059-21']
    datum_stations_list = []

    sig0_mugal = 1  # a priori standard deviaiton of unit wheigt of observations [µGal]
    sd_cond_scaling_factor = 1e-1  # Scaling factor for condition "observation2"
    g_sd_defaut_mugal = 10  # Default value for the SD of obsrvations

    # If Sd id not provided in the observation file, use the default value:
    obs_df.loc[obs_df.g_sd_mugal.isna(), 'g_sd_mugal'] = g_sd_defaut_mugal

    stat_df['is_datum'] = False
    obs_df['is_datum'] = False
    if not datum_stations_list:  # is not empty
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

    # obs_df = obs_df.sort_values('is_datum')  # Sort: [obs of new stations; obs of datum stations]
    obs_df = obs_df.reset_index()  # enumeration as index instead of time-tag
    # stat_df = stat_df.sort_values('is_datum')  # Sort: [obs of new stations; obs of datum stations]

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
    num_est_para_per_instr = pol_degree + 1  # number of estimated parameters per gravimeter
    num_est_instr_para = num_est_para_per_instr * num_gravimeter  # total number of estimated instrument parameters
    num_obs = obs_df.shape[0]
    num_stations = stat_df.shape[0]
    num_datum_stations = stat_df[stat_df['is_datum'] == True].shape[0]  # all ÖSGN stations
    num_new_stations = num_stations - num_datum_stations
    num_estimates = num_est_instr_para + num_stations
    degree_of_freedom = num_obs + num_datum_stations - num_estimates  # num_datum_stations = num of conditions

    # Preallocate matices and vectors:
    G = np.zeros([num_obs, num_stations])
    D = np.zeros([num_obs, num_est_instr_para])
    Z = np.zeros([num_datum_stations, num_est_instr_para])
    H = np.zeros([num_datum_stations, num_stations])
    L = np.empty(num_obs)
    L = L.reshape(-1, 1)
    L[:] = np.nan
    h = np.empty(num_datum_stations)
    h = h.reshape(-1, 1)
    h[:] = np.nan

    # Weight matrix:
    P = np.diag((sig0_mugal / obs_df['g_sd_mugal']) ** 2)

    # Set up individual matrices:
    for i_obs, values in obs_df.iterrows():
        # print(values['punktnummer'])
        idx_stat = stat_num_dict[values['punktnummer']]  # Get index of observed station
        G[i_obs, idx_stat] = 1
        L[i_obs] = values['g_red_mugal']
        for grav_type, idx_start in gravimeter_type_dict.items():
            for i_deg in range(0, pol_degree + 1):
                # print(i_deg)
                D[i_obs, idx_start * num_est_para_per_instr + i_deg] = values['dt_h'] ** i_deg  # a_k

    # Conditions:
    i_col = 0
    h_sig_mugal = np.empty(num_datum_stations)
    h_sig_mugal = h_sig_mugal.reshape(-1, 1)
    h_sig_mugal[:] = np.nan
    sd_obs_min_mugal = obs_df['g_sd_mugal'].min()
    for i_stat, values in stat_df.iterrows():
        if values['is_datum']:
            H[i_col, stat_num_dict[values['punktnummer']]] = 1
            h[i_col, 0] = values['g_oesgn_mugal']
            h_sig_mugal[i_col, 0] = sd_obs_min_mugal * sd_cond_scaling_factor
            i_col = i_col + 1

    HP = np.diagflat((sig0_mugal / h_sig_mugal) ** 2)
    P = block_diag(P, HP)
    L = np.concatenate((L, h), axis=0)  # vertical

    A1 = np.concatenate((G, -D), axis=1)  # horizontal
    A2 = np.concatenate((H, Z), axis=1)  # horizontal
    A = np.concatenate((A1, A2), axis=0)  # vertical
    ATP = A.transpose().dot(P)
    N = ATP.dot(A)

    print(' - Defree of Freedom {}'.format(degree_of_freedom))
    print(' - Number of estimates: {}'.format(num_estimates))
    print(' - Rank of N: {}'.format(np.linalg.matrix_rank(N)))
    print(' - Condition of N: {}'.format(np.linalg.cond(N)))
    print(' - Determinate of N: {}'.format(np.linalg.det(N)))
    print(' - Rank of A: {}'.format(np.linalg.matrix_rank(A)))
    print(' - Size of N: {} x {}'.format(N.shape[0], N.shape[0]))
    print(' - Rangdefizit: {}'.format(N.shape[0] - np.linalg.matrix_rank(N)))

    b = ATP.dot(L)
    N_inv = np.linalg.inv(N)
    xd = N_inv.dot(b)  # Ausgeglichene Unbekannte

    v = A.dot(N_inv).dot(ATP).dot(L) - L  # Verbessungeruen, lt. AG1 Equ. (2.40)
    Ld = L + v  # Ausgeglichenen Beobachtungen, lt. AG1 Equ. (2.6)

    # Stochastisches Modell a posteriori:
    Qxdxd = N_inv
    QLdLd = A.dot(N_inv).dot(A.transpose())

    # Empirische Varianz der Gewichtseinheit (a posteriori):
    s02 = (v.transpose().dot(P).dot(v)) / degree_of_freedom  # AG1 (2.99)

    # Kovarianz-Matritzen a posteriori:
    Cxdxd = s02 * Qxdxd
    CLdLd = s02 * QLdLd

    sig_xd = np.sqrt(np.diag(Cxdxd))
    sig_LD = np.sqrt(np.diag(CLdLd))

    # Gewichtsreziproke nach Ansermet AG1 (2.74)
    u = np.trace(P.dot(QLdLd))  # Muss Anzahl der Unbekannten ergeben!
    if not u.round(10) == num_estimates:
        print('WARNING: Gewichtsreziproke nach Ansermet nicht erfüllt (u = {})!'.format(u))

    # Hauptprobe:
    # Funktionaler Zusammenhang zwischen den ausgeglichenen Beobachteungen und den ausgeglichenen Unbekannten genau
    # genug erfüllt?
    w_probe = A.dot(xd) - Ld
    print('Hauptprobe: Max. Widerspuch (Absolutwert) = {}'.format(abs(w_probe).max()))

    print('Empirische Varianz der Gewichtseinheit (a posteriori) = {}'.format(s02.item()))

    # Ausgeglichene Beobachtungen (inkl. Std.Abw. a posteriori):
    print('Ausgeglichene Beobachtungen (Std.Abw.):')
    for i_obs, values in obs_df.iterrows():
        if values.is_datum:
            print(' - {:5.2f}h @ {}: {:12.3f} ({:6.3f}) µGal; v = {:8.3f} µGal *'.format(values['dt_h'],
                                                                                         values['punktnummer'],
                                                                                         Ld[i_obs].item(),
                                                                                         sig_LD[i_obs].item(),
                                                                                         v[i_obs].item()))
        else:
            print(' - {:5.2f}h @ {}: {:12.3f} ({:6.3f}) µGal; v = {:8.3f} µGal'.format(values['dt_h'],
                                                                                       values['punktnummer'],
                                                                                       Ld[i_obs].item(),
                                                                                       sig_LD[i_obs].item(),
                                                                                       v[i_obs].item()))

    # Retrive and print results:
    print('Geschätzte Schwere an Stationen (Std.Abw.):')
    for stat, num in stat_num_dict.items():
        print(' - {}: {:10.3f} µGal ({:6.3f} µGal)'.format(stat, xd[num].item(), sig_xd[num].item()))
    print('Sonstige Parameter (Std.Abw.):')
    for i in range(num_stations, len(xd)):
        print(' - a{}: {:15.3f} ({:6.3f})'.format(i - num_stations, xd[i].item(), sig_xd[i].item()))

    pass


def gauss_markoff_omc(obs_df, stat_df, pol_degree=1):
    """
    Adjustment of gravity meter observations using a Gauss-Markoff model.

    WITH A RRIORI Vales for the estimates
    Estimated parameters:

     - offset and linear drift per instrument
     - gravity at observed stations

    :return:
    """

    # Options:
    # datum_stations_list = ['0-059-21']
    datum_stations_list = []

    sig0_mugal = 1  # a priori standard deviaiton of unit wheigt of observations [µGal]
    sd_cond_scaling_factor = 1e-0  # Scaling factor for conditions ("fictional observations")
    g_sd_defaut_mugal = 10  # Default value for the SD of observations

    # If Sd id not provided in the observation file, use the default value:
    obs_df.loc[obs_df.g_sd_mugal.isna(), 'g_sd_mugal'] = g_sd_defaut_mugal

    stat_df['is_datum'] = False
    obs_df['is_datum'] = False
    if not datum_stations_list:  # is not empty
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

    # obs_df = obs_df.sort_values('is_datum')  # Sort: [obs of new stations; obs of datum stations]
    obs_df = obs_df.reset_index()  # enumeration as index instead of time-tag
    # stat_df = stat_df.sort_values('is_datum')  # Sort: [obs of new stations; obs of datum stations]

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
    num_est_para_per_instr = pol_degree + 1  # number of estimated parameters per gravimeter
    num_est_instr_para = num_est_para_per_instr * num_gravimeter  # total number of estimated instrument parameters
    num_obs = obs_df.shape[0]
    num_stations = stat_df.shape[0]
    num_datum_stations = stat_df[stat_df['is_datum'] == True].shape[0]  # all ÖSGN stations
    num_new_stations = num_stations - num_datum_stations
    num_estimates = num_est_instr_para + num_stations
    degree_of_freedom = num_obs + num_datum_stations - num_estimates  # num_datum_stations = num of conditions

    # Preallocate matices and vectors:
    G = np.zeros([num_obs, num_stations])
    D = np.zeros([num_obs, num_est_instr_para])
    Z = np.zeros([num_datum_stations, num_est_instr_para])
    H = np.zeros([num_datum_stations, num_stations])
    L = np.empty(num_obs)
    L = L.reshape(-1, 1)
    L[:] = np.nan
    h = np.empty(num_datum_stations)
    h = h.reshape(-1, 1)
    h[:] = np.nan

    # ### Model computed observations: ###

    # A priori gravity at observed stations:
    stat_df['g_0_mugal'] = stat_df['g_oesgn_mugal']  # Auf ÖSGN Werte setzen, wenn möglich
    # - Wenn kein apriori g-Wert vorhanden => Auf mittelwert der beobachteten ÖSGN-Punkte setzen
    #   - Besser wäre: Interpoliren, g-Wert vom nächsten Punkt mit VG, etc.
    stat_df.loc[stat_df['g_0_mugal'].isna(), 'g_0_mugal'] = stat_df['g_oesgn_mugal'].mean()
    g0 = stat_df.g_0_mugal.to_numpy()
    g0 = g0.reshape(-1, 1)  # vertical vector

    # A priori coefficients of drit polynomial:
    c0 = np.zeros([pol_degree, 1])

    # a priori gravimeter reading offset a0:
    # - Difference between gravimeter reading and a priori gravity at first observed Datum station
    first_datum_obs = obs_df[obs_df['is_datum'] == True].iloc[0]
    a0 = stat_df[stat_df['punktnummer'] == first_datum_obs.punktnummer].g_oesgn_mugal - first_datum_obs.g_red_mugal
    a0 = a0.to_numpy().reshape(-1, 1)  # vertical vector

    # Vector of a priori parameters X0:
    X0 = np.concatenate((g0, a0, c0), axis=0)  # vertical

    # Weight matrix:
    P = np.diag((sig0_mugal / obs_df['g_sd_mugal']) ** 2)

    # Set up individual matrices:
    for i_obs, values in obs_df.iterrows():
        # print(values['punktnummer'])
        idx_stat = stat_num_dict[values['punktnummer']]  # Get index of observed station
        G[i_obs, idx_stat] = 1
        L[i_obs] = values['g_red_mugal']
        for grav_type, idx_start in gravimeter_type_dict.items():
            for i_deg in range(0, pol_degree + 1):
                # print(i_deg)
                D[i_obs, idx_start * num_est_para_per_instr + i_deg] = values['dt_h'] ** i_deg  # a_k

    # Conditions:
    i_col = 0
    h_sig_mugal = np.empty(num_datum_stations)
    h_sig_mugal = h_sig_mugal.reshape(-1, 1)
    h_sig_mugal[:] = np.nan
    sd_obs_min_mugal = obs_df['g_sd_mugal'].min()
    for i_stat, values in stat_df.iterrows():
        if values['is_datum']:
            H[i_col, stat_num_dict[values['punktnummer']]] = 1
            h[i_col, 0] = values['g_oesgn_mugal']
            # h[i_col, 0] = 0
            h_sig_mugal[i_col, 0] = sd_obs_min_mugal * sd_cond_scaling_factor
            i_col = i_col + 1

    HP = np.diagflat((sig0_mugal / h_sig_mugal) ** 2)
    P = block_diag(P, HP)
    L = np.concatenate((L, h), axis=0)  # vertical

    # Modell-Matrix inkl. Bedingungen (fiktive Beobachtungen) in A2
    A1 = np.concatenate((G, -D), axis=1)  # horizontal
    A2 = np.concatenate((H, Z), axis=1)  # horizontal
    A = np.concatenate((A1, A2), axis=0)  # vertical

    # O-C vector (AG1 (2.31)):
    L0 = A.dot(X0)
    l = L - L0  # AG2: (2.12)

    ATP = A.transpose().dot(P)
    N = ATP.dot(A)

    print(' - Defree of Freedom {}'.format(degree_of_freedom))
    print(' - Number of estimates: {}'.format(num_estimates))
    print(' - Rank of N: {}'.format(np.linalg.matrix_rank(N)))
    print(' - Condition of N: {}'.format(np.linalg.cond(N)))
    print(' - Determinate of N: {}'.format(np.linalg.det(N)))
    print(' - Rank of A: {}'.format(np.linalg.matrix_rank(A)))
    print(' - Size of N: {} x {}'.format(N.shape[0], N.shape[0]))
    print(' - Rangdefizit: {}'.format(N.shape[0] - np.linalg.matrix_rank(N)))

    b = ATP.dot(l)  # AG1 (2.37)
    N_inv = np.linalg.inv(N)
    xd = N_inv.dot(b)  # Ausgeglichene Unbekannte # AG1 (2.37)
    Xd = X0 + xd  # AG1 (2.7)

    v = A.dot(N_inv).dot(ATP).dot(l) - l  # Verbessungeruen, lt. AG1 Equ. (2.40)
    Ld = L + v  # Ausgeglichenen Beobachtungen, lt. AG1 Equ. (2.6)

    # Stochastisches Modell a posteriori:
    Qxdxd = N_inv
    QLdLd = A.dot(N_inv).dot(A.transpose())

    # Empirische Varianz der Gewichtseinheit (a posteriori):
    s02 = (v.transpose().dot(P).dot(v)) / degree_of_freedom  # AG1 (2.99)

    # Kovarianz-Matritzen a posteriori:
    Cxdxd = s02 * Qxdxd
    CLdLd = s02 * QLdLd

    sig_xd = np.sqrt(np.diag(Cxdxd))
    sig_Ld = np.sqrt(np.diag(CLdLd))

    # Gewichtsreziproke nach Ansermet AG1 (2.74)
    u = np.trace(P.dot(QLdLd))  # Muss Anzahl der Unbekannten ergeben!
    if not u.round(10) == num_estimates:
        print('WARNING: Gewichtsreziproke nach Ansermet nicht erfüllt (u = {})!'.format(u))

    # Hauptprobe:
    # Funktionalr Zusammenhang zwischen den ausgeglichenen Beobachteungen und den ausgeglichenen Unbekannten genau
    # genug erfüllt?
    w_probe = A.dot(Xd) - Ld
    print('Hauptprobe: Max. Widerspuch (Absolutwert) = {}'.format(abs(w_probe).max()))

    print('Empirische Varianz der Gewichtseinheit (a posteriori) = {}'.format(s02.item()))

    # Ausgeglichene Beobachtungen (inkl. Std.Abw. a posteriori):
    print('Ausgeglichene Beobachtungen (Std.Abw.):')
    for i_obs, values in obs_df.iterrows():
        if values.is_datum:
            print(' - {:5.2f}h @ {}: {:12.3f} ({:6.3f}) µGal; v = {:8.3f} µGal *'.format(values['dt_h'],
                                                                                         values['punktnummer'],
                                                                                         Ld[i_obs].item(),
                                                                                         sig_Ld[i_obs].item(),
                                                                                         v[i_obs].item()))
        else:
            print(' - {:5.2f}h @ {}: {:12.3f} ({:6.3f}) µGal; v = {:8.3f} µGal'.format(values['dt_h'],
                                                                                       values['punktnummer'],
                                                                                       Ld[i_obs].item(),
                                                                                       sig_Ld[i_obs].item(),
                                                                                       v[i_obs].item()))

    # Retrive and print results:
    print('Geschätzte Schwere an Stationen (Std.Abw.):')
    for stat, num in stat_num_dict.items():
        print(' - {}: {:10.3f} µGal ({:6.3f} µGal)'.format(stat, Xd[num].item(), sig_xd[num].item()))
    print('Sonstige Parameter (Std.Abw.):')
    for i in range(num_stations, len(xd)):
        print(' - a{}: {:15.3f} ({:6.3f})'.format(i - num_stations, Xd[i].item(), sig_xd[i].item()))

    pass


# def vermittelnde_beob_mit_bedingungen(obs_df, stat_df, pol_degree=1):
#     """
#     Adjustment of gravity meter observations using a vermittelnde Beobachtungen mit Bedingungsgleichungen
#     See: AG1, pp.16-17
#     Linear problem => Without linearization - no a priori g values required at non-datum stations.
#     Estimated parameters:
#      - offset and linear drift per instrument
#      - gravity at observed stations
#     :return:
#     """
#     pass


def main(path_oesgn_table, name_oesgn_table, path_obs_file, name_obs_file, out_path):

    obs_df, stat_df, df_oesgn, obs_info_dict = init.read_and_prep_data(path_oesgn_table,
                                                                       name_oesgn_table,
                                                                       path_obs_file,
                                                                       name_obs_file)

    # Gauss-Markoff model
    gauss_markoff(obs_df, stat_df, 1)

    # Gauss - Markoff model with computed observatons and linearization:
    gauss_markoff_omc(obs_df, stat_df, 1)

    # Vermittelnde Beobachtungen mit Bedgingungsgleichungen:
    # adjust_vermittelnde_beob_mit_bedingungen(obs_df, stat_df, 1)

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
