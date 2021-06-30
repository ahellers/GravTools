"""
Legacy processing pipeline for gravity observations at BEV, emulating the legacy fortran routines by Meurers & Ruess.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

import sys
import numpy as np

# Imports from other gravtools modules:
from bev_legacy import settings
from bev_legacy import init, drift_mlr, output, plots


def calc_abs_g(stat_df):
    """
    Berechnung absluter Schwerewerte bezogen auf das ÖSGN, ausgehend von Drift-korrigierten Lesungen eines
    Relativgravimeters. Für die Auswertung einzelner Tagesmessungen.
    Berechnung analog zu Fortran Code SCHWAUS2016.FOR von D. Ruess ()
    :param stat_df:
    :return: stat_df (mit Ergebnissen)
    """

    # ### Parameters: ###
    # Lesungsabhängiger Fehler (als konstant angenommen, unabhängig vom Betrag der Ablesung und vom Gravimeter-Typ):
    sig_reading = 2  # [µGal]

    # df, only containingÖSGN stations:
    df_stat_oesgn = stat_df.loc[stat_df['is_oesgn']].copy()  # Ensure to make a copy instead of a sliced view of stat_df

    # Differences between drift-corrected readings and OESGN g-values:
    df_stat_oesgn['g0_mugal'] = df_stat_oesgn['g_oesgn_mugal'] - df_stat_oesgn['g_est_mugal']
    df_stat_oesgn['sig_g0_mugal'] = np.sqrt(df_stat_oesgn['sig_g_est_mugal'] ** 2
                                            + df_stat_oesgn['g_sig_oesgn_mugal'] ** 2
                                            + sig_reading ** 2)

    # Calculation of the weighted mean of all go:
    df_stat_oesgn['go_p'] = 1 / (df_stat_oesgn['sig_g0_mugal'] ** 2)
    xm_mugal = sum(df_stat_oesgn['go_p'] * df_stat_oesgn['g0_mugal']) / sum(
        df_stat_oesgn['go_p'])  # Weighted mean of g0 values (weights: go_p)
    if df_stat_oesgn.shape[0] == 1:  # Only one OEGSN stations measured
        sig_xm_mugal = df_stat_oesgn['g_sig_oesgn_mugal']
    else:
        n = df_stat_oesgn.shape[0]  # Number of observed ÖSGN stations
        sig_xm_mugal = np.sqrt(((sum(df_stat_oesgn['g0_mugal'] ** 2) - ((sum(df_stat_oesgn['g0_mugal'])) ** 2 / n)) /
                                (n * (n - 1))) + (1 / (sum(df_stat_oesgn['go_p']))))

    # Calculation of the absolute g values at all observed stations:
    stat_df['g_abs_mugal'] = stat_df['g_est_mugal'] + xm_mugal
    stat_df['sig_g_abs_mugal'] = np.sqrt(sig_xm_mugal ** 2 + stat_df['sig_g_est_mugal'] ** 2)
    stat_df['verb_mugal'] = stat_df['g_abs_mugal'] - stat_df['g_oesgn_mugal']
    stat_df['g_abs_full_mugal'] = stat_df['g_abs_mugal'] + 9.8e8  # Get absolute gravity value at station [µGal]

    # Add g0 and sig_go to stat_df:
    stat_df.update(df_stat_oesgn[['g0_mugal', 'sig_g0_mugal']])

    return stat_df, xm_mugal, sig_xm_mugal


def drift_schwaus(obs_df, stat_df, polynomial_degree, instrument_id, name_obs_file, out_path='',
                  skalierung=None, date_str='', timezone_str=''):
    """
    Nachbau aller Brechnungen auf DRIFT2011.FOR und SCHWAUS2016.FOR.

    :param obs_dict:
    :param obs_df:
    :param stat_df:
    :return:
    """

    # ### Drift-Korrektur ###
    # Analog zu Fortran Code DRIFT2011.FOR
    # Berechnungen mittels multipler linearer Regression
    # - 1.) Polynom vom Grad n=1,2,3 an die Lesungen anpassen
    # - 2.) Drift-korrigierte Lesungen bestimmen

    # ### Polynom zur Korrektur der Instrumentendrift vom Grad n anpassen  ###
    if settings.VERBOSE:
        print('Schätzung der Drift mit Polynom vom Grad {}:'.format(polynomial_degree))
    results_dict = drift_mlr.calc_drift_corr_mlr(obs_df, stat_df, polynomial_degree)
    stat_df = results_dict['stat_df']
    if settings.VERBOSE:
        # Print results:
        print(' - Polynom Koeffizeinten (sig = {:5.2f} µGal):'.format(results_dict['pol_coef_sig_mugal']))
        for degree, value in results_dict['pol_coef'].items():
            print('    b{deg} = {value:9.5f} µGal/h^{deg} ({value2:9.5f} mGal/Tag^{deg})'.format(deg=degree, value=value, value2=value*24/1000))
        print(' - Korrigierte Gravimeter-Lesungen je Station:')
        for i, values in stat_df.iterrows():
            print('    {}: {:12.2f} µGal (sig = {:5.2f} µGal)'.format(values['punktnummer'], values['g_est_mugal'],
                                                                      values['sig_g_est_mugal']))

    # ### Schwere bezogen auf das ÖSGN brechnen ###
    # - Analog zu Fortan Code SCHWAUS2016.FOR
    # - Berechnete absolute Schwere an den beobachteten Stationen bezogen auf das Höheniveau des Festpunktes!
    if settings.VERBOSE:
        print('Berechnung absoluter Schwere-Werte:')
        print(' - Lagerung der Drift-korrigierten Lesungen auf {} ÖSGN Stationen'.format(stat_df[stat_df['is_oesgn']
                                                                                                 == True].shape[0]))
    stat_df, g0_mean_mugal, sig_g0_mean_mugal = calc_abs_g(stat_df)
    if settings.VERBOSE:
        # Print results:
        print(' - Ergebnisse:')
        for i, values in stat_df.iterrows():
            print('    {}: g = {:12.2f} µGal (sig = {:5.2f} µGal), verb. = {:7.2f} µGal'.format(values['punktnummer'],
                                                                                                values[
                                                                                                    'g_abs_full_mugal'],
                                                                                                values[
                                                                                                    'sig_g_abs_mugal'],
                                                                                                values['verb_mugal']))

    # ### Write file for NSDB input ###
    path_name_nsb_file = out_path + name_obs_file + '.nsb'
    if settings.VERBOSE:
        print('NSDB Input Datei schreiben ({})'.format(path_name_nsb_file))
    output.write_nsb_file(stat_df, instrument_id, path_name_nsb_file)

    # ### Create Plot ###
    path_name_drift_plot = out_path + name_obs_file
    if settings.VERBOSE:
        print('Drift-Plot erstellen')
        if settings.FLAG_SAVE_DRIFT_PLOT_PDF:
            print(' - Speichern unter: {}_drift.pdf'.format(path_name_drift_plot))

    plots.create_drift_plot(obs_df, stat_df, results_dict['pol_coef'],
                            save_pdf=settings.FLAG_SAVE_DRIFT_PLOT_PDF,
                            session_name=name_obs_file,
                            path_save_file=out_path)

    # ### Write protocol file ###
    if settings.FLAG_CREATE_SCHWAUS_PROTOCOL:
        if settings.VERBOSE:
            print('Berechnungsprotokoll erstellen ({})'.format(out_path + name_obs_file + '_prot.txt'))
        output.write_schwaus_protcol(obs_df=obs_df,
                                     stat_df=stat_df,
                                     pol_coef=results_dict['pol_coef'],
                                     pol_coef_sig_mugal=results_dict['pol_coef_sig_mugal'],
                                     session_name=name_obs_file,
                                     path_save_file=out_path,
                                     instrument_id=instrument_id,
                                     skalierung=skalierung,
                                     g0_mean_mugal=g0_mean_mugal,
                                     sig_g0_mean_mugal=sig_g0_mean_mugal)


def main(path_oesgn_table, name_oesgn_table, path_obs_file, name_obs_file, out_path):

    obs_df, stat_df, df_oesgn, obs_info_dict = init.read_and_prep_data(path_oesgn_table,
                                                                       name_oesgn_table,
                                                                       path_obs_file,
                                                                       name_obs_file)

    # ### DRIFT2011 und SCHWAUS2016 ###
    drift_schwaus(obs_df, stat_df,
                  polynomial_degree=obs_info_dict['polynomial_degree'],
                  instrument_id=obs_info_dict['instrument_id'],
                  name_obs_file=obs_info_dict['filename'],
                  out_path=out_path,
                  skalierung=obs_info_dict['skalierung'],
                  date_str=obs_info_dict['date_str'],
                  timezone_str=obs_info_dict['timezone_str'])

    print('...fertig!')



# Run as standalone program:
if __name__ == "__main__":

    if len(sys.argv) == 1:
        print('Eingangsparameter aus options.py bezogen (Default-Parameter).')
        path_oesgn_table = settings.PATH_OESGN_TABLE
        name_oesgn_table = settings.NAME_OESGN_TABLE
        path_obs_file = settings.PATH_OBS_FILE_BEV
        name_obs_file = settings.NAME_OBS_FILE_BEV
        out_path = settings.OUT_PATH

    elif len(sys.argv) == 2:  # Name of observation file as input argument
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
