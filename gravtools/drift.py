"""
Programm zur Schätzung der Gerätedrift bei Relativgravimetern
Autor: Andreas Hellerschmied
Date: 2020-05-08
Benützung: python3 drift.py <Beobachungsdatei (Feldbuch) in Sub-Odner ./data>
"""

# Imports:

import sys
import pandas as pd
import numpy as np
from scipy.linalg import block_diag
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ####################
# # Options          #
# ####################
dic_gravimeter_hoehenbezug_korrektur_m = {
    'CG3': -0.211,
    'CG5': -0.211,
}
# Instrumenten-IDs in der Messdatei einem Instrument zuweisen:
dict_gravimeter_id_obs_file = {
    '3': 'CG3',
    '5': 'CG5',
}
# Instrumenten-IDs in der Messdatei einem Instrumenten-KZ zuweisen:
dict_gravimeter_KZ_obs_file = {
    '3': 'C',
    '5': 'C',
    '9': 'D',
    '90': 'D',
    '900': 'D',
    '51': 'D',
    '510': 'D',
    '500': 'W',
}

vg_default = 308.6  # muGal/m
verbous = True
flag_save_drift_plot_pdf = True


# #####################
# ##### Functions #####
# #####################

def read_oesgn_table(filename, filepath=''):
    """
    Read ÖSGN Table and return pandas dataframe containing the data.
    :return:
    """
    widths = [
        10,  # Punktnummer im ÖSGN
        24,  # Anmerkung
        8,  # Geograph. Breite [deg] 34
        8,  # Geograph. Länge [deg] 42
        8,  # Höhe [m] 50
        7,  # g [??] 58
        3,  # mittlerer Fehler von g [?] 65
        4,  # Vertikalgradient der Schwere [µGal/m] 68
        6,  # Datum 73
        12,  # Punktidentität 79
    ]

    column_names = [
        'punktnummer',
        'anmerkungen',
        'breite_deg',
        'laenge_deg',
        'hoehe_m',
        'g_oesgn_mugal',
        'g_sig_oesgn_mugal',
        'vg_mugal',
        'date',
        'identitaet'
    ]

    df = pd.read_fwf(filepath + filename, widths=widths, header=None)
    df.columns = column_names

    return df


def read_obs_file(filename, filepath=''):
    """
    Einlesen der (konvertierten) Messdateien mit dem folgenden Format:
    Zeile 1: Für skalierte bzw. nicht skalierte Werte ('Y' bzw. 'N')
    Zeile 2: Gerätenummer (3 bzw. 5 für CG3 bzw. CG5)
    Zeile 3: Grad des geschätzten Polynoms (1, 2 oder 3)
    Zeile 4: Datum (YYYY MM DD)
    Zeile 5: Zeitzone (UTC, MEZ, OEZ)
    Zeile 6 bis n: Punktnummer(x-xxx-xx), Messzeit (hh.mm), g [mGal], Breite [°] (bb.b), Länge [°] (ll.l)
    Zeile n+1: 'end'

    return: dict with observation data
    """

    obs_data = dict()

    # Check, if the input file contains an additional columns with the tidal correction:
    flag_tide_corr = False
    if filename.split('_')[-1] == 'tideCorr':
        flag_tide_corr = True

    # Check, if the input file contains an additional columns with the tidal correction:
    flag_sd = False
    if filename.split('_')[-1] == 'sd':
        flag_sd = True


    # Read 5 header lines
    num_of_header_lines = 5
    with open(filepath + filename) as myfile:
        head = [next(myfile) for x in range(num_of_header_lines)]
    skalierung = head[0][:-1]
    instrument_id = head[1].split()[0]
    polynomial_degree = int(head[2][:-1])
    date_str = head[3][:-1]
    timezone_str = head[4][:-1]

    # Read obs data => pandas dataframe
    widths = [
        11,  # Punktnummer im ÖSGN
        5,  # Zeit
        9,  # g [mGal]
        7,  # Geograph. Breite [deg]
        6,  # Geograph. Länge [deg]
    ]
    column_names = [
        'punktnummer',
        'zeit',
        'g_mugal',
        'dhb_m',
        'dhf_m',
    ]

    if flag_tide_corr:
        widths.append(8)
        column_names.append('tide_mugal')

    if flag_sd:
        widths.append(6)
        column_names.append('g_sd_mugal')



    df = pd.read_fwf(filepath + filename, widths=widths, header=None, names=column_names, skiprows=5, skipfooter=1,
                     dtype={'zeit': object})
    if not flag_sd:
        df['g_sd_mugal'] = np.nan
    if not flag_tide_corr:
        df['tide_mugal'] = np.nan

    df['zeit'] = pd.to_datetime(date_str + ' ' + df['zeit'] + ' ' + timezone_str, format='%Y %m %d %H.%M %Z')
    df = df.set_index('zeit')  # set 'zeit' as index
    df['g_mugal'] = df['g_mugal'] * 1e3  # Conversion from mGal to muGal
    df['tide_mugal'] = df['tide_mugal'] * 1e3  # Conversion from mGal to muGal
    df['g_sd_mugal'] = df['g_sd_mugal'] * 1e3  # Conversion from mGal to muGal
    df['dhb_m'] = df['dhb_m'] * 1e-2  # Conversion from cm to m
    df['dhf_m'] = df['dhf_m'] * 1e-2  # Conversion from cm to m

    df['gravimeter_typ'] = dict_gravimeter_id_obs_file[instrument_id]

    # Store results in dict:
    obs_data['skalierung'] = skalierung
    obs_data['instrument_id'] = instrument_id
    obs_data['polynomial_degree'] = polynomial_degree
    obs_data['date_str'] = date_str
    obs_data['timezone_str'] = timezone_str
    obs_data['obs_df'] = df

    return obs_data


def prep_polyval_coef(pol_coef_dict):
    """
    Vorbereitung der Polynomkoeffizienten für np.polyval().
    :param pol_coef_dict: Dictionary with polynomial coefficients from mlr.
    :return: poly_coef: list of polynomial coeficients suitable for np.polyval()
    """
    poly_coef = list()
    for degree, value in pol_coef_dict.items():
        poly_coef.append(value)
    poly_coef.reverse()
    poly_coef.append(0)
    return poly_coef


def prep_obs_df(obs_df, df_oesgn):
    """
    Dataframe mit Beobachtungen aufbereiten und mit Informationen aus dm ÖSGN erweitern.
    :param obs_df: dataframe mit Messungen
    :param df_oesgn: dataframe mit ÖSGN Daten
    :return: dataframe mit Messungen (erweitert)
    """

    # Mit ÖSGN Tabelle mergen und die Messepoche als Index setzen:
    obs_df = obs_df.reset_index().merge(df_oesgn.drop(columns=['anmerkungen', 'identitaet', 'date', 'laenge_deg',
                                                               'hoehe_m', 'breite_deg']), on='punktnummer',
                                        how='left').set_index('zeit')  # Index = Zeit

    # Initialize columns that are to be filled later on:
    obs_df['g_red_mugal'] = np.nan  # Auf Festpunkt reduzierter Schwerewert [µGal]
    obs_df['corr_drift_mugal'] = np.nan  # Auf Drift-Korrektur für Lesung [µGal]
    obs_df['abw_mugal'] = np.nan  # Abweichung zwischen Drift korrigierter Lesung und Schätzwert g_est am Punkt

    # Zeitliche Referenz der Messungen:
    obs_df['unix_time'] = pd.to_numeric(obs_df.index)  # Unix Zeit (absolut) = ns since 1970
    obs_df['dt_h'] = (obs_df.index - obs_df.index.min()).seconds / 3600.0  # Zeitspannt zw. Messung i = 1 bis n und
    # erster Messung i=0

    # Falls VG nicht über die ÖSGN Tabelle gegeben, default Wert verwenden:
    obs_df.loc[obs_df.vg_mugal.isna(), 'vg_mugal'] = vg_default

    return obs_df


def corr_ref_heights_instrument(obs_df):
    """
    Korrektur der Höhenoffsets zwischen Instrument-Oberkante (bei Messung notiert) und Boden bzw. Festpunkt (dHF, dHB)
    um den Höhenoffset zwischen Instrumenten-Oberkante und Messwerk.
    :param obs_df: dataframe mit Messungen
    :return: dataframe mit Messungen (erweitert)
    """

    # ### Je nach Gerätetyp (CG5, CG3) die Bezugs-Höhen (dhb, dhf) anpassen ###:
    for grav_typ, corr_m in dic_gravimeter_hoehenbezug_korrektur_m.items():
        obs_df.loc[obs_df['gravimeter_typ'] == grav_typ, 'dhb_m'] = \
            obs_df.loc[obs_df['gravimeter_typ'] == grav_typ, 'dhb_m'] + corr_m
        obs_df.loc[obs_df['gravimeter_typ'] == grav_typ, 'dhf_m'] = \
            obs_df.loc[obs_df['gravimeter_typ'] == grav_typ, 'dhf_m'] + corr_m
    return obs_df


def red_g_obs_vg(obs_df):
    """
    Gravimeter-Lesungen mittels Vertikalgadienten auf die Höhe des Festpunktes reduzieren (mittels dHF)
    :param obs_df: dataframe mit Messungen
    :return: obs_df: dataframe mit Messungen (mit auf die Festpunkte reduzierten Schweremessungen)
    """
    obs_df['g_red_mugal'] = obs_df['g_mugal'] + obs_df['vg_mugal'] * obs_df['dhf_m']
    return obs_df


def get_station_df(obs_df, df_oesgn):
    """
    Erzeugt einen dataframe mit allen Stationen um stationsspezifische Daten und Ergebnisse zu speichern.
    :param obs_df: dataframe mit Messungen
    :param df_oesgn: dataframe mit ÖSGN Daten
    :return: stat_df: dict with observation data
    """

    # get unique station list and number of observations per station:
    tmp_series = obs_df.punktnummer.value_counts()
    stat_df = tmp_series.to_frame(name='num_of_obs')
    stat_df = stat_df.reset_index().rename(columns={'index': 'punktnummer'})

    # Initialize columns that are to be filled later on:
    stat_df['is_drift_point'] = stat_df.num_of_obs > 1  # More than 1 obs => Drif point!
    stat_df['g_est_mugal'] = np.nan
    stat_df['sig_g_est_mugal'] = np.nan
    stat_df['dhb_m'] = np.nan  # Depends on the setup of each measurement => Only the difference (dhb-dhf) is const.!
    stat_df['dhf_m'] = np.nan  # Depends on the setup of each measurement => Only the difference (dhb-dhf) is const.!
    stat_df['g_abs_mugal'] = np.nan
    stat_df['sig_g_abs_mugal'] = np.nan
    stat_df['verb_mugal'] = np.nan
    stat_df['g_abs_full_mugal'] = np.nan  # absolute value in µGal incl. 9.8e8

    # Add ÖSGN data:
    stat_df = stat_df.merge(df_oesgn.drop(columns=['anmerkungen', 'identitaet', 'date']), on='punktnummer', how='left')
    stat_df['is_oesgn'] = ~stat_df.g_oesgn_mugal.isna()
    stat_df['date'] = obs_df.index[0].date()  # date for observation

    # Assign dhb and dhf:
    for stat in stat_df.punktnummer:
        stat_df.loc[stat_df['punktnummer'] == stat, 'dhb_m'] = \
            obs_df[obs_df['punktnummer'] == stat].dhb_m[-1]  # Assign dhb of the last obervation
        stat_df.loc[stat_df['punktnummer'] == stat, 'dhf_m'] = \
            obs_df[obs_df['punktnummer'] == stat].dhf_m[-1]  # Assign dhb of the last obervation

    return stat_df


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

    # Create obs_df with categorical dummy variables for each Pointname:
    prefix = 'is_pkt'
    df_dummy = pd.get_dummies(df_short, prefix=[prefix])
    categorical_variables = [prefix + '_' + str for str in stat_df.punktnummer]

    # Set up target parameters accoding to polynomial degree:
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

    pol_coef_dict = {i: mlr.coef_[i] for i in range(0, pol_degree)}  # Create dict with polynomila coefficients

    # Store corrected reading:
    tmp_list = mlr.coef_[-stat_df.shape[0]:]  #
    for i, value in enumerate(tmp_list):
        # print(i, value)
        stat_df.iloc[i, 3] = value

    # Calculate drift correction for all observations:
    poly_coef = prep_polyval_coef(pol_coef_dict)
    np.polyval(poly_coef, obs_df.dt_h)
    obs_df['corr_drift_mugal'] = np.polyval(poly_coef, obs_df.dt_h)

    # Calculate statistics:
    for stat in stat_df.punktnummer:

        # Estimate for this station:
        g_est_mugal = stat_df.loc[stat_df['punktnummer'] == stat].g_est_mugal.values[0]  # Schätzwert an der Station

        # Calculate difference between drift corrected reading and the estimates gravity raeding at this point (abw):
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
    df_stat_oesgn = stat_df[stat_df['is_oesgn'] == True]

    # Differences between drift-corrected readings and OESGN g-values:
    df_stat_oesgn['go_mugal'] = df_stat_oesgn['g_oesgn_mugal'] - df_stat_oesgn['g_est_mugal']
    df_stat_oesgn['sig_go_mugal'] = np.sqrt(df_stat_oesgn['sig_g_est_mugal'] ** 2
                                            + df_stat_oesgn['g_sig_oesgn_mugal'] ** 2
                                            + sig_reading ** 2)

    # Calculation of the weighted mean of all go:
    df_stat_oesgn['go_p'] = 1 / (df_stat_oesgn['sig_go_mugal'] ** 2)
    xm_mugal = sum(df_stat_oesgn['go_p'] * df_stat_oesgn['go_mugal']) / sum(
        df_stat_oesgn['go_p'])  # Weighted mean of g0 values (weights: go_p)
    if df_stat_oesgn.shape[0] == 1:  # Only one OEGSN stations measured
        sig_xm_mugal = df_stat_oesgn['g_sig_oesgn_mugal']
    else:
        n = df_stat_oesgn.shape[0]  # Number of observed ÖSGN stations
        sig_xm_mugal = np.sqrt(((sum(df_stat_oesgn['go_mugal'] ** 2) - ((sum(df_stat_oesgn['go_mugal'])) ** 2 / n)) /
                                (n * (n - 1))) + (1 / (sum(df_stat_oesgn['go_p']))))

    # Calculation of the absolute g values at all observed stations:
    stat_df['g_abs_mugal'] = stat_df['g_est_mugal'] + xm_mugal
    stat_df['sig_g_abs_mugal'] = np.sqrt(sig_xm_mugal ** 2 + stat_df['sig_g_est_mugal'] ** 2)
    stat_df['verb_mugal'] = stat_df['g_abs_mugal'] - stat_df['g_oesgn_mugal']
    stat_df['g_abs_full_mugal'] = stat_df['g_abs_mugal'] + 9.8e8  # Get absolute gravity value at station [µGal]

    return stat_df


def write_nsb_file(stat_df, instrument_id, path_name_nsb_file):
    """
    Schreiben der Input Datei für die NSDB (<session_name>.nsb).
    :return:
    """
    fobj_out = open(path_name_nsb_file, 'w')

    # Loop over stations:
    for i, values in stat_df.iterrows():
        out_str = '{:10s} {:8s}  {:9.0f}{:4.0f} {:1s} {:4s}{:5.0f}{:5.0f}\n'\
            .format(values.punktnummer,
                    values.date.strftime('%Y%m%d'),
                    values.g_abs_full_mugal,
                    values.sig_g_abs_mugal,
                    dict_gravimeter_KZ_obs_file[instrument_id],
                    instrument_id,
                    values.dhb_m * 100,
                    values.dhf_m * 100,
                    )
        fobj_out.write(out_str)
    fobj_out.close()


def create_drift_plot(obs_df, stat_df, poly_coef_dict, save_pdf=True, path_name_save_file=''):
    """
    Erstellen eines Drift-Plots.
    :param obs_df: observation dataframe
    :param stat_df: station dataframe
    :param poly_coef_dict: dict with coefficients of drift polanomial
    :param save_pdf: flag, save plot as PDF?
    :param path_name_save_file: path and filename for saved plot (e.g. 'data/n1234')
    :return:
    """
    # Check input arguments
    if save_pdf and (len(path_name_save_file) == 0):
        print('ERROR: Name for output PDF file (drift plot) is not defined!')
        exit()

    # Evaluate drift polynomial:
    pol_degree = len(poly_coef_dict)
    poly_coef = prep_polyval_coef(poly_coef_dict)
    dt_h = np.linspace(obs_df.dt_h.min(), obs_df.dt_h.max(), 100)
    yy_mugal = np.polyval(poly_coef, dt_h)

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
        ax.plot(df_tmp.dt_h, df_tmp.g_red_mugal - g_est, 'o', label=label_str)

    # - Legend and labels:
    plt.legend(loc='best')
    # ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    ax.grid()
    plt.title('Drift-Auswertung: {}'.format(name_obs_file))
    plt.xlabel('Zeit [h]')
    plt.ylabel('Lesung [µGal]')

    # Save as PDF file:
    if save_pdf:
        plt.savefig(path_name_save_file + '_drift.pdf')

    plt.show()


def drift_schwaus(obs_dict, obs_df, stat_df):
    """
    Nachbau alelr Brechnungen auf DRIFT2011.FOR und SCHWAUS2016.FOR.
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
    if verbous:
        print('Schätzung der Drift mit Polynom vom Grad {}:'.format(obs_dict['polynomial_degree']))
    results_dict = calc_drift_corr_mlr(obs_df, stat_df, obs_dict['polynomial_degree'])
    stat_df = results_dict['stat_df']
    if verbous:
        # Print results:
        print(' - Polynom Koeffizeinten (sig = {:5.2f} µGal):'.format(results_dict['pol_coef_sig_mugal']))
        for degree, value in results_dict['pol_coef'].items():
            print('    b{deg} = {value:9.5f} µGal/h^{deg}'.format(deg=degree + 1, value=value))
        print(' - Korrigierte Gravimeter-Lesungen je Station:')
        for i, values in stat_df.iterrows():
            print('    {}: {:12.2f} µGal (sig = {:5.2f} µGal)'.format(values['punktnummer'], values['g_est_mugal'],
                                                                      values['sig_g_est_mugal']))

    # ### Schwere bezogen auf das ÖSGN brechnen ###
    # - Analog zu Fortan Code SCHWAUS2016.FOR
    # - Berechnete absolute Schwere an den beobachteten Stationen bezogen auf das Höheniveau des Festpunktes!
    if verbous:
        print('Berechnung absoluter Schwere-Werte:')
        print(' - Lagerung der Drift-korrigierten Lesungen auf {} ÖSGN Stationen'.format(stat_df[stat_df['is_oesgn']
                                                                                                 == True].shape[0]))
    stat_df = calc_abs_g(stat_df)
    if verbous:
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
    if verbous:
        print('NSDB Input Datei schreiben ({})'.format(path_name_nsb_file))
    write_nsb_file(stat_df, obs_dict['instrument_id'], path_name_nsb_file)

    # ### Create Plot ###
    path_name_drift_plot = out_path + name_obs_file
    if verbous:
        print('Drift-Plot erstellen')
        if flag_save_drift_plot_pdf:
            print(' - Speichern unter: {}_drift.py'.format(path_name_drift_plot))

    create_drift_plot(obs_df, stat_df, results_dict['pol_coef'],
                      save_pdf=flag_save_drift_plot_pdf, path_name_save_file=path_name_drift_plot)


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
            T[i_obs, idx_start*num_est_para_per_instr] = 1  # a_k
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
        print(' - a{}: {:15.3f}'.format(i-num_stations+1, x[i].item()))


def adjust_gauss_markoff(obs_df, stat_df, pol_degree=1):
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
    num_est_para_per_instr = pol_degree+1  # number of estimated parameters per gravimeter
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
    P = np.diag((sig0_mugal / obs_df['g_sd_mugal'])**2)

    # Set up individual matrices:
    for i_obs, values in obs_df.iterrows():
        # print(values['punktnummer'])
        idx_stat = stat_num_dict[values['punktnummer']]  # Get index of observed station
        G[i_obs, idx_stat] = 1
        L[i_obs] = values['g_red_mugal']
        for grav_type, idx_start in gravimeter_type_dict.items():
            for i_deg in range(0, pol_degree+1):
                # print(i_deg)
                D[i_obs, idx_start*num_est_para_per_instr + i_deg] = values['dt_h']**i_deg  # a_k

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

    HP = np.diagflat((sig0_mugal / h_sig_mugal)**2)
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
        print(' - a{}: {:15.3f} ({:6.3f})'.format(i-num_stations, xd[i].item(), sig_xd[i].item()))

    pass


def adjust_gauss_markoff_modelled_obs(obs_df, stat_df, pol_degree=1):
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
    num_est_para_per_instr = pol_degree+1  # number of estimated parameters per gravimeter
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
''
    # Weight matrix:
    P = np.diag((sig0_mugal / obs_df['g_sd_mugal'])**2)

    # Set up individual matrices:
    for i_obs, values in obs_df.iterrows():
        # print(values['punktnummer'])
        idx_stat = stat_num_dict[values['punktnummer']]  # Get index of observed station
        G[i_obs, idx_stat] = 1
        L[i_obs] = values['g_red_mugal']
        for grav_type, idx_start in gravimeter_type_dict.items():
            for i_deg in range(0, pol_degree+1):
                # print(i_deg)
                D[i_obs, idx_start*num_est_para_per_instr + i_deg] = values['dt_h']**i_deg  # a_k

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

    HP = np.diagflat((sig0_mugal / h_sig_mugal)**2)
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
        print(' - a{}: {:15.3f} ({:6.3f})'.format(i-num_stations, Xd[i].item(), sig_xd[i].item()))

    pass

def adjust_vermittelnde_beob_mit_bedingungen(obs_df, stat_df, pol_degree=1):
    """
    Adjustment of gravity meter observations using a vermittelnde Beobachtungen mit Bedingungsgleichungen
    See: AG1, pp.16-17
    Linear problem => Without linearization - no a priori g values required at non-datum stations.
    Estimated parameters:
     - offset and linear drift per instrument
     - gravity at observed stations
    :return:
    """
    pass

# ### Anmerkung ###
# - Ld (ausgeglichene Beobachtungen), mit Drift-Parameter auf Epoche T=0 gerechnet, entspricht genau den
#   Drift-Ergebnissen der MLR mit Pgm DRIFT2011!!! EXAKT (wenn SD der Beob. als ident angenommen!)!
#    - D.h.: Unteschiede in der berechneten absoluten ergeben sich durch die Lagerungs-Bedingung (und deren Gewichtung)!
#  - Wenn nur ein Punkt zur Lagerung und SD der Beob. ident (gleiches Gewicht) => Drift gleich wie bei DRIFT2011!

# ### To Do: ###
#  - Vermittelnder Ausgleich mit Bedingungsgleichungen
#    - Idente Erfebnisse wie mit fiktiven Beobachtungen?
#  - "Datumsfestlgung" nach Thaller
#    - Idente Ergebniss zu Vermittelnder Ausgeleich mit Bedingungsgleichungen?
#  - Bedingungsgleichung definieren, um Ergebnis von SCHWAUS2016.FOR nachzubilden!
#    - Zweck: Genauigkeitsinformation der Stationskooridnaten nutzen!
#  - Wie a priori SD der Datumsstationen nutzen?
#    - Gewichtung bei der Lagerung
#    - Über "Fehlerfortpflanzung" Genauigkeit a posteriori bestimmen!
#  - Drift-Plot erzeugen
#  - Korrelationen visualisieren (Plots erstellen) und analysieren!
#  - Recherche: Kondition einer Matrix, etc. ...


# #########################
# ##### Main function #####
# #########################
def main(path_oesgn_table, name_oesgn_table, path_obs_file, name_obs_file, out_path):

    # ##### Load data #####
    # ÖSGN Tabelle laden:
    if verbous:
        print('Datei einlesen: {}'.format(path_oesgn_table + name_oesgn_table))
    df_oesgn = read_oesgn_table(path_oesgn_table + name_oesgn_table)
    if verbous:
        print(' - {} OESGN Stationen geladen.'.format(df_oesgn.shape[0]))
        # print(df_oesgn.info())
    # Messdatei laden:
    if verbous:
        print('Datei einlesen: {}'.format(path_obs_file + name_obs_file))
    obs_dict = read_obs_file(path_obs_file + name_obs_file)
    if verbous:
        print(' - {} Beobachtungen geladen.'.format(obs_dict['obs_df'].shape[0]))
        # print(obs_dict['obs_df'].info())
    # dataframe für Berechnung erstellen:
    obs_df = obs_dict['obs_df']

    # ### Korrektur der Höhenoffsets vom Instrument zum Boden bzw. Festpunkt (dHF, dHB): ###
    obs_df = corr_ref_heights_instrument(obs_df)

    # ### Dataframe mit Stationsliste erzeugen: ###
    stat_df = get_station_df(obs_df, df_oesgn)
    if verbous:
        print('Station info:')
        print(' - Drift Stationen:   {}'.format(stat_df[stat_df['is_drift_point'] == True].shape[0]))
        print(' - ÖSGN Stationen:    {}'.format(stat_df[stat_df['is_oesgn'] == True].shape[0]))
        print(' - Stationen gesamt:  {}'.format(stat_df.shape[0]))
        print(' - Mittlere Breite:   {:6.3f}°'.format(stat_df.breite_deg.mean()))
        print(' - Mittlere Länge:    {:6.3f}°'.format(stat_df.laenge_deg.mean()))

    # ### Dataframe mit Messungen aufbereiten: ###
    obs_df = prep_obs_df(obs_df, df_oesgn)
    # if verbous:
    # print(obs_df.info())

    # ### Messwert mittels Vertikalgadienten auf die Höhe des Festpunktes reduzieren (mittls dHF) ###
    obs_df = red_g_obs_vg(obs_df)
    if verbous:
        print('Gravimeter-Lesungen mittels VG auf Festpunkt-Niveau reduziert.')

    # ### DRIFT2011 und SCHWAUS2016 ###
    # drift_schwaus(obs_dict, obs_df, stat_df)

    # Adjustment
    # adjust_reilly1970(obs_df, stat_df)

    # Gauss-Markoff model
    adjust_gauss_markoff(obs_df, stat_df, 1)

    # Gauss - Markoff model with computed observatons and linearization:
    adjust_gauss_markoff_modelled_obs(obs_df, stat_df, 1)

    # Vermittelnde Beobachtungen mit Bedgingungsgleichungen:
    # adjust_vermittelnde_beob_mit_bedingungen(obs_df, stat_df, 1)

    print('...fertig!')


# Start main()
if __name__ == "__main__":

    path_oesgn_table = '../data/'
    path_obs_file = '../data/'
    name_oesgn_table = 'OESGN.TAB'
    # name_obs_file = '20200527_tideCorr'
    # name_obs_file = '20200527_sd'
    # name_obs_file = '20200527_2'
    # name_obs_file = '20200527'
    name_obs_file = 'n20200701_1'
    out_path = ''
    # name_obs_file = 'n191021_2'

    if len(sys.argv) == 1:
        print('Eingangsparameter aus dem Python-File bezogen.')
    elif len(sys.argv) == 2:  # Name of obervation file as input argument
        name_obs_file = sys.argv[1]
    else:
        print('Error: Invalid number of input arguments!')
        exit()

    # Start Calculations:
    main(path_oesgn_table, name_oesgn_table, path_obs_file, name_obs_file, out_path)

else:
    # not run as standalone program, but as module
    pass
