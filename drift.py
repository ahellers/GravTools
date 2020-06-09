"""
Programm zur Schätzung der Gerätedrift bei Relativgravimetern
Autor: Andreas Hellerschmied
Date: 2020-05-08
Benützung:

"""

# Imports:
import pandas as pd
import numpy as np
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

vg_default = 308.6  # muGal/m

verbous = True


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
    df = pd.read_fwf(filepath + filename, widths=widths, header=None, names=column_names, skiprows=5, skipfooter=1,
                     dtype={'zeit': object})
    df['zeit'] = pd.to_datetime(date_str + ' ' + df['zeit'] + ' ' + timezone_str, format='%Y %m %d %H.%M %Z')
    df = df.set_index('zeit')  # set 'zeit' as index
    df['g_mugal'] = df['g_mugal'] * 1e3  # Conversion from mGal to muGal
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
        obs_df.loc[obs_df['gravimeter_typ'] == grav_typ, 'dhb_m'] = obs_df.loc[obs_df[
                                                                                   'gravimeter_typ'] == grav_typ, 'dhb_m'] + corr_m
        obs_df.loc[obs_df['gravimeter_typ'] == grav_typ, 'dhf_m'] = obs_df.loc[obs_df[
                                                                                   'gravimeter_typ'] == grav_typ, 'dhf_m'] + corr_m
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
    # stat_df['dhb_m'] = np.nan  # Depends on the setup of each measurement => May not be const. for a station!
    # stat_df['dhf_m'] = np.nan  # Depends on the setup of each measurement => May not be const. for a station!
    stat_df['verb_mugal'] = np.nan

    # Add ÖSGN data:
    stat_df = stat_df.merge(df_oesgn.drop(columns=['anmerkungen', 'identitaet', 'date']), on='punktnummer', how='left')
    stat_df['is_oesgn'] = ~stat_df.g_oesgn_mugal.isna()

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
                np.sqrt((sum(obs_df.loc[obs_df['punktnummer'] == stat, 'abw_mugal']**2) / (num_of_obs-1)) + 4)

    pol_coef_sig_mugal = np.sqrt(
        sum(obs_df.abw_mugal ** 2) / (obs_df.shape[0] - 1 - sum(stat_df.is_drift_point == False)))

    # Sigma for all estimats with only one observation:
    stat_df.loc[~stat_df['is_drift_point'], 'sig_g_est_mugal'] = np.sqrt(pol_coef_sig_mugal**2 + 25)

    results_dict = dict()
    results_dict['pol_coef'] = pol_coef_dict
    results_dict['pol_coef_sig_mugal'] = pol_coef_sig_mugal
    results_dict['stat_df'] = stat_df
    results_dict['obs_df'] = obs_df  # Erweitert mit Verbessungen und Drift-Korekturen je Beobachtung

    return results_dict


# #########################
# ##### Main function #####
# #########################
def main():
    # Options:
    path_oesgn_table = './data/'
    path_obs_file = './data/'
    name_oesgn_table = 'OESGN.TAB'
    # name_obs_file = '20200527'
    name_obs_file = 'n191021_2'

    # ##### Load data #####
    # ÖSGN Tabelle laden:
    if verbous:
        print('Read file: {}'.format(path_oesgn_table + name_oesgn_table))
    df_oesgn = read_oesgn_table(path_oesgn_table + name_oesgn_table)
    if verbous:
        print(' - {} OESGN Stationen geladen.'.format(df_oesgn.shape[0]))
        # print(df_oesgn.info())
    # Messdatei laden:
    if verbous:
        print('Read file: {}'.format(path_obs_file + name_obs_file))
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
        print(' - ÖSGN Stationen:    {}'.format(stat_df[stat_df['is_drift_point'] == True].shape[0]))
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

    # ### Drift-Korrektur ###
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
        print(' - Korrigierte Gravimeter-Lesungen:')
        tmp_df = results_dict['stat_df']
        for i, values in stat_df.iterrows():
            print('    {}: {:12.2f} µGal (sig = {:5.2f} µGal)'.format(values['punktnummer'], values['g_est_mugal'],
                                                            values['sig_g_est_mugal']))

    # ### Create Plot ###

    # Evaluate drift polynomial:
    poly_coef = prep_polyval_coef(results_dict['pol_coef'])
    dt_h = np.linspace(obs_df.dt_h.min(), obs_df.dt_h.max(), 100)
    yy_mugal = np.polyval(poly_coef, dt_h)

    fig, ax = plt.subplots()
    ax.plot(dt_h, yy_mugal, '--', label='Polynom (n={})'.format(obs_dict['polynomial_degree']))
    for pkt_num in stat_df.punktnummer:
        # print(pkt_num)
        df_tmp = obs_df.loc[obs_df['punktnummer'] == pkt_num]
        # print(df_tmp)
        g_est = stat_df.loc[stat_df['punktnummer'] == pkt_num].g_est_mugal.values[0]
        ax.plot(df_tmp.dt_h, df_tmp.g_red_mugal - g_est, 'o', label=pkt_num)

    # - Legend and labels:
    plt.legend(loc='best')
    # ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    ax.grid()
    plt.title('Drift-Auswertung: {}'.format(name_obs_file))
    plt.xlabel('Zeit [h]')
    plt.ylabel('Lesung [µGal]')
    plt.show()

    print('...finished!')


# Start main()
if __name__ == "__main__":
    main()

    # Test:

else:
    # not run as standalone program, but as module
    pass
