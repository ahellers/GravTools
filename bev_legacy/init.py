"""
Reading data from files and preparing the data for further processing is handled here.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

import sys

import pandas as pd
import numpy as np

# Imports from other gravtools modules:
from bev_legacy import const
from bev_legacy import settings
from bev_legacy import utils


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
        7,  # g [µGal] 58
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
        7,  # dhb [cm]
        6,  # dhf [cm]
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

    df['gravimeter_typ'] = settings.GRAVIMETER_ID_BEV[instrument_id]

    # Store results in dict:
    obs_data['skalierung'] = skalierung
    obs_data['instrument_id'] = instrument_id
    obs_data['polynomial_degree'] = polynomial_degree
    obs_data['date_str'] = date_str
    obs_data['timezone_str'] = timezone_str
    obs_data['obs_df'] = df
    obs_data['filename'] = filename

    return obs_data


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
    obs_df.loc[obs_df.vg_mugal.isna(), 'vg_mugal'] = const.VG_DEFAULT

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
    stat_df['g0_mugal'] = np.nan  # Difference between g_ÖSGN and estimated g_rel
    stat_df['sig_g0_mugal'] = np.nan  # Std.Dev. of difference between g_ÖSGN and estimated g_rel

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



def read_and_prep_data(path_oesgn_table, name_oesgn_table, path_obs_file, name_obs_file):
    # ##### Load data #####
    # ÖSGN Tabelle laden:
    if settings.VERBOSE:
        print('Datei einlesen: {}'.format(path_oesgn_table + name_oesgn_table))
    df_oesgn = read_oesgn_table(path_oesgn_table + name_oesgn_table)
    if settings.VERBOSE:
        print(' - {} OESGN Stationen geladen.'.format(df_oesgn.shape[0]))
        # print(df_oesgn.info())
    # Messdatei laden:
    if settings.VERBOSE:
        print('Datei einlesen: {}'.format(path_obs_file + name_obs_file))
    obs_dict = read_obs_file(name_obs_file, path_obs_file)
    # dataframe für Berechnung erstellen:
    obs_df = obs_dict['obs_df']
    del obs_dict['obs_df']
    obs_info_dict = obs_dict
    del obs_dict

    if settings.VERBOSE:
        print(' - {} Beobachtungen geladen.'.format(obs_df.shape[0]))
        # print(obs_dict['obs_df'].info())

    # ### Korrektur der Höhenoffsets vom Instrument zum Boden bzw. Festpunkt (dHF, dHB): ###
    obs_df = utils.corr_instrument_ref_heights(obs_df)

    # ### Dataframe mit Stationsliste erzeugen: ###
    stat_df = get_station_df(obs_df, df_oesgn)
    if settings.VERBOSE:
        print('Stations info:')
        print(' - Drift Stationen:   {}'.format(stat_df[stat_df['is_drift_point'] == True].shape[0]))
        print(' - ÖSGN Stationen:    {}'.format(stat_df[stat_df['is_oesgn'] == True].shape[0]))
        print(' - Stationen gesamt:  {}'.format(stat_df.shape[0]))
        print(' - Mittlere Breite:   {:6.3f}°'.format(stat_df.breite_deg.mean()))
        print(' - Mittlere Länge:    {:6.3f}°'.format(stat_df.laenge_deg.mean()))

    # ### Dataframe mit Messungen aufbereiten: ###
    obs_df = prep_obs_df(obs_df, df_oesgn)

    # ### Messwert mittels Vertikalgadienten auf die Höhe des Festpunktes reduzieren (mittls dHF) ###
    obs_df = utils.red_g_obs_vg(obs_df)
    if settings.VERBOSE:
        print('Gravimeter-Lesungen mittels VG auf Festpunkt-Niveau reduziert.')

    return obs_df, stat_df, df_oesgn, obs_info_dict


# Run as standalone program:
if __name__ == "__main__":

    if len(sys.argv) == 1:
        print('Eingangsparameter aus options.py bezogen (Default-Parameter).')
        path_oesgn_table = settings.PATH_OESGN_TABLE
        name_oesgn_table = settings.NAME_OESGN_TABLE
        path_obs_file = settings.PATH_OBS_FILE_BEV
        name_obs_file = settings.NAME_OBS_FILE_BEV

    elif len(sys.argv) == 2:  # Name of obervation file as input argument
        name_obs_file = sys.argv[1]
        # Default-parameters from options-file:
        path_oesgn_table = settings.PATH_OESGN_TABLE
        name_oesgn_table = settings.NAME_OESGN_TABLE
        path_obs_file = settings.PATH_OBS_FILE_BEV

    else:
        print('Error: Invalid number of input arguments!')
        exit()

    # Run init module:
    obs_df, stat_df, oesgn_df, obs_info_dict = read_and_prep_data(path_oesgn_table,
                                                                  name_oesgn_table,
                                                                  path_obs_file,
                                                                  name_obs_file)
    # Print dataframe info:
    print('\n--- Observation dataframe ---')
    obs_df.info()
    print('\n--- Stations dataframe ---')
    stat_df.info()
    print('\n--- ÖSGN dataframe ---')
    oesgn_df.info()

else:
    # not run as standalone program, but as module
    pass
