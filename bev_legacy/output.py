"""
Output functions of gravtools.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

import datetime as dt

# Imports from other gravtools modules:
from bev_legacy import settings


def write_nsb_file(stat_df, instrument_id, path_name_nsb_file):
    """
    Schreiben der Input Datei für die NSDB (<session_name>.nsb).

    :param stat_df: Station dataframe
    :param path_name_nsb_file: path and name of output file
    :param instrument_id: instrument ID of the used gravimeter
    :return:
    """
    fobj_out = open(path_name_nsb_file, 'w')

    # Loop over stations:
    for i, values in stat_df.iterrows():
        out_str = '{:10s} {:8s}  {:9.0f}{:4.0f} {:1s} {:4s}{:5.0f}{:5.0f}\n'.format(
            values.punktnummer,
            values.date.strftime('%Y%m%d'),
            values.g_abs_full_mugal,
            values.sig_g_abs_mugal,
            settings.GRAVIMETER_KZ_BEV[instrument_id],
            instrument_id,
            values.dhb_m * 100,
            values.dhf_m * 100,
        )
        fobj_out.write(out_str)
    fobj_out.close()


def write_schwaus_protcol(obs_df, stat_df, pol_coef, pol_coef_sig_mugal, session_name, path_save_file,
                          instrument_id, skalierung=None, g0_mean_mugal=None, sig_g0_mean_mugal=None):
    """
    Create processing and results protocol for the analyzed gravity observations.

    :param obs_df:
    :param stat_df:
    :param pol_coef:
    :param pol_coef_sig_mugal:
    :param session_name:
    :param path_save_file:
    :param skalierung:
    :return:
    """

    file_path_name_str =  path_save_file + session_name + '_prot.txt'
    td = (obs_df.index.max() - obs_df.index.min())  # Session duration
    num_of_obs = obs_df.shape[0]
    gravimeter_str = settings.GRAVIMETER_ID_BEV[instrument_id]
    grav_height_corr_cm = settings.GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[gravimeter_str] * 1e2  # in [cm]

    with open(file_path_name_str, 'w') as f:
        f.write('### Protokoll zur Gravimetrie Auswertung ({}) ###\n'.format(dt.datetime.now().strftime('%Y-%m-%d, '
                                                                                                      '%H:%M')))
        f.write('\n')
        f.write(' - Messdatei:          {}\n'.format(session_name))
        f.write('   - Start der Messung:  {} [{}]\n'.format(obs_df.index.min().strftime('%Y-%m-%d, %H:%M'),
                                                            obs_df.index.min().tzname()))
        f.write('   - Ende der Messung:   {} [{}]\n'.format(obs_df.index.max().strftime('%Y-%m-%d, %H:%M'),
                                                            obs_df.index.min().tzname()))
        f.write('   - Messdauer:          {} Tage, {} h, {} min\n'.format(td.days, td.seconds//3600,
                                                                       (td.seconds//60) % 60))
        f.write('   - Anz. Beobachtungen: {} ({:3.1f}/h)\n'.format(num_of_obs, num_of_obs / (td.seconds/3600)))
        f.write('   - Skalierung:         {}\n'.format(skalierung))
        f.write('   - Annahme: Beobachtungen bereits Gezeiten-korrigiert!\n')
        f.write('   - Beobachtungen mittels VG (und dhf, dhb) auf Festpunkt-Niveau reduziert.')
        f.write(' - Gravimeter:         {}\n'.format(settings.GRAVIMETER_ID_BEV[instrument_id]))
        f.write('   - Bezugshöhen (dhf, dhf) korrigiert um: {:+5.1f} cm\n'.format(grav_height_corr_cm))
        f.write(' - Stationen gesamt:   {}\n'.format(stat_df.shape[0]))
        f.write('   - Drift Stationen:    {}\n'.format(stat_df[stat_df['is_drift_point'] == True].shape[0]))
        f.write('   - ÖSGN Stationen:     {}\n'.format(stat_df[stat_df['is_oesgn'] == True].shape[0]))
        f.write('   - Mittlere Breite:    {:6.3f}°\n'.format(stat_df.breite_deg.mean()))
        f.write('   - Mittlere Länge:     {:6.3f}°\n'.format(stat_df.laenge_deg.mean()))
        f.write('\n')
        f.write('\n')
        f.write('##### Drift-Korrektur #####\n')
        f.write(' - Anpassung eines Drift-Polynoms (Grad {}) mittels multipler linearer Regression.\n'.format(
            len(pol_coef)))
        f.write('   - Polynomkoeffizienten:\n')
        for degree, coef in pol_coef.items():
            f.write('     - A{} = {:8.6f} µGal/h ({:8.6f} mGal/Tag)\n'.format(degree, coef, coef*24/1e3))
        # Mittere quadratische Abweichung der Drift-korrigierten Lesungen vom Lesungs-Schätzwert an der Station
        f.write('     - Mittlere quadr. Abweichung = {:8.6f} µGal\n'.format(pol_coef_sig_mugal))

        f.write('\n')
        f.write('### (Reduzierte) Lesungen und Schätzwerte an den Stationen ###\n')
        f.write('Zeit   Lesung   dhf   VG       Red. Lesung   Abw.\n')
        f.write('[{}]  [µGal]   [cm]  [µGal/m] [µGal]        [µGal]\n'.format(obs_df.index.tz.zone))
        f.write('\n')
        for i, stats in stat_df.iterrows():
            stat_name = stats.punktnummer

            for j, obs in obs_df.loc[obs_df['punktnummer'] == stat_name].iterrows():
                f.write('{:2d}:{:2d}  {:7.0f}  {:4.1f}  {:5.1f}    {:9.1f}   {:5.1f}\n'.format(
                    obs.name.hour,
                    obs.name.minute,
                    obs.g_mugal,
                    obs.dhf_m * 100,
                    obs.vg_mugal,
                    obs.g_red_mugal,
                    obs.abw_mugal))
                pass

            f.write(' => {:10s} : {:9.1f} µGal (+/- {:5.1f} µGal)\n'.format(stat_name,
                                                                            stats.g_est_mugal,
                                                                            stats.sig_g_est_mugal))
            if stats.is_oesgn:
                f.write('    Koordinaten ÖSGN: {:8.4f}°, {:8.4f}° \n'.format(stats.breite_deg, stats.laenge_deg))
            f.write('\n')

        f.write('\n')
        f.write('##### Absolute Schwere mit Bezug zum ÖSGN #####\n')
        f.write(' - Mittelwert der Abweichungen an ÖSGN-Stationen:\n')
        f.write('     go_mean = {:9.1f} +/- {:5.1f} µGal\n'.format(g0_mean_mugal, sig_g0_mean_mugal))
        f.write('\n')
        f.write('PunktNr.    g_ÖSGN       g_relativ         g0                 g_abs            Verb.\n'.format(obs_df.index.tz.zone))
        f.write('            [µGal]       [µGal]            [µGal]             [µGal]           [µGal]\n')
        # f.write('\n')
        for i, stats in stat_df.iterrows():
            if stats.is_oesgn:
                f.write('{:10s}  {:6.0f}+-{:3.0f}  {:9.1f}+-{:5.1f}  {:9.1f}+-{:5.1f}  {:8.1f}+-{:5.1f}  ' \
                        '{:5.1f}\n'.format(stats.punktnummer,
                                              stats.g_oesgn_mugal,
                                              stats.g_sig_oesgn_mugal,
                                              stats.g_est_mugal,
                                              stats.sig_g_est_mugal,
                                              stats.g0_mugal,
                                              stats.sig_g0_mugal,
                                              stats.g_abs_mugal,
                                              stats.sig_g_abs_mugal,
                                              stats.verb_mugal))
            else:
                f.write('{:10s}               {:9.1f}+-{:5.1f}                     {:8.1f}+-{:5.1f}\n'.format(
                                           stats.punktnummer,
                                           stats.g_est_mugal,
                                           stats.sig_g_est_mugal,
                                           stats.g_abs_mugal,
                                           stats.sig_g_abs_mugal,))

    pass

