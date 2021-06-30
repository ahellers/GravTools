"""
Utility functions for gravtools.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

# Imports from other gravtools modules:
from bev_legacy import settings


def corr_instrument_ref_heights(obs_df):
    """
    Korrektur der Höhenoffsets zwischen Instrument-Oberkante (bei Messung notiert) und Boden bzw. Festpunkt (dHF, dHB)
    um den Höhenoffset zwischen Instrumenten-Oberkante und Messwerk zu berücksichtigen. Diese Korrektur ist abhängig vom
    Intrumenten-Typ.
    :param obs_df: dataframe mit Messungen
    :return: dataframe mit Messungen (erweitert)
    """

    # ### Je nach Gerätetyp (CG5, CG3) die Bezugs-Höhen (dhb, dhf) anpassen ###:
    for grav_typ, corr_m in settings.GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m.items():
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