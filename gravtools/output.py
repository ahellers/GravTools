"""
Output functions of gravtools.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

# Imports from other gravtools modules:
from gravtools import options


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
            options.dict_gravimeter_KZ_obs_file[instrument_id],
            instrument_id,
            values.dhb_m * 100,
            values.dhf_m * 100,
        )
        fobj_out.write(out_str)
    fobj_out.close()
