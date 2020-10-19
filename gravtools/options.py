"""
User-defined options for gravtools.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

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
verbose = True
flag_save_drift_plot_pdf = True
flag_create_schwaus_protocol = True

# Default names and fielpaths of input files:
path_oesgn_table = '/home/heller/pyProjects/gravtools/data/'
path_obs_file = '/home/heller/pyProjects/gravtools/data/'
name_oesgn_table = 'OESGN.TAB'
# name_obs_file = '20200527_tideCorr'
# name_obs_file = '20200527_sd'
# name_obs_file = '20200527_2'
# name_obs_file = '20200527'
name_obs_file = 'n20200701_1'
name_obs_file = 'e201001'
out_path = '/home/heller/pyProjects/gravtools/out'
