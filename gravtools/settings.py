"""
User-defined settings for gravtools.
This module is part of the gravtools package.
Author: Andreas Hellerschmied
"""

GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m = {
    'CG3': -0.211,
    'CG5': -0.211,
}

# Instrumenten-IDs in der Messdatei einem Instrument zuweisen:
GRAVIMETER_ID_BEV = {
    '3': 'CG3',
    '5': 'CG5',
}

# Instrumenten-IDs in der Messdatei einem Instrumenten-KZ zuweisen:
GRAVIMETER_KZ_BEV = {
    '3': 'C',
    '5': 'C',
    '9': 'D',
    '90': 'D',
    '900': 'D',
    '51': 'D',
    '510': 'D',
    '500': 'W',
}

VERBOSE = True
FLAG_SAVE_DRIFT_PLOT_PDF = True
FLAG_CREATE_SCHWAUS_PROTOCOL = True

# Default names and fielpaths of input files:
PATH_OESGN_TABLE = '/home/heller/pyProjects/gravtools/data/'
PATH_OBS_FILE_BEV = '/home/heller/pyProjects/gravtools/data/'
NAME_OESGN_TABLE = 'OESGN.TAB'
# NAME_OBS_FILE_BEV = '20200527_tideCorr'
# NAME_OBS_FILE_BEV = '20200527_sd'
# NAME_OBS_FILE_BEV = '20200527_2'
# NAME_OBS_FILE_BEV = '20200527'
NAME_OBS_FILE_BEV = 'n20200701_1'
NAME_OBS_FILE_BEV = 'e201001'
OUT_PATH = '/home/heller/pyProjects/gravtools/out'
