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

# SCHWAUS and DRIFT settings:
FLAG_SAVE_DRIFT_PLOT_PDF = True
FLAG_CREATE_SCHWAUS_PROTOCOL = True
VERBOSE = True

# Output directory
OUT_PATH = '/home/heller/pyProjects/gravtools/out'

# Default names and paths of input files:

# ÖSGN Table:
PATH_OESGN_TABLE = '/home/heller/pyProjects/gravtools/data/'
NAME_OESGN_TABLE = 'OESGN.TAB'

# BEV observation files:
PATH_OBS_FILE_BEV = '/home/heller/pyProjects/gravtools/data/BEV/'
# NAME_OBS_FILE_BEV = '20200527_tideCorr'
# NAME_OBS_FILE_BEV = '20200527_sd'
# NAME_OBS_FILE_BEV = '20200527_2'
# NAME_OBS_FILE_BEV = '20200527'
NAME_OBS_FILE_BEV = 'n20200701_1'
# NAME_OBS_FILE_BEV = 'e201001'

# CG-5 observation files (text)
PATH_OBS_FILE_CG5 = '/home/heller/pyProjects/gravtools/data/CG5/'
NAME_OBS_FILE_CG5 = '2020-06-18_DACH.TXT'
#NAME_OBS_FILE_CG5 = '20200907_test.TXT'
