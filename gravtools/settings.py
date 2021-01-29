"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
User-defined settings for gravtools.
"""

SURVEY_DATA_SOURCE_TYPES = {
    'cg5_obs_file_txt': 'Scintrex CG5 observation file (text format)',
    'bev_obs_file': 'Simple observation file format used by BEV',
}

STATION_DATA_SOURCE_TYPES = {
    'oesgn_table': 'Control points of the Austrian gravity base network (OESGN).',
    'obs_file': 'From an observation file'
}

TIDE_CORRECTION_TYPES = {
    'cg5_longman1959': 'Instrument-implemented tidal correction of the Scintrex CG-5',
    'no_tide_corr': 'No tide correction applied',
    'unknown': 'Unknown whether a tide correction was applied',
}

REFERENCE_HEIGHT_TYPE = {
    'sensor_height': 'The gravity reverence point at the station is at the sensor height',
    'instrument_top': 'The gravity value refers to the height og the instrument top',
    'ground': 'The gravity value refers to the ground point at the station',
    'control_point': 'The gravity value refers to the control point at the station',
}

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

#  Lookup table for matching gravimeter IDs and the tidal corrections that are applied per default in the BEV legacy
#  observation files:
BEV_GRAVIMETER_TIDE_CORR_LOOKUP = {
    '5': 'cg5_longman1959'
}

DEFAULT_GRAVIMETER_ID_CG5_SURVEY = '5'

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
# NAME_OBS_FILE_CG5 = '20200907_test.TXT'
