"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
User-defined settings for gravtools.
"""

# Additive constant for the determination of the full absolute gravity [µGal] from observed values:
ADDITIVE_CONST_ABS_GRTAVITY = 9.8e8

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

# Reference heights (sensor heights) of different gravimeter types:
GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m = {
    'CG3': -0.211,
    'CG5': -0.211,
}

# Valid gravimeter types and description.
GRAVIMETER_TYPES = {
    'CG5': 'Sctintrex CG5',
    'CG3': 'Sctintrex CG3',
}

# Default Gavimeter Type when loading data from an CG5 observation file
DEFAULT_GRAVIMETER_TYPE_CG5_SURVEY = 'CG5'

# Lookuptable to convert gravimeter type to Kennzeichen used at BEV (in the database NSDB):
GRAVIMETER_TYPES_KZG_LOOKUPTABLE = {
    'CG5': 'C',
    'CG3': 'C',
}

# Valid Gravimeter serial numbers (S/N and type)
GRAVIMETER_SERIAL_NUMBERS = {
    '40601': 'CG5'
}

# Lookuptable to convert the gravimeter S/N to the IDs written to the database NSDB:
GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE = {
    '40601': '5'
}


#--------------

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

#--------------------

#  Lookup table for matching gravimeter IDs and the tidal corrections that are applied per default in the BEV legacy
#  observation files:
BEV_GRAVIMETER_TIDE_CORR_LOOKUP = {
    '5': 'cg5_longman1959'
}

# Available adjustment methods
ADJUSTMENT_METHODS = {
    'LSM_diff': 'LSM (differential observations)',  # Least-squares adjustment of differential observations
    'MLR_BEV': 'MLR (BEV legacy processing)',  # Least-squares adjustment of differential observations
}

# Treshold for the "Gewichtsreziprokenprobe nach Ansermet" (see Skriptum AG1, p. 136, Eq. (6.86))
ANSERMET_DIFF_TRESHOLD = 1e-3

# Treshold for the redundancy component of an observation in order to apply a pope test for outlier detection:
R_POPE_TEST_TRESHOLD = 1e-6

# Only consider active observations for the determination of the reference epochs, e.g. for the drift polynomial. The
# Reference epochs are determined based on the first (active only or active/inactive) observations in the campaign or
# in each individual survey, depending on the settings.
ACTIVE_OBS_ONLY_FOR_REF_EPOCH = True

# GUI and Program options:
CALCULATE_REDUCED_OBS_WHEN_LOADING_DATA = True  # Calculate reduced observations when loading observation data.

# SCHWAUS and DRIFT settings:
FLAG_SAVE_DRIFT_PLOT_PDF = True
FLAG_CREATE_SCHWAUS_PROTOCOL = True
VERBOSE = True

# Output directory
OUT_PATH = '/home/heller/pyProjects/GravTools/out/'

# Default names and paths of input files:

# ÖSGN Table:
PATH_OESGN_TABLE = '/home/heller/pyProjects/GravTools/data/'
NAME_OESGN_TABLE = 'OESGN_f200701.tab'

# BEV observation files:
PATH_OBS_FILE_BEV = '/home/heller/pyProjects/GravTools/data/BEV/'
# NAME_OBS_FILE_BEV = '20200527_tideCorr'
# NAME_OBS_FILE_BEV = '20200527_sd'
# NAME_OBS_FILE_BEV = '20200527_2'
# NAME_OBS_FILE_BEV = '20200527'
#NAME_OBS_FILE_BEV = 'n20200701_1'
#NAME_OBS_FILE_BEV = 'f200701_1'
NAME_OBS_FILE_BEV = 'f200701_3'
# NAME_OBS_FILE_BEV = 'e201001'

# CG-5 observation files (text)
PATH_OBS_FILE_CG5 = '/home/heller/pyProjects/GravTools/data/CG5/'
# NAME_OBS_FILE_CG5 = '2020-06-18_DACH.TXT'
NAME_OBS_FILE_CG5 = '20200907_test.TXT'
