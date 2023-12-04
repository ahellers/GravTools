"""GravTools settings defined by the user.

Copyright (C) 2021  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

# Additive constant for the determination of the full absolute gravity [µGal] from observed values:
ADDITIVE_CONST_ABS_GRTAVITY = 9.8e8

SURVEY_DATA_SOURCE_TYPES = {
    'cg5_obs_file_txt': 'Scintrex CG5 observation file (text format)',
    'bev_obs_file': 'Simple observation file format used by BEV',
}

STATION_DATA_SOURCE_TYPES = {
    'oesgn_table': 'Austrian gravity base network stations (OESGN).',
    'csv_stat_file': 'CVS station file',
    'obs_file': 'From an observation file',
}

TIDE_CORRECTION_TYPES = {
    'cg5_longman1959': 'Instrument-implemented tidal correction of the Scintrex CG-5',
    'longman1959': 'Tidal corrections by the model of Longman (1959)',
    'no_tide_corr': 'No tide correction applied',
    'from_time_series': 'Interpolated from correction time series',
    'unknown': 'Unknown whether a tide correction was applied',
}

ATM_PRES_CORRECTION_TYPES = {
    'no_atm_pres_corr': 'No atmospheric pressure correction applied',
    'iso_2533_1975': 'Atmopheric pressure correction using ISO 2533:1975 normal air pressure',
}

# Default admittance factor for the calculation of atmospheric pressure correction using the ISO 2533:1975 normal air pressure
ATM_PRES_CORRECTION_ADMITTANCE_DEFAULT = 0.3

# Conversion factors of different units of gravity to µGal:
# - Keys: Units as used e.g. in TSF files by Tsoft
# - Values: Multiplicative conversion factors to µGal
UNIT_CONVERSION_TO_MUGAL = {
    'nm/s^2': 1e-1  # 1 nm/s^2 * 0.1 = 1 µGal
}

# Interpolation methods provided by scipy.interpolation.interp1, e.g. used for the interpolation of correction time
# series data.
SCIPY_INTERP1_INTERPOLATION_METHODS = {
    'linear': 'Linear interpolation.',
    'nearest': '"Snaps" to the nearest data point. Rounds down at half-integers (e.g. 0.5, 1.5)',
    'nearest-up': '"Snaps" to the nearest data point. Rounds up at half-integers (e.g. 0.5, 1.5)',
    'zero': 'Zero order spline. It`s value at any point is the last raw value seen.',
    'slinear': 'First order spline interpolation (similar to "linear").',
    'quadratic': 'Second order spline interpolation.',
    'cubic': 'Third order spline interpolation.',
    'previous': 'Return previous value.',
    'next': 'Return next value.',
}
# Default interpolation method:
SCIPY_INTERP1_INTERPOLATION_DEFAULT_METHOD = 'quadratic'

REFERENCE_HEIGHT_TYPE = {
    'sensor_height': 'The gravity values refer to the sensor height',
    'instrument_top': 'The gravity values refers to the instrument top',
    'ground': 'The gravity values refers to the ground point of the station',
    'control_point': 'The gravity values refers to the control point of the station',
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
    '40601': 'CG5',
    '40236': 'CG5',
}

# Lookuptable to convert the gravimeter S/N to the IDs written to the database NSDB:
GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE = {
    '40601': '5',
    '40236': '5',
}

#  Lookup table for matching gravimeter IDs and the tidal corrections that are applied per default in the BEV legacy
#  observation files:
BEV_GRAVIMETER_TIDE_CORR_LOOKUP = {
    '5': 'cg5_longman1959'
}

# Methods for calculation of setup observations:
SETUP_CALC_METHODS = {
    'variance_weighted_mean': 'Variance weighted mean of observations within a setup.',
    'individual_obs': 'Each observation is treated as setup.'
}

# Methods for calculation of standard deviations (SD) of setup observations:
SETUP_SD_METHODS = {
    'sd_from_obs_file': 'Apply the standard deviations from observation files.',
    'sd_default_per_obs': 'Apply the default standard deviation to all individual observations.',
    'sd_default_per_setup': 'Apply the default standard deviation to setup obsertvations.'
}

# Available adjustment methods:
ADJUSTMENT_METHODS = {
    'LSM_non_diff': 'LSM (non-differential observations)',  # Least-squares adjustment of non-differential observations
    'LSM_diff': 'LSM (differential observations)',  # Least-squares adjustment of differential observations
    'MLR_BEV': 'MLR (BEV legacy processing)',  # Least-squares adjustment of differential observations
    'VG_LSM_nondiff': 'VG LSM (non-differential observations)',  # Vertical Gravity Gradient estimation based on Least-squares adjustment of non-differential observations
}

# List of LSM methods where export of nsb files is allowed:
LSM_METHODS_NSD_FILE_EXPORT = [
    'LSM_diff',
    'LSM_non_diff',
    'MLR_BEV',
]

# List of LSM methods where export of gis files is allowed:
LSM_METHODS_GIS_EXPORT = [
    'LSM_diff',
    'LSM_non_diff',
]

# Default subdirectory in the campaign's output directory for exporting the results as shapefiles:
GIS_RESULTS_OUTPUT_SUBDIR = 'GIS_results'

# List of LSM methods where VG plots are created (and should be saved as PNG files):
LSM_METHODS_VG_PLOT = [
    'VG_LSM_nondiff',
]

# Available iteration approaches for scaling the SD of setup observations:
ITERATION_APPROACHES = {
    'Multiplicative': 'Multiplicative iteration approach',
    'Additive': 'Additive iteration approach'
}

# Treshold for the "Gewichtsreziprokenprobe nach Ansermet" (see Skriptum AG1, p. 136, Eq. (6.86))
ANSERMET_DIFF_TRESHOLD = 1e-3

# Treshold for the redundancy component of an observation in order to apply a pope test for outlier detection:
R_POPE_TEST_TRESHOLD = 1e-6

# GUI and Program options:
CALCULATE_REDUCED_OBS_WHEN_LOADING_DATA = True  # Calculate reduced observations when loading observation data.
INIT_OESGN_STATION_AS_DATUM = False  # Initialize OESGN stations as datum stations when loading rom an OESGN file?


# ----- GUI appearance and plot settings -----

# Methods for the determination of histogram bin edges:
NUMPY_HISTOGRAM_BIN_EDGES_OPTIONS = {
    'auto': 'Maximum of the ‘sturges’ and ‘fd’ estimators. Provides good all around performance.',
    'fd': 'Freedman Diaconis Estimator: Robust (resilient to outliers) estimator that takes into account data variability and data size.',
    'doane': 'An improved version of Sturges’ estimator that works better with non-normal datasets.',
    'scott': 'Estimator based on leave-one-out cross-validation estimate of the integrated squared error. Can be regarded as a generalization of Scott’s rule.',
    'rice': 'Estimator does not take variability into account, only data size. Commonly overestimates number of bins required.',
    'sturges': 'R’s default method, only accounts for data size. Only optimal for gaussian data and underestimates number of bins for large non-gaussian datasets.',
    'sqrt': 'Square root (of data size) estimator, used by Excel and other programs for its speed and simplicity.',
    'Num. of bins': 'User defined number of bins (min. = 2, max. = 1000).'
}

# Time label format for y-ticks in time-series plots:
# - Format string for datetime.strftime()
Y_TICK_DATETIME_FORMAT = '%Y-%m-%d %H:%M:%S'

# --- Gerneral color settings ---
DATUM_STATION_COLOR = (255, 204, 204)

# --- Drift plots in the results tab: ---
# Number of plot items for plotting the drift function (polynomial):
DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION = 100
# Options for the observation data points in the drift plot:
DRIFT_PLOT_SCATTER_PLOT_SYMBOL_SIZE = 10
DRIFT_PLOT_SCATTER_PLOT_PEN_WIDTH = 1
DRIFT_PLOT_SCATTER_PLOT_PEN_COLOR = 'k'

# --- Correlation matrix in the results tab: ---
# Background color scheme for the correlation coefficient [0 - 1]
# - Get HEX colors from here: https://colorbrewer2.org/
# - https://www.pythonguis.com/tutorials/qtableview-modelviews-numpy-pandas/
# green (low correlation) => red (high correlation)
CORRELATION_COEF_COLORS = ['#006837', '#1a9850', '#66bd63', '#a6d96a', '#d9ef8b', '#fee08b', '#fdae61', '#f46d43', '#d73027', '#a50026']
# Background color for diagonal elements (=1):
CORRELATION_COEF_DIAG_ELEMENTS = '#bababa'  # light grey

# --- VG plot in the results tab: ---
# Min. and max. height for the horizontal plot axis [m]:
# - "None" indicates that the min/max of the horizontal axis is adjusted to the actual setup height (dhf) of the
# observations. In this case the plotted height range is widened by "VG_PLOT_HEIGHT_DELTA_M".
VG_PLOT_MIN_HEIGHT_M = 0.0
VG_PLOT_MAX_HEIGHT_M = 1.8
VG_PLOT_HEIGHT_DELTA_M = 0.1
# Number of datapoints plotted in the range betwen min. and max. height:
VG_PLOT_NUM_ITEMS_VG_POLYNOMIAL = 100
# Min. Y-range (default: [VG_PLOT_MIN_LOWER_L_RANGE, VG_PLOT_MIN_UPPER_Y_RANGE] = [-1, +1]):
VG_PLOT_MIN_UPPER_Y_RANGE = 1.0
VG_PLOT_MIN_LOWER_L_RANGE = -1.0
# Define pen for scatter plot of residuals:
VG_PLOT_SCATTER_PLOT_PEN_COLOR = 'k'
VG_PLOT_SCATTER_PLOT_PEN_WIDTH = 1
VG_PLOT_SCATTER_PLOT_SYMBOL_SIZE = 10

# ----- Data export options: -----
# List of columns in the `obs_df` dataframe that are written to the exported observation list CSV file:
EXPORT_OBS_LIST_COLUMNS = ['survey_name', 'station_name', 'obs_epoch', 'keep_obs']
# Maximum allowed SD of the estimated gravity at stations when exporting data to a nsb file!
#  => This is important as only 3 characters are reserved in the nsb file for the SD!
MAX_SD_FOR_EXPORT_TO_NSB_FILE = 999.0
# Choices for 5-character comment written to the nsb file (appears as Operat-G in the NSDB):
# - One letter instrument-ID: 'inst_id'
# - Official 5-character serial number of the Scintrex CG-5 Gravitimeters: 'cg5_serial_number'
# - Version of Gravtools: 'gravtools_version'
# !! Warning: Actually only the gravtools version makes sense, because multiple surveys observed with different
# instruments can be combindes
WRITE_COMMENT_TO_NSB = 'gravtools_version'

# ----- Program and Software settings -----
# Specify the protocol version of pickle used when saving the campaign data to a pickle file.
# Python 3.6.8 (used at BEV) does not support version 5, only up to version 4
# Options (type: int):
# - <int number> => Directly specify the protocol version
#                   (4 should be a good choice, as it is compatible with python 3.6)
# - 999 => Use the highest version available at the installation (pickle.HIGHEST_VERSION)
PICKLE_PROTOCOL_VERSION = 4


# ----- GIS data export settings: -----
DEFAULT_EPSG_CODE = 4312  # 4312: MGI, lat/lon, Greenwich
# Default filenames for shapefile export from the results tab:
DEFUALT_FILENAME_OBERVATION_RESULTS_SHP = 'obs_results_'  # + <LSM run method>.shp
DEFUALT_FILENAME_STATION_RESULTS_SHP = 'stat_results_'  # + <LSM run method>.shp



# ----- SCHWAUS and DRIFT settings (legacy code) -----

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

FLAG_SAVE_DRIFT_PLOT_PDF = True
FLAG_CREATE_SCHWAUS_PROTOCOL = True
VERBOSE = True

# Output directory
OUT_PATH = '/home/heller/pyProjects/GravTools/out'

# Default names and paths of input files:

# ÖSGN Table:
PATH_OESGN_TABLE = '/home/heller/pyProjects/GravTools/data/'
NAME_OESGN_TABLE = 'OESGN.TAB'

# BEV observation files:
PATH_OBS_FILE_BEV = '/home/heller/pyProjects/GravTools/data/BEV/'
NAME_OBS_FILE_BEV = 'f200701_1'

# CG-5 observation files (text)
PATH_OBS_FILE_CG5 = '/home/heller/pyProjects/GravTools/data/CG5/'
NAME_OBS_FILE_CG5 = '20200907_test.TXT'

# ##### Debug options #####
DEBUG_TIME_IT = False