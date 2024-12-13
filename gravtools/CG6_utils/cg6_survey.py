"""Models for surveys with the Scintrex CG-6 gravity meters.

This module contains all models that are specifically required to
handle observation with the Scintrex CG-5 gravity meter.

Copyright (C) 2024  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

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
import pandas as pd
import numpy as np
import re
import datetime as dt

from numpy import datetime64

from gravtools.models.exceptions import InvaliFileContentError
from gravtools import settings


class DataCursor:
    """Data cursor for matplotlib plot. X and Y coordinates are printed life to the plot canvas.

    From: https://stackoverflow.com/questions/4652439/is-there-a-matplotlib-equivalent-of-matlabs-datacursormode
    """
    # text_template = 'x: %0.2f\ny: %0.2f'
    x, y = 0.0, 0.0
    xoffset, yoffset = -20, 20
    text_template = 'x: %0.2f\ny: %0.2f'

    def __init__(self, ax):
        self.ax = ax
        self.annotation = ax.annotate(self.text_template,
                                      xy=(self.x, self.y), xytext=(self.xoffset, self.yoffset),
                                      textcoords='offset points', ha='right', va='bottom',
                                      bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                                      arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
                                      )
        self.annotation.set_visible(False)

    def __call__(self, event):
        self.event = event
        # xdata, ydata = event.artist.get_data()
        # self.x, self.y = xdata[event.ind], ydata[event.ind]
        self.x, self.y = event.mouseevent.xdata, event.mouseevent.ydata
        if self.x is not None:
            self.annotation.xy = self.x, self.y
            self.annotation.set_text(self.text_template % (self.x, self.y))
            self.annotation.set_visible(True)
            event.canvas.draw()


class CG6Survey:
    """CG-6 survey data.

    Beside the actual observations class instances contain alls the header information provided in the CG-6 observation
    files.

    Attributes
    ----------
    obs_filename : str
        Name (and path) of the CG-6 observation file.
    obs_file_type : str (`cg6_obs_file_lynx` or `cg6_obs_file_solo`)
        Observation file type as described and listed in `settings.CG6_SURVEY_DATA_SOURCE_TYPES`.
    survey_name : str
        Name of the survey.
    serial_number : str
        Gravimeter serial number.
    created_datetime : datetime object
        Creation time and date of the observation file.
    operator : str
        Instrument operator.
    gcal1 : float
        Instrumental calibration factor GCAL1 of the Scintrex CG6 gravimeter as used for the survey [mGal].
    goff : float
        Gravity offset of the Scintrex CG6 used for the survey [ADU].
    gref : float
        Gravity reference value of the Scintrex CG6 used for the survey [mGal].
    x_scale : float
        X scaling factor of the Scintrex CG6 used for the survey [arcsec/ADU].
    y_scale : float
        Y scaling factor of the Scintrex CG6 used for the survey [arcsec/ADU].
    x_offset : float
        X offset value of the Scintrex CG6 used for the survey [ADU].
    y_offset : float
        X offset value of the Scintrex CG6 used for the survey [ADU].
    temp_corr_coeff : float
        Temperature correction coefficient of the Scintrex CG6 used for the survey [mGal/mK].
    temp_sensor_scale : float
        Temperature sensor scale value of the Scintrex CG6 used for the survey [mK/ADU].
    drift_rate : float
        Long-term drift value [mGal/day]
    drift_ref_datetime : datetime object
        Drift zero time and date.
    tilt_corr : bool
        Flag that indicates whether tilt correction was enabled.
    tidal_corr : bool
        Flag that indicates whether tidal correction was enabled.
    ocean_loading_corr : bool
        Flag that indicates whether ocean loading correction was enabled. Only available when using LynxLG.
    temp_corr : bool
        Flag that indicates whether temperature correction was enabled.
    drift_corr : bool
        Flag that indicates whether drift correction was enabled.
    obs_df : pandas data frame
       Contains the actual observation data records with the columns defined in `self._OBS_DF_COLUMN_NAMES`.
    temp_sensor_offset : float, optional (default = numpy.nan)
        Temperature sensor offset. LynxLG only.
    firmware_version : str, optional (default = '')
        Firmware version of the instrument. Only available in CG6-solo observation files.
    gravity_filter_type : str, optional (default = '')
        Type of the filter used to calculate filtered gravity values. LynxLG only.
    gravity_filter_period : float, optional (default = numpy.nan)
        Length of filter period [s]. LynxLG only.
    date_time_source : str, optional (default = '')
        Source of date and time. LynxLG only.
    tidal_corr_model : str, optional (default = '')
        Name of the tidal correction model. LynxLG only.
    lynxlg_version : str, optional (default = '')
        LynxLG software version. LynxLG only.
    """
    # Column names of the dataframe containing tha actual observation data:
    _OBS_DF_COLUMNS = {
        'ref_time': np.datetime64,
        'station': str,
        'occupation': str,
        'line': str,
        'user_lat_deg': float,
        'user_lon_deg': float,
        'user_height_m': float,
        'g_corr_mugal': float,
        'se_mugal': float,
        'sd_mugal': float,
        'tilt_x_arcsec': float,
        'tilt_y_arcsec': float,
        'sensor_temp_mk': float,
        'corr_tilt_mugal': float,
        'corr_tide_mugal': float,
        'corr_oceanload_mugal': float,
        'corr_temp_mugal': float,
        'corr_drift_mugal': float,
        'gps_lat_deg': float,
        'gps_lon_deg': float,
        'gps_fix_quality': float,
        'gps_satellites': float,
        'gps_hdop': float,
        'gps_h_m': float,
        'duration': float,
        'num_rejected': float,
        'instr_height_m': float,
        'gradient_mugalm': float,
        'g_raw_mugal': float,
        'setup_id': float,
        'note': str,
    }
    _OBS_DF_COLUMN_NAMES = _OBS_DF_COLUMNS.keys()

    def __init__(self,
                 obs_filename: str,
                 obs_file_type: str,
                 survey_name: str,
                 serial_number: str,
                 created_datetime: dt.datetime,
                 operator: str,
                 gcal1: float,
                 goff: float,
                 gref: float,
                 x_scale: float,
                 y_scale: float,
                 x_offset: float,
                 y_offset: float,
                 temp_corr_coeff: float,
                 temp_sensor_scale: float,
                 drift_rate: float,
                 drift_ref_datetime: dt.datetime,
                 tilt_corr: bool,
                 tidal_corr: bool,
                 ocean_loading_corr: bool,
                 temp_corr: bool,
                 drift_corr: bool,
                 obs_df: pd.DataFrame,
                 temp_sensor_offset: float = np.nan,
                 firmware_version: str = '',
                 gravity_filter_type: str = '',
                 gravity_filter_period: float = np.nan,
                 date_time_source: str = '',
                 tidal_corr_model: str = '',
                 lynxlg_version: str = '',
                 ):
        """
        Parameters
        ----------
        obs_filename : str, optional (default = '')
            Name (and path) of CG-6 observation file.
        obs_file_type : str, optional (default = '')
            Observation file type as described and listed in `settings.CG6_SURVEY_DATA_SOURCE_TYPES`.
        survey_name : str
            Name of the survey.
        serial_number : str
            Gravimeter serial number.
        created_datetime : datetime object
            Creation time and date of the observation file.
        operator : str
            Instrument operator.
        gcal1 : float
            Instrumental calibration factor GCAL1 of the Scintrex CG6 gravimeter as used for the survey [mGal].
        goff : float
            Gravity offset of the Scintrex CG6 used for the survey [ADU].
        gref : float
            Gravity reference value of the Scintrex CG6 used for the survey [mGal].
        x_scale : float
            X scaling factor of the Scintrex CG6 used for the survey [arcsec/ADU].
        y_scale : float
            Y scaling factor of the Scintrex CG6 used for the survey [arcsec/ADU].
        x_offset : float
            X offset value of the Scintrex CG6 used for the survey [ADU].
        y_offset : float
            X offset value of the Scintrex CG6 used for the survey [ADU].
        temp_corr_coeff : float
            Temperature correction coefficient of the Scintrex CG6 used for the survey [mGal/mK].
        temp_sensor_scale : float
            Temperature sensor scale value of the Scintrex CG6 used for the survey [mK/ADU].
        drift_rate : float
            Long-term drift value [mGal/day]
        drift_ref_datetime : datetime object
            Drift zero time and date.
        tilt_corr : bool
            Flag that indicates whether tilt correction was enabled.
        tidal_corr : bool
            Flag that indicates whether tidal correction was enabled.
        ocean_loading_corr : bool
            Flag that indicates whether ocean loading correction was enabled. Only available when using LynxLG.
        temp_corr : bool
            Flag that indicates whether temperature correction was enabled.
        drift_corr : bool
            Flag that indicates whether drift correction was enabled.
        obs_df : pandas data frame
            Contains the actual observation data records with the columns defined in `self._OBS_DF_COLUMN_NAMES`.
        temp_sensor_offset : float, optional (default = numpy.nan)
            Temperature sensor offset. LynxLG only.
        firmware_version : str, optional (default = '')
            Firmware version of the instrument. Only available in CG6-solo observation files.
        gravity_filter_type : str, optional (default = '')
            Type of the filter used to calculate filtered gravity values. LynxLG only.
        gravity_filter_period : float, optional (default = numpy.nan)
            Length of filter period [s]. LynxLG only.
         date_time_source : str, optional (default = '')
            Source of date and time. LynxLG only.
        tidal_corr_model : str, optional (default = '')
            Name of the tidal correction model. LynxLG only.
        lynxlg_version : str, optional (default = '')
            LynxLG software version. LynxLG only.
        """
        # obs_filename
        if not isinstance(obs_filename, str):
            raise RuntimeError('"obs_filename" needs to be a string.')
        self.obs_filename = obs_filename

        # obs_file_type
        if obs_file_type not in settings.CG6_SURVEY_DATA_SOURCE_TYPES:
            raise RuntimeError(f'Invalid observation file type.')
        self.obs_file_type = obs_file_type

        # survey_name
        if not isinstance(survey_name, str):
            raise RuntimeError('"survey_name" needs to be a string.')
        self.survey_name = survey_name

        # serial_number
        if not isinstance(serial_number, str):
            raise RuntimeError('"serial_number" needs to be a string.')
        self.serial_number = serial_number

        # created_datetime
        if not isinstance(created_datetime, dt.datetime):
            raise RuntimeError('"created_datetime" needs to be a datetime object.')
        self.created_datetime = created_datetime

        # operator
        if not isinstance(operator, str):
            raise RuntimeError('"operator" needs to be a string.')
        self.operator = operator

        # gcal1
        if not isinstance(gcal1, float):
            raise RuntimeError('"gcal1" needs to be a float.')
        self.gcal1 = gcal1

        # goff
        if not isinstance(goff, float):
            raise RuntimeError('"goff" needs to be a float.')
        self.goff = goff

        # gref
        if not isinstance(gref, float):
            raise RuntimeError('"gref" needs to be a float.')
        self.gref = gref

        # x_scale
        if not isinstance(x_scale, float):
            raise RuntimeError('"x_scale" needs to be a float.')
        self.x_scale = x_scale

        # y_scale
        if not isinstance(y_scale, float):
            raise RuntimeError('"y_scale" needs to be a float.')
        self.y_scale = y_scale

        # x_offset
        if not isinstance(x_offset, float):
            raise RuntimeError('"x_offset" needs to be a float.')
        self.x_offset = x_offset

        # y_offset
        if not isinstance(y_offset, float):
            raise RuntimeError('"y_offset" needs to be a float.')
        self.y_offset = y_offset

        # temp_corr_coeff
        if not isinstance(temp_corr_coeff, float):
            raise RuntimeError('"temp_corr_coeff" needs to be a float.')
        self.temp_corr_coeff = temp_corr_coeff

        # temp_sensor_scale
        if not isinstance(temp_sensor_scale, float):
            raise RuntimeError('"temp_sensor_scale" needs to be a float.')
        self.temp_sensor_scale = temp_sensor_scale

        # drift_rate
        if not isinstance(drift_rate, float):
            raise RuntimeError('"drift_rate" needs to be a float.')
        self.drift_rate = drift_rate

        # drift_ref_datetime
        if not isinstance(drift_ref_datetime, dt.datetime):
            raise RuntimeError('"drift_ref_datetime" needs to be a datetime object.')
        self.drift_ref_datetime = drift_ref_datetime

        # tilt_corr
        if not isinstance(tilt_corr, bool):
            raise RuntimeError('"tilt_corr" needs to be a bool type variable.')
        self.tilt_corr = tilt_corr

        # tidal_corr
        if not isinstance(tidal_corr, bool):
            raise RuntimeError('"tidal_corr" needs to be a bool type variable.')
        self.tidal_corr = tidal_corr

        # ocean_loading_corr
        if not isinstance(ocean_loading_corr, bool):
            raise RuntimeError('"ocean_loading_corr" needs to be a bool type variable.')
        self.ocean_loading_corr = ocean_loading_corr

        # temp_corr
        if not isinstance(temp_corr, bool):
            raise RuntimeError('"temp_corr" needs to be a bool type variable.')
        self.temp_corr = temp_corr

        # drift_corr
        if not isinstance(drift_corr, bool):
            raise RuntimeError('"drift_corr" needs to be a bool type variable.')
        self.drift_corr = drift_corr

        # check obs_df:
        if not isinstance(obs_df, pd.DataFrame):
            raise RuntimeError('"obs_df" needs to be a pandas DataFrame.')
        self._check_obs_df_columns(obs_df)
        obs_df = self._obs_df_reorder_columns(obs_df)

        self.obs_df = obs_df

        # Optional parameters:

        # temp_sensor_offset
        if not isinstance(temp_sensor_offset, float) and (temp_sensor_offset is not np.nan):
            raise RuntimeError('"temp_sensor_offset" needs to be a float.')
        self.temp_sensor_offset = temp_sensor_offset

        # firmware_version
        if not isinstance(firmware_version, str) and (firmware_version != ''):
            raise RuntimeError('"firmware_version" needs to be a string.')
        self.firmware_version = firmware_version

        # gravity_filter_type
        if not isinstance(gravity_filter_type, str) and (gravity_filter_type != ''):
            raise RuntimeError('"gravity_filter_type" needs to be a string.')
        self.gravity_filter_type = gravity_filter_type

        # gravity_filter_period
        if not isinstance(gravity_filter_period, float) and (gravity_filter_period is not np.nan):
            raise RuntimeError('"gravity_filter_period" needs to be a float.')
        self.gravity_filter_period = gravity_filter_period

        # date_time_source
        if not isinstance(date_time_source, str) and (date_time_source != ''):
            raise RuntimeError('"date_time_source" needs to be a string.')
        self.date_time_source = date_time_source

        # tidal_corr_model
        if not isinstance(tidal_corr_model, str) and (tidal_corr_model != ''):
            raise RuntimeError('"tidal_corr_model" needs to be a string.')
        self.tidal_corr_model = tidal_corr_model

        # lynxlg_version
        if not isinstance(lynxlg_version, str) and (lynxlg_version != ''):
            raise RuntimeError('"lynxlg_version" needs to be a string.')
        self.lynxlg_version = lynxlg_version

    @staticmethod
    def read_cg6_file(filename, comment_marker: str = '#') -> tuple[str, list[str]]:
        """Returns the content of the input file as string and as a list of line strings.

        Notes
        -----
        Each line of the input file is returned with exactly one EOL symbol in the return string.

        Parameters
        ----------
        filename : str
            Name of the input file.
        comment_marker : str, optional (default = '#')
            All lines that begin with `comment_marker` are supposed to be comments and are skipped.

        Returns
        -------
            tuple[str, list[str]] : Content of the input file as one string and as list of line strings.
        """
        # Read in file and ignore comment lines:
        file_handle = open(filename, 'r')
        lines = [file_handle.readline().strip()]
        if lines[0].startswith('\ufeff'):
            # Check byte order mark and change codec. Reopen file if required with the correct codec.
            file_handle.close()
            file_handle = open(filename, 'r', encoding='utf-8-sig')
            lines = []
        for line in file_handle:
            line_tmp = line.strip()
            if not line_tmp.startswith(comment_marker):
                lines.append(line_tmp)
        file_handle.close()
        str_obs_file = '\n'.join(lines)

        # Remove all or add one end-of-line symbols from end of string:
        # - Last character of string has to be a \n so that regex works correctly!
        number_of_trailing_eol_symbols = 0
        str_idx = -1
        while str_obs_file[str_idx] == '\n':
            number_of_trailing_eol_symbols += 1
            str_idx -= 1
        if number_of_trailing_eol_symbols == 0:
            str_obs_file += '\n'
        elif number_of_trailing_eol_symbols > 1:
            str_obs_file = str_obs_file[0:str_idx + 2]

        return str_obs_file, lines

    @classmethod
    def from_lynxlg_file(cls, filename: str, obs_file_type: str, expect_notes: bool,
                         dt_setup_sec: float = None, verbose: bool = True):
        """Constractor for LynxLG formatted observation files (version 1 and version 2).

        Notes
        -----
        The difference to the V1 observation file format is, that the V2 format supports notes taken per setup. Notes
        are saved in the file right before the block of observations at a station. Be aware that the note for the first
        setup is written to the file header.

        Subsequent setups are supposed to be separated by note entries in the observation file or by station names.
        Optionally, a time span may be defined that is used to distinguish between consecutive setups. In the latter
        case observations are divided into different setups, if the time gap between them is larger than the defined
        value.

        Parameters
        ----------
        filename : str
            Name and path of the observation file.
        obs_file_type : str
            Type of observation file format. All valid identifiers are listed as keys in
            `settings.CG6_SURVEY_DATA_SOURCE_TYPES`.
        expect_notes : bool
            `True` implies that notes are mandatory in the observation file. Notes are only available in the LYnxLG
            format version 2.
        dt_setup_sec : float, optional (default = `None`)
            Observations are split up into separate setups, if the time gap between them is larger than the provided
            value (in seconds).
        verbose : bool, optional (default = `True`)
            True indicates that status messages are written to the command line.
        """

        _HEADER_LINES = {
            # <variable_name str>: [<regex_expr str>, <field_name str>, <flag_optional bool>, <default var>, <dtype>]
            'created_date': [r'\/ Date: (?P<input_str>\S+)\s*\n', 'Date', False, '', 'str'],
            'created_time': [r'\/ Time: (?P<input_str>\S+)\s*\n', 'Time', False, '', 'str'],
            'meter': [r'\/ Meter: (?P<input_str>\S+)\s*\n', 'Meter', False, '', 'str'],
            'operator': [r'\/ Operator: (?P<input_str>\S+)\s*\n', 'Operator', False, '', 'str'],
            'gcal1': [r'\/ GCAL1: (?P<input_str>\S+)\s*\n', 'GCAL1', False, np.nan, 'float'],
            'goff': [r'\/ Gravity Offset: (?P<input_str>\S+)\s*\n', 'Gravity Offset', False, np.nan, 'float'],
            'gref': [r'\/ Gravity Reference: (?P<input_str>\S+)\s*\n', 'Gravity Reference', False, np.nan, 'float'],
            'temp_sensor_scale': [r'\/ Sensor Temperature Gain: (?P<input_str>\S+)\s*\n', 'Sensor Temperature Gain',
                                  False, np.nan, 'float'],
            'temp_sensor_offset': [r'\/ Sensor Temperature Offset: (?P<input_str>\S+)\s*\n',
                                   'Sensor Temperature Offset', False, np.nan, 'float'],
            'x_scale': [r'\/ X Level Gain: (?P<input_str>\S+)\s*\n', 'X Level Gain', False, np.nan, 'float'],
            'x_offset': [r'\/ X Level Offset: (?P<input_str>\S+)\s*\n', 'X Level Offset', False, np.nan, 'float'],
            'y_scale': [r'\/ Y Level Gain: (?P<input_str>\S+)\s*\n', 'Y Level Gain', False, np.nan, 'float'],
            'y_offset': [r'\/ Y Level Offset: (?P<input_str>\S+)\s*\n', 'Y Level Offset', False, np.nan, 'float'],
            'gravity_filter_type': [r'\/ Gravity Filter: (?P<input_str>.+)\n', 'Gravity Filter', False, '', 'str'],
            'gravity_filter_period': [r'\/ Gravity Filter Period: (?P<input_str>\S+)\s*\n', 'Gravity Filter Period',
                                      False, np.nan, 'float'],
            'date_time_source': [r'\/ Date/Time Source: (?P<input_str>\S+)\s*\n', 'Date/Time Source', False, '', 'str'],
            'tilt_corr': [r'\/ Level Correction: (?P<input_str>\S+)\s*\n', 'Level Correction', False, '', 'str'],
            'tidal_corr': [r'\/ Tidal Correction: (?P<input_str>\S+)\s*\n', 'Tidal Correction', False, '', 'str'],
            'ocean_loading_corr': [r'\/ Ocean Load Correction: (?P<input_str>\S+)\s*\n', 'Ocean Load Correction', False,
                                   '', 'str'],
            'temp_corr': [r'\/ Temperature Correction: (?P<input_str>\S+)\s*\n', 'Temperature Correction', False, '',
                          'str'],
            'drift_corr': [r'\/ Drift Correction: (?P<input_str>\S+)\s*\n', 'Drift Correction', False, '', 'str'],
            'tidal_corr_model': [r'\/ Tidal Model: (?P<input_str>\S+)\s*\n', 'Tidal Model', False, '', 'str'],
            'temp_corr_coeff': [r'\/ Temperature Correction Coefficient \(mGals\/mK\): (?P<input_str>\S+)\s*\n',
                                'Temperature Correction Coefficient (mGals/mK)', False, np.nan, 'float'],
            'drift_rate': [r'\/ Drift Rate \(mGals\/day\): (?P<input_str>\S+)\s*\n', 'Drift Rate (mGals/day)', False,
                           np.nan, 'float'],
            'drift_ref_date': [r'\/ Drift Zero Date: (?P<input_str>\S+)\s*\n', 'Drift Zero Date', False, '', 'str'],
            'drift_ref_time': [r'\/ Drift Zero Time: (?P<input_str>\S+)\s*\n', 'Drift Zero Time', False, '', 'str'],
        }

        if verbose:
            print(f'Load CG6 observation file (LYnxLG, V2): {filename}')

        # Read file:
        file_str, lines = cls.read_cg6_file(filename)

        # Read header items:
        # - This item is treated separately, because it can occure more than once in the obs. file.
        expr = r'\/ Software Version: (?P<lynxlg_version>\S+)\s*\n'
        # Read setup parameters blocks from obs file string (only one block allowed!):
        lynxlg_version_list = []
        for block in re.finditer(expr, file_str):
            block_dict = block.groupdict()
            lynxlg_version_list.append(block_dict['lynxlg_version'])
        if len(set(lynxlg_version_list)) > 1:
            raise InvaliFileContentError(f'The file {filename} contains non-unique software version numbers!')
        lynxlg_version = lynxlg_version_list[0]

        # Read header strings:
        for key, item in _HEADER_LINES.items():
            expr = item[0]
            field_name = item[1]
            flag_optional = item[2]
            data_type = item[4]
            block_count = 0
            for block in re.finditer(expr, file_str):
                block_dict = block.groupdict()
                block_count += 1
            if block_count == 1:
                if data_type == 'str':
                    item[3] = block_dict['input_str']
                elif data_type == 'float':
                    item[3] = float(block_dict['input_str'])
                elif data_type == 'int':
                    item[3] = int(block_dict['input_str'])
                else:
                    raise RuntimeError(f'Unknown data type "{data_type}": Conversion of field "{field_name}" failed '
                                       f'when loading file {filename}.')
            elif block_count == 0:
                if not flag_optional:
                    raise InvaliFileContentError(f'The file {filename} does not contain a "{field_name}" header line.')
            if block_count > 1:
                raise InvaliFileContentError(f'The file {filename} contains more than one "{field_name}" header line.')
            # _HEADER_LINES[key] = item  # Not required due to pass by reference

        # Convert strings to bools:
        bool_lookup_dict = {'Enabled': True, 'Disabled': False}
        convert_to_bool_items = ['tilt_corr', 'tidal_corr', 'ocean_loading_corr', 'temp_corr', 'drift_corr']
        for item in convert_to_bool_items:
            header_str = _HEADER_LINES[item][3]
            try:
                _HEADER_LINES[item][3] = bool_lookup_dict[header_str]
            except KeyError:
                raise RuntimeError(f'Invalid value in header field "{_HEADER_LINES[item][1]}": {header_str}.')

        # Convert strings to datetime:
        datetime_str = _HEADER_LINES['drift_ref_date'][3] + ' ' + _HEADER_LINES['drift_ref_time'][3]
        drift_ref_datetime = dt.datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')
        datetime_str = _HEADER_LINES['created_date'][3] + ' ' + _HEADER_LINES['created_time'][3]
        created_datetime = dt.datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')

        # Read observations:
        ref_time_list = []
        station_list = []
        occupation_list = []
        line_list = []
        user_lat_deg_list = []
        user_lon_deg_list = []
        user_height_m_list = []
        g_corr_mugal_list = []
        se_mugal_list = []
        sd_mugal_list = []
        tilt_x_arcsec_list = []
        tilt_y_arcsec_list = []
        sensor_temp_mk_list = []
        corr_tilt_mugal_list = []
        corr_tide_mugal_list = []
        corr_oceanload_mugal_list = []
        corr_temp_mugal_list = []
        corr_drift_mugal_list = []
        gps_lat_deg_list = []
        gps_lon_deg_list = []
        gps_fix_quality_list = []
        gps_satellites_list = []
        gps_hdop_list = []
        gps_h_m_list = []
        duration_list = []
        num_rejected_list = []
        instr_height_m_list = []
        gradient_mugalm_list = []
        g_raw_mugal_list = []
        setup_id_list = []
        note_list = []
        survey_name_list = []

        if expect_notes:
            flag_first_note = False
            flag_new_setup = False
        else:
            flag_first_note = True
            flag_new_setup = True
            note = ''  # Reset note
        num_setups = 0
        num_obs = 0
        for line in lines:
            # print(line)
            if line.startswith('/ Notes:'):
                flag_new_setup = True
                flag_first_note = True
                note = line.split('/ Notes:')[1].strip()
            elif not line.startswith('/') and line and flag_first_note:  # obs. line
                expr = r'(?P<year>[0-9]{4}) (?P<month>[ 0-9][0-9]) (?P<day>[ 0-9][0-9]) (?P<hour>[ 0-9][0-9]) (?P<minute>[ 0-9][0-9]) (?P<second>[ 0-9][0-9]) (?P<millisecond>[ 0-9][ 0-9][0-9]) (?P<survey>\S+) (?P<station>\S+) (?P<occupation>\S*) (?P<line>\S*) (?P<user_lat_deg>-?[0-9.]+) (?P<user_lon_deg>-?[0-9.]+) (?P<user_height_m>-?[0-9.]+) (?P<corr_g_mgal>[0-9.]+) (?P<se_g_mgal>[0-9.]+) (?P<sd_g_mgal>[0-9.]+) (?P<tilt_x_arcsec>-?[0-9.]+) (?P<tilt_y_arcsec>-?[0-9.]+) (?P<sensor_temp_mk>-?[0-9.]+) (?P<corr_tilt_mgal>-?[0-9.]+) (?P<corr_tide_mgal>-?[0-9.]+) (?P<corr_oceanload_mgal>-?[0-9.]+) (?P<corr_temp_mgal>-?[0-9.]+) (?P<corr_drift_mgal>-?[0-9.]+) (?P<gps_lat_deg>-?[0-9.]+) (?P<gps_lon_deg>-?[0-9.]+) +(?P<gps_fix_quality>[0-9]+) +(?P<gps_satellites>[0-9]+) (?P<gps_hdop>-?[0-9.]+) (?P<gps_h_m>-?[0-9.]+) +(?P<duration>[0-9]+) +(?P<num_rejected>[0-9]+) (?P<instr_height_m>[0-9.]+) (?P<gradient_mgalcm>-?[0-9.]+)'
                block_count = 0
                for block in re.finditer(expr, line):
                    block_dict = block.groupdict()
                    block_count += 1
                if block_count == 0:
                    raise InvaliFileContentError(
                        f'Invalid observation line in file {filename}: {line}. Could not match with regex.')
                elif block_count > 1:
                    raise InvaliFileContentError(
                        f'Invalid observation line in file {filename}: {line}. Multiple matches with regex.')

                # One match => OK:
                num_obs += 1

                t_ref = dt.datetime(year=int(block_dict['year']), month=int(block_dict['month']),
                                    day=int(block_dict['day']),
                                    hour=int(block_dict['hour']), minute=int(block_dict['minute']),
                                    second=int(block_dict['second']), microsecond=int(block_dict['millisecond']) * 1000)
                ref_time_list.append(t_ref)

                station_list.append(block_dict['station'])
                occupation_list.append(block_dict['occupation'])
                line_list.append(block_dict['line'])
                user_lat_deg_list.append(float(block_dict['user_lat_deg']))
                user_lon_deg_list.append(float(block_dict['user_lon_deg']))
                user_height_m_list.append(float(block_dict['user_height_m']))
                g_corr_mugal_list.append(float(block_dict['corr_g_mgal']) * 1e3)
                se_mugal_list.append(float(block_dict['se_g_mgal']) * 1e3)
                sd_mugal_list.append(float(block_dict['sd_g_mgal']) * 1e3)
                tilt_x_arcsec_list.append(float(block_dict['tilt_x_arcsec']))
                tilt_y_arcsec_list.append(float(block_dict['tilt_y_arcsec']))
                sensor_temp_mk_list.append(float(block_dict['sensor_temp_mk']))
                corr_tilt_mugal_list.append(float(block_dict['corr_tilt_mgal']) * 1e3)
                corr_tide_mugal_list.append(float(block_dict['corr_tide_mgal']) * 1e3)
                corr_oceanload_mugal_list.append(float(block_dict['corr_oceanload_mgal']) * 1e3)
                corr_temp_mugal_list.append(float(block_dict['corr_temp_mgal']) * 1e3)
                corr_drift_mugal_list.append(float(block_dict['corr_drift_mgal']) * 1e3)
                gps_lat_deg_list.append(float(block_dict['gps_lat_deg']))
                gps_lon_deg_list.append(float(block_dict['gps_lon_deg']))
                gps_fix_quality_list.append(float(block_dict['gps_fix_quality']))
                gps_satellites_list.append(float(block_dict['gps_satellites']))
                gps_hdop_list.append(float(block_dict['gps_hdop']))
                gps_h_m_list.append(float(block_dict['gps_h_m']))
                duration_list.append(float(block_dict['duration']))
                num_rejected_list.append(float(block_dict['num_rejected']))
                instr_height_m_list.append(float(block_dict['instr_height_m']))
                gradient_mugalm_list.append(float(block_dict[
                                                      'gradient_mgalcm']) * 1e5)  # TODO: Unit mgal/cm makes no sense...
                g_raw_mugal_list.append(np.nan)
                survey_name_list.append(block_dict['survey'])

                # Check conditions for new setups:
                if num_obs > 1:  # Minimum 2nd observation in survey
                    if station_list[num_obs - 1] != station_list[num_obs - 2]:  # Check station names to create setups
                        flag_new_setup = True
                        note = ''  # Reset note
                    elif dt_setup_sec is not None:  # Check time gap
                        if (ref_time_list[num_obs - 1] - ref_time_list[num_obs - 2]).seconds > dt_setup_sec:
                            flag_new_setup = True
                            note = ''  # Reset note

                # # Check time gap between observations to create setups:
                # if dt_setup_sec is not None:
                #     if not flag_new_setup:  # At least one obs. after a "/ Notes" line
                #         if num_obs > 0:
                #             if (ref_time_list[num_obs - 1] - ref_time_list[num_obs - 2]).seconds > dt_setup_sec:
                #                 flag_new_setup = True
                #                 note = ''  # Reset note

                # New setup id:
                if flag_new_setup:
                    flag_new_setup = False
                    setup_id = int(dt.datetime.timestamp(t_ref))
                    num_setups += 1
                setup_id_list.append(setup_id)
                note_list.append(note)

        # Check, whether data was loaded from file:
        if num_obs == 0:
            raise RuntimeError(f'No observations loaded from CG6 file ({filename}). Check, whether the correct file '
                               f'format was selected!')

        # Get survey name from observation lines:
        if len(set(survey_name_list)) > 1:
            raise RuntimeError(f'The survey name is not unique in the observation file {filename}!')
        survey_name = set(survey_name_list).pop()

        # Create pandas dataframe:
        obs_df = pd.DataFrame(list(zip(ref_time_list,
                                       station_list,
                                       occupation_list,
                                       line_list,
                                       user_lat_deg_list,
                                       user_lon_deg_list,
                                       user_height_m_list,
                                       g_corr_mugal_list,
                                       se_mugal_list,
                                       sd_mugal_list,
                                       tilt_x_arcsec_list,
                                       tilt_y_arcsec_list,
                                       sensor_temp_mk_list,
                                       corr_tilt_mugal_list,
                                       corr_tide_mugal_list,
                                       corr_oceanload_mugal_list,
                                       corr_temp_mugal_list,
                                       corr_drift_mugal_list,
                                       gps_lat_deg_list,
                                       gps_lon_deg_list,
                                       gps_fix_quality_list,
                                       gps_satellites_list,
                                       gps_hdop_list,
                                       gps_h_m_list,
                                       duration_list,
                                       num_rejected_list,
                                       instr_height_m_list,
                                       gradient_mugalm_list,
                                       g_raw_mugal_list,
                                       setup_id_list,
                                       note_list)),
                              columns=cls._OBS_DF_COLUMN_NAMES)

        return cls(obs_filename=filename,
                   obs_file_type=obs_file_type,
                   survey_name=survey_name,  # From observation block. Check for equality!
                   serial_number=_HEADER_LINES['meter'][3],
                   created_datetime=created_datetime,
                   operator=_HEADER_LINES['operator'][3],
                   gcal1=_HEADER_LINES['gcal1'][3],
                   goff=_HEADER_LINES['goff'][3],
                   gref=_HEADER_LINES['gref'][3],
                   x_scale=_HEADER_LINES['x_scale'][3],
                   y_scale=_HEADER_LINES['y_scale'][3],
                   x_offset=_HEADER_LINES['x_offset'][3],
                   y_offset=_HEADER_LINES['y_offset'][3],
                   temp_corr_coeff=_HEADER_LINES['temp_corr_coeff'][3],
                   temp_sensor_scale=_HEADER_LINES['temp_sensor_scale'][3],
                   drift_rate=_HEADER_LINES['drift_rate'][3],
                   drift_ref_datetime=drift_ref_datetime,
                   tilt_corr=_HEADER_LINES['tilt_corr'][3],
                   tidal_corr=_HEADER_LINES['tidal_corr'][3],
                   ocean_loading_corr=_HEADER_LINES['ocean_loading_corr'][3],
                   temp_corr=_HEADER_LINES['temp_corr'][3],
                   drift_corr=_HEADER_LINES['drift_corr'][3],
                   obs_df=obs_df,
                   temp_sensor_offset=_HEADER_LINES['temp_sensor_offset'][3],
                   # firmware_version: str, optional(default='')
                   gravity_filter_type=_HEADER_LINES['gravity_filter_type'][3],
                   gravity_filter_period=_HEADER_LINES['gravity_filter_period'][3],
                   date_time_source=_HEADER_LINES['date_time_source'][3],
                   tidal_corr_model=_HEADER_LINES['tidal_corr_model'][3],
                   lynxlg_version=lynxlg_version)

    @classmethod
    def from_lynxlg_file_v1(cls, filename: str, dt_setup_sec: float = None, verbose: bool = True):
        """Constractor for LynxLG formatted observation files (version 1, without notes).

        Notes
        -----
        The difference to the V1 observation file format is, that the V2 format supports notes taken per setup. Notes
        are saved in the file right before the block of observations at a station. Be aware that the note for the first
        setup is written to the file header.

        Subsequent setups at the same station can be separated based on the time difference between observations. If the
        separation of the reference time of two observations is larger than the value given by `dt_setup_sec` in seconds
        they are assumed to belong to two consecutive setups.

        Parameters
        ----------
        filename : str
            Name and path of the observation file.
        dt_setup_sec : float, optional (default = `None`)
            Minimum time gap between two consecutive observations [sec] in order to create separate setups. `None`
            implies that time gaps are not considered in order to create setups.
        verbose : bool, optional (default = True)
            True indicates that status messages are written to the command line.
        """
        return cls.from_lynxlg_file(filename, expect_notes=False, obs_file_type='cg6_obs_file_lynx_v1',
                                    dt_setup_sec=dt_setup_sec, verbose=verbose)

    @classmethod
    def from_lynxlg_file_v2(cls, filename: str, dt_setup_sec: float = None, verbose: bool = True):
        """Constractor for LynxLG formatted observation files (version 2, with notes).

        Notes
        -----
        The difference to the V1 observation file format is, that the V2 format supports notes taken per setup. Notes
        are saved in the file right before the block of observations at a station. Be aware that the note for the first
        setup is written to the file header.

        Subsequent setups are supposed to be separated by note entries in the observation file and by station names.
        Optionally, a time span may be defined that is used to distinguish between consecutive setups. In the latter
        case observations are divided into different setups, if the time gap between them is larger than the defined
        value.

        Parameters
        ----------
        filename : str
            Name and path of the observation file.
        dt_setup_sec : float, optional (default = `None`)
            Minimum time gap between two consecutive observations [sec] in order to create separate setups. `None`
            implies that time gaps are not considered in order to create setups.
        verbose : bool, optional (default = `True`)
            True indicates that status messages are written to the command line.
        """
        return cls.from_lynxlg_file(filename, obs_file_type='cg6_obs_file_lynx_v2', expect_notes=True,
                                    dt_setup_sec=dt_setup_sec, verbose=verbose)

    @classmethod
    def from_cg6solo_file(cls, filename: str, dt_setup_sec: float = None, verbose: bool = True):
        """Constractor for CG6-solo formatted observation files (version 2).

        Notes
        -----
        Optionally, subsequent setups at the same station can be separated based on the time difference between
        observations. If the separation of the reference time of two observations is larger than the value given by
        `dt_setup_sec` in seconds they are assumed to belong to two consecutive setups.

        Parameters
        ----------
        filename : str
            Name and path of the observation file.
        dt_setup_sec : float, optional (default = `None`)
            Minimum time gap between two consecutive observations [sec] in order to create separate setups. `None`
            implies that time gaps are not considered in order to create setups.
        verbose : bool, optional (default = True)
            True indicates that status messages are written to the command line.
        """
        _HEADER_LINES = {
            # <variable_name str>: [<regex_expr str>, <field_name str>, <flag_optional bool>, <default var>, <dtype>]
            'survey_name': [r'\/\t\tSurvey Name:\t(?P<input_str>\S+)\n', 'Survey name', False, '', 'str'],
            'serial_number': [r'\/\t\tInstrument Serial Number:\t(?P<input_str>\S+)\n', 'Instrument Serial Number',
                              False, '', 'str'],
            'created_datetime': [r'\/\t\tCreated:\t(?P<input_str>[0-9 :-]+)\n', 'Created', False, '', 'str'],
            'operator': [r'\/\t\tOperator:\t(?P<input_str>\S+)\n', 'Operator', False, '', 'str'],
            'gcal1': [r'\/\t\tGcal1 \[mGal\]:\t(?P<input_str>\S+)\n', 'Gcal1 [mGal]', False, np.nan, 'float'],
            'goff': [r'\/\t\tGoff \[ADU\]:\t(?P<input_str>\S+)\n', 'Goff [ADU]', False, np.nan, 'float'],
            'gref': [r'\/\t\tGref \[mGal\]:\t(?P<input_str>\S+)\n', 'Gref [mGal]', False, np.nan, 'float'],
            'x_scale': [r'\/\t\tX Scale \[arc-sec\/ADU\]:\t(?P<input_str>\S+)\n', 'X Scale [arc-sec/ADU]', False,
                        np.nan, 'float'],
            'y_scale': [r'\/\t\tY Scale \[arc-sec\/ADU\]:\t(?P<input_str>\S+)\n', 'Y Scale [arc-sec/ADU]', False,
                        np.nan, 'float'],
            'x_offset': [r'\/\t\tX Offset \[ADU\]:\t(?P<input_str>\S+)\n', 'X Offset [ADU]', False, np.nan, 'float'],
            'y_offset': [r'\/\t\tY Offset \[ADU\]:\t(?P<input_str>\S+)\n', 'Y Offset [ADU]', False, np.nan, 'float'],
            'temp_corr_coeff': [r'\/\t\tTemperature Coefficient \[mGal\/mK\]:\t(?P<input_str>\S+)\n',
                                'Temperature Coefficient [mGal/mK]', False, np.nan, 'float'],
            'temp_sensor_scale': [r'\/\t\tTemperature Scale \[mK\/ADU\]:\t(?P<input_str>\S+)\n',
                                  'Temperature Scale [mK/ADU]',
                                  False, np.nan, 'float'],
            'drift_rate': [r'\/\t\tDrift Rate \[mGal\/day\]:\t(?P<input_str>\S+)\n', 'Drift Rate [mGal/day]', False,
                           np.nan, 'float'],
            'drift_ref_datetime': [r'\/\t\tDrift Zero Time:\t(?P<input_str>[0-9 :-]+)\n', 'Drift Zero Time', False, '',
                                   'str'],
            'firmware_version': [r'\/\t\tFirmware Version:\t(?P<input_str>\S+)\n', 'Firmware Version', False, '', 'str'],
        }

        if verbose:
            print(f'Load CG6 observation file (CG6 solo): {filename}')
            if dt_setup_sec is not None:
                print(f' - A time gap of {dt_setup_sec} seconds or more between observations creates a new setup.')

        # Read file:
        file_str, lines = cls.read_cg6_file(filename)

        # Read header strings:
        for key, item in _HEADER_LINES.items():
            expr = item[0]
            field_name = item[1]
            flag_optional = item[2]
            data_type = item[4]
            block_count = 0
            for block in re.finditer(expr, file_str):
                block_dict = block.groupdict()
                block_count += 1
            if block_count == 1:
                if data_type == 'str':
                    item[3] = block_dict['input_str']
                elif data_type == 'float':
                    item[3] = float(block_dict['input_str'])
                elif data_type == 'int':
                    item[3] = int(block_dict['input_str'])
                else:
                    raise RuntimeError(
                        f'Unknown data type "{data_type}": Conversion of field "{field_name}" failed '
                        f'when loading file {filename}.')
            elif block_count == 0:
                if not flag_optional:
                    raise InvaliFileContentError(
                        f'The file {filename} does not contain a "{field_name}" header line.')
            if block_count > 1:
                raise InvaliFileContentError(
                    f'The file {filename} contains more than one "{field_name}" header line.')
            # _HEADER_LINES[key] = item  # Not required due to pass by reference

        # Convert strings to datetime:
        drift_ref_datetime = dt.datetime.strptime(_HEADER_LINES['drift_ref_datetime'][3], '%Y-%m-%d %H:%M:%S')
        created_datetime = dt.datetime.strptime(_HEADER_LINES['created_datetime'][3], '%Y-%m-%d %H:%M:%S')

        # Read observations:
        station_list = []
        ref_time_list = []
        g_corr_mugal_list = []
        line_list = []
        sd_mugal_list = []
        se_mugal_list = []
        g_raw_mugal_list = []
        tilt_x_arcsec_list = []
        tilt_y_arcsec_list = []
        sensor_temp_mk_list = []
        corr_tide_mugal_list = []
        corr_tilt_mugal_list = []
        corr_temp_mugal_list = []
        corr_drift_mugal_list = []
        duration_list = []
        instr_height_m_list = []
        user_lat_deg_list = []
        user_lon_deg_list = []
        user_height_m_list = []
        gps_lat_deg_list = []
        gps_lon_deg_list = []
        gps_h_m_list = []
        correction_flags_list = []
        occupation_list = []
        corr_oceanload_mugal_list = []
        gradient_mugalm_list = []
        gps_fix_quality_list = []
        gps_satellites_list = []
        gps_hdop_list = []
        num_rejected_list = []
        note_list = []
        setup_id_list = []

        # expr = r'(?P<station>\S+)\t(?P<date>[0-9-]{10})\t(?P<time>[0-9:]{8})\t(?P<corr_g_mgal>[0-9.]+)\t(?P<line>\S*)\t(?P<sd_g_mgal>[0-9.]+)\t(?P<se_g_mgal>[0-9.]+)\t(?P<raw_g_mgal>[0-9.]+)\t(?P<tilt_x_arcsec>-?[0-9.]+)\t(?P<tilt_y_arcsec>-?[0-9.]+)\t(?P<sensor_temp_mk>-?[0-9.]+)\t(?P<corr_tide_mgal>-?[0-9.]+)\t(?P<corr_tilt_mgal>-?[0-9.]+)\t(?P<corr_temp_mgal>-?[0-9.]+)\t(?P<corr_drift_mgal>-?[0-9.]+)\t(?P<duration>[0-9]+)\t(?P<instr_height_m>-?[0-9.]+)\t(?P<user_lat_deg>-?[0-9.]+)\t(?P<user_lon_deg>-?[0-9.]+)\t(?P<user_height_m>-?[0-9.]+)\t(?P<gps_lat_deg>-?[0-9.]+)\t(?P<gps_lon_deg>-?[0-9.]+)\t(?P<gps_height_m>-?[0-9.]+)\t(?P<correction_flags>[01]{5})\s*'
        expr = r'(?P<station>\S+)\t(?P<date>[0-9-]{10})\t(?P<time>[0-9:]{8})\t(?P<corr_g_mgal>[0-9.]+)\t(?P<line>\S*)\t(?P<sd_g_mgal>[0-9.]+)\t(?P<se_g_mgal>[0-9.]+)\t(?P<raw_g_mgal>[0-9.]+)\t(?P<tilt_x_arcsec>-?[0-9.]+)\t(?P<tilt_y_arcsec>-?[0-9.]+)\t(?P<sensor_temp_mk>-?[0-9.]+)\t(?P<corr_tide_mgal>-?[0-9.]+)\t(?P<corr_tilt_mgal>-?[0-9.]+)\t(?P<corr_temp_mgal>-?[0-9.]+)\t(?P<corr_drift_mgal>-?[0-9.]+)\t(?P<duration>[0-9]+)\t(?P<instr_height_m>-?[0-9.]+)\t(?P<user_lat_deg>-?[0-9.]+)\t(?P<user_lon_deg>-?[0-9.]+)\t(?P<user_height_m>-?[0-9.]+)\t(?P<gps_lat_deg>-?[0-9.-]+)\t(?P<gps_lon_deg>-?[0-9.-]+)\t(?P<gps_height_m>-?[0-9.-]+)\t(?P<correction_flags>[01]{5})\s*'
        obs_count = 0
        for block in re.finditer(expr, file_str):
            block_dict = block.groupdict()
            station_list.append(block_dict['station'])
            datetime_str = block_dict['date'] + ' ' + block_dict['time']
            ref_time_list.append(dt.datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S'))
            g_corr_mugal_list.append(float(block_dict['corr_g_mgal']) * 1000)
            line_list.append(block_dict['line'])
            sd_mugal_list.append(float(block_dict['sd_g_mgal']) * 1000)
            se_mugal_list.append(float(block_dict['se_g_mgal']) * 1000)
            g_raw_mugal_list.append(float(block_dict['raw_g_mgal']) * 1000)
            tilt_x_arcsec_list.append(float(block_dict['tilt_x_arcsec']))
            tilt_y_arcsec_list.append(float(block_dict['tilt_y_arcsec']))
            sensor_temp_mk_list.append(float(block_dict['sensor_temp_mk']))
            corr_tide_mugal_list.append(float(block_dict['corr_tide_mgal']) * 1000)
            corr_tilt_mugal_list.append(float(block_dict['corr_tilt_mgal']) * 1000)
            corr_temp_mugal_list.append(float(block_dict['corr_temp_mgal']) * 1000)
            corr_drift_mugal_list.append(float(block_dict['corr_drift_mgal']) * 1000)
            duration_list.append(float(block_dict['duration']))  # TODO: [sec]?
            instr_height_m_list.append(float(block_dict['instr_height_m']))
            user_lat_deg_list.append(float(block_dict['user_lat_deg']))
            user_lon_deg_list.append(float(block_dict['user_lon_deg']))
            user_height_m_list.append(float(block_dict['user_height_m']))
            try:
                gps_lat_deg_list.append(float(block_dict['gps_lat_deg']))
                gps_lon_deg_list.append(float(block_dict['gps_lon_deg']))
                gps_h_m_list.append(float(block_dict['gps_height_m']))
            except ValueError:  # In case no GPS signal was available/recorded, e.g. at indoor use
                gps_lat_deg_list.append(float('nan'))
                gps_lon_deg_list.append(float('nan'))
                gps_h_m_list.append(float('nan'))
            correction_flags_list.append(block_dict['correction_flags'])
            # Initialize fields that are not available in the obs. file:
            occupation_list.append('')
            corr_oceanload_mugal_list.append(np.nan)
            gradient_mugalm_list.append(np.nan)
            gps_fix_quality_list.append(np.nan)
            gps_satellites_list.append(np.nan)
            gps_hdop_list.append(np.nan)
            num_rejected_list.append(np.nan)
            note_list.append('')

            obs_count += 1
        if obs_count == 0:
            raise RuntimeError(f'No observations found in file {filename}.')

        # Separate instrument setups based on the time gaps and station names and determine setup IDs:
        num_setups = 0
        flag_new_setup = False
        for count_idx, ref_time in enumerate(ref_time_list):
            if count_idx == 0:  # First observation
                flag_new_setup = True
            else:
                if station_list[count_idx] != station_list[count_idx - 1]:  # New station?
                    flag_new_setup = True
                elif dt_setup_sec is not None:
                    if (ref_time_list[count_idx] - ref_time_list[count_idx - 1]).seconds > dt_setup_sec:  # Time gap?
                        flag_new_setup = True
            if flag_new_setup:
                setup_id = int(dt.datetime.timestamp(ref_time))
                num_setups += 1
                flag_new_setup = False
            setup_id_list.append(setup_id)

        # Create pandas dataframe:
        obs_df = pd.DataFrame(list(zip(ref_time_list,
                                       station_list,
                                       occupation_list,
                                       line_list,
                                       user_lat_deg_list,
                                       user_lon_deg_list,
                                       user_height_m_list,
                                       g_corr_mugal_list,
                                       se_mugal_list,
                                       sd_mugal_list,
                                       tilt_x_arcsec_list,
                                       tilt_y_arcsec_list,
                                       sensor_temp_mk_list,
                                       corr_tilt_mugal_list,
                                       corr_tide_mugal_list,
                                       corr_oceanload_mugal_list,
                                       corr_temp_mugal_list,
                                       corr_drift_mugal_list,
                                       gps_lat_deg_list,
                                       gps_lon_deg_list,
                                       gps_fix_quality_list,
                                       gps_satellites_list,
                                       gps_hdop_list,
                                       gps_h_m_list,
                                       duration_list,
                                       num_rejected_list,
                                       instr_height_m_list,
                                       gradient_mugalm_list,
                                       g_raw_mugal_list,
                                       setup_id_list,
                                       note_list)),
                              columns=cls._OBS_DF_COLUMN_NAMES)

        # Handle correction flags:
        if not all(item == correction_flags_list[0] for item in correction_flags_list):
            raise RuntimeError(f'The correction flags in the file {filename} are not equal for all observations!')
        drift_corr = correction_flags_list[0][0] == '1'
        temp_corr = correction_flags_list[0][1] == '1'
        tidal_corr = correction_flags_list[0][3] == '1'
        tilt_corr = correction_flags_list[0][4] == '1'

        return cls(obs_filename=filename,
                   obs_file_type='cg6_obs_file_solo',
                   survey_name=_HEADER_LINES['survey_name'][3],
                   serial_number=_HEADER_LINES['survey_name'][3],
                   created_datetime=created_datetime,
                   operator=_HEADER_LINES['operator'][3],
                   gcal1=_HEADER_LINES['gcal1'][3],
                   goff=_HEADER_LINES['goff'][3],
                   gref=_HEADER_LINES['gref'][3],
                   x_scale=_HEADER_LINES['x_scale'][3],
                   y_scale=_HEADER_LINES['y_scale'][3],
                   x_offset=_HEADER_LINES['x_offset'][3],
                   y_offset=_HEADER_LINES['y_offset'][3],
                   temp_corr_coeff=_HEADER_LINES['temp_corr_coeff'][3],
                   temp_sensor_scale=_HEADER_LINES['temp_sensor_scale'][3],
                   drift_rate=_HEADER_LINES['drift_rate'][3],
                   drift_ref_datetime=drift_ref_datetime,
                   tilt_corr=tilt_corr,
                   tidal_corr=tidal_corr,
                   ocean_loading_corr=False,
                   temp_corr=temp_corr,
                   drift_corr=drift_corr,
                   obs_df=obs_df,
                   # temp_sensor_offset=_HEADER_LINES['temp_sensor_offset'][3],  # Not available
                   firmware_version=_HEADER_LINES['firmware_version'][3]
                   # gravity_filter_type=_HEADER_LINES['gravity_filter_type'][3],  # Not available
                   # gravity_filter_period=_HEADER_LINES['gravity_filter_period'][3],  # Not available
                   # date_time_source=_HEADER_LINES['date_time_source'][3],  # Not available
                   # tidal_corr_model=_HEADER_LINES['tidal_corr_model'][3],  # Not available
                   # lynxlg_version=lynxlg_version  # Not available
                   )

    @classmethod
    def _check_obs_df_columns(cls, obs_df):
        """Check the columns of the `obs_df` DataFrame."""
        # Check for missing columns:
        missing_columns = set(cls._OBS_DF_COLUMN_NAMES) - set(obs_df.columns)
        if len(missing_columns):
            col_names = ', '.join(i for i in missing_columns)
            raise RuntimeError(f'Missing columns in CG6Survey.obs_df: {col_names}')

    @classmethod
    def _obs_df_reorder_columns(cls, obs_df):
        """Change order of columns of obs_df to the order specified in cls._OBS_DF_COLUMN_NAMES."""
        obs_df = obs_df[list(cls._OBS_DF_COLUMNS)]
        return obs_df

    @staticmethod
    def resolve_station_name(station_name_in):
        """Convert station name from Scintrex observation file
        (as Note) to the naming convention used in the output
        file (BEV conventions).

        Parameters
        ----------
        station_name_in : str
            Station name string as written to the observation file
            as note.

        Returns
        -------
        Corrected station name : str
        """
        station_name_in = station_name_in.upper()
        # Check first letter of name in order detect the station type:
        if station_name_in[0] == 'S':
            station_name_out = station_name_in
        elif station_name_in[0] == 'P':
            station_name_out = 'P  ' + station_name_in[1:]
        elif station_name_in[0] == 'T' and '.' in station_name_in:
            [str1_tmp, str2_tmp] = station_name_in.split('.')
            station_name_out = 'T{0:>4} {1:>3}'.format(str1_tmp[1:], str2_tmp)
        elif station_name_in[0] == 'N':
            station_name_out = station_name_in
        else:
            station_name_out = station_name_in.replace('.', '-')
        return station_name_out.upper()

    @staticmethod
    def get_dhb_dhf(dh_str):
        """Convert dhb and dhf from notes in the CG-5 observation file
        to an actual number. In the notes '.' is used instead of '-'.

        Parameters
        ----------
        dh_str : str
            Height difference [m] with '.' instead of '-'.

        Returns
        -------
        Height difference [m] : float
        """
        if dh_str.startswith('.'):
            dh_str = '-' + dh_str[1:]
        return float(dh_str)

    @property
    def number_of_observations(self):
        """Return the number of observations in the survey."""
        return len(self.obs_df)

    @property
    def number_of_setups(self):
        """Return the number of setups in the survey."""
        return len(set(self.obs_df['setup_id']))

    def __str__(self):
        if self.obs_df is None:
            return 'Empty CG-6 survey.'
        else:
            return f'Survey {self.survey_name} with {self.number_of_observations} observations in {self.number_of_setups} setups.'

    def plot_g_values(self, station_names=None):
        """Plot g-values of selected or all stations in the df.

        Notes
        -----
        This method requires matplotlib as optional dependency!

        Parameters
        ----------
        station_names : list of str, optional
            List of names of stations for which the observations will be plotted.
            The default value is None which implied that the data of all stations
            is plotted.

        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('Package matplotlib not installed.')
            return

            # Check if obs data is available first:
        if self.obs_df is not None:
            # Get list of station names and loop over them:
            if station_names is None:
                station_names = self.obs_df['station'].unique()
            fig, ax = plt.subplots()
            for station_name in station_names:
                x = self.obs_df[self.obs_df['station'] == station_name].ref_time
                y = self.obs_df[self.obs_df['station'] == station_name].g_corr_mugal
                ax.scatter(x, y, label=station_name)
            ax.legend()
            ax.grid(True)
            ax.set_xlabel('Time')
            ax.set_ylabel('g [µGal]')
            ax.set_title('survey: ' + self.survey_name)

            # Data cursor:
            fig.canvas.mpl_connect('pick_event', DataCursor(plt.gca()))
            plt.show()


# Run as standalone program:
if __name__ == "__main__":
    s_solo = CG6Survey.from_cg6solo_file('../../data/CG6_data/Heligrav_2023/CG-6_0265_HELIGRAV_2023-09-11.dat',
                                         300,
                                         verbose=True)
    print(s_solo)
    s_solo.plot_g_values()
else:
    pass
