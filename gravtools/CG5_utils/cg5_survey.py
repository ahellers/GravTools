"""Models for surveys with the Scintrex CG-5 gravity meters.

This module contains all models that are specifically required to
handle observation with the Scintrex CG-5 gravity meter.

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

import pandas as pd
import numpy as np
import re
import datetime as dt

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


class CG5SurveyParameters:
    """CG-5 Survey parameters.

    Scintrex CG-5 survey parameters from the 'CG-5 SURVEY' block in
    the observation file (text format).
    The class is initialized either by the :py:meth:`.__init__`,
    method by passing the attributes as keyword arguments directly,
    or by the class method :py:meth:`.populate_from_obs_file_string`
    that parses the content of a CG-5 observation file (txt).

    Attributes
    ----------
    survey_name : str
        Name of gravity survey.
    client : str
        Client of survey.
    operator : str
        Name of the survey's operator.
    long_deg : float
        Geographical longitude at begin of survey [°].
    lat_deg : float
        Geographical latitude at begin of survey [°].
    zone : str
        Timezone of all time records.
    date_time : datetime object
        Start epoch of the survey.
     """

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        *args
            Variable length argument list.
        **kwargs : dict
            Keyword arguments that are parsed to class attributes.
        """
        self.survey_name = kwargs.get('survey_name', '')  # string
        self.instrument_sn = kwargs.get('instrument_sn', '')  # string
        self.client = kwargs.get('client', '')  # string
        self.operator = kwargs.get('operator', '')  # string
        self.long_deg = kwargs.get('long_deg', np.nan)  # float
        self.lat_deg = kwargs.get('lat_deg', np.nan)  # float
        self.zone = kwargs.get('zone', '')  # string
        self.date_time = kwargs.get('date_time', None)  # datetime object (timezone aware)

    @classmethod
    def create_from_obs_file_string(cls, str_obs_file):
        """ Create instance of class :py:class:`.CG5SurveyParameters` by
        parsing a CG-5 observation file string.

        Parameters
        ----------
        str_obs_file : str
            Content of the CG-5 observation file (.txt).

        Returns
        -------
        Object: :py:class:`.CG5SurveyParameters`
            Initialized class instance.
        """

        # Parse observation file string:
        # expr = r'(\/\tCG-5 SURVEY\s*\n\/\s+Survey name:\s*(?P<survey_name>\S+)\s*\n\/\tInstrument S\/N:\s*(?P<instrument_sn>\S+)\s*\n\/\tClient:\s*(?P<client>\S+)\s*\n\/\tOperator:\s*(?P<operator>\S+)\s*\n\/\tDate:\s*(?P<date_year>\d{4})\/\s*(?P<date_month>\d{1,2})\/\s*(?P<date_day>\d{1,2})\s*\n\/\tTime:\s*(?P<time_hour>\d{2}):(?P<time_minu>\d{2}):(?P<time_sec>\d{2})\s*\n\/\tLONG:\s*(?P<long_num>\d*\.?\d*)\s+(?P<long_dir>[E|N|S|W])\s*\n\/\tLAT:\s*(?P<lat_num>\d*\.?\d*)\s+(?P<lat_dir>[E|N|S|W])\s*\n\/\tZONE:\s*(?P<zone>\d*)\s*\n\/\tGMT DIFF.:\s*(?P<gmt_diff>\d*\.?\d*))+\s*\n'
        expr = r'(\/\tCG-5 SURVEY\s*\n\/\s+Survey name:\s*(?P<survey_name>\S+)\s*\n' \
               r'\/\tInstrument S\/N:\s*(?P<instrument_sn>\S+)\s*\n' \
               r'\/\tClient:\s*(?P<client>\S+)\s*\n' \
               r'\/\tOperator:\s*(?P<operator>\S+)\s*\n' \
               r'\/\tDate:\s*(?P<date_year>\d{4})\/\s*(?P<date_month>\d{1,2})\/\s*(?P<date_day>\d{1,2})\s*\n' \
               r'\/\tTime:\s*(?P<time_hour>\d{2}):(?P<time_minu>\d{2}):(?P<time_sec>\d{2})\s*\n' \
               r'\/\tLONG:\s*(?P<long_num>\d*\.?\d*)\s+(?P<long_dir>[E|N|S|W])\s*\n' \
               r'\/\tLAT:\s*(?P<lat_num>\d*\.?\d*)\s+(?P<lat_dir>[E|N|S|W])\s*\n' \
               r'\/\tZONE:\s*(?P<zone>\d*)\s*\n' \
               r'\/\tGMT DIFF.:\s*(?P<gmt_diff>\d*\.?\d*))+\s*\n'

        # Read survey blocks from obs file string (only one block allowed!):
        survey_count = 0  # number of survey blocks in obs file string
        for survey_block in re.finditer(expr, str_obs_file):
            survey_dict = survey_block.groupdict()
            survey_count += 1

        if survey_count == 1:  # OK => Parse data in string:
            # Handle geographic locations:
            longitude_deg = float(survey_dict['long_num'])
            if survey_dict['long_dir'] == "W":
                longitude_deg = -longitude_deg
            latitude_deg = float(survey_dict['lat_num'])
            if survey_dict['lat_dir'] == "S":
                latitude_deg = -latitude_deg
            instrument_sn = survey_dict['instrument_sn']  # instrument serial number

            return cls(survey_name=survey_dict['survey_name'],
                       instrument_sn=instrument_sn,
                       client=survey_dict['client'],
                       operator=survey_dict['operator'],
                       long_deg=longitude_deg,
                       lat_deg=latitude_deg,
                       zone=survey_dict['zone'],
                       date_time=dt.datetime(int(survey_dict['date_year']),
                                             int(survey_dict['date_month']),
                                             int(survey_dict['date_day']),
                                             int(survey_dict['time_hour']),
                                             int(survey_dict['time_sec']),
                                             tzinfo=dt.timezone(dt.timedelta(hours=float(survey_dict['gmt_diff'])))
                                             )
                       )
        elif survey_count == 0:  # Not available
            return cls()  # Initialize with default values
        else:  # More than 1 block found => Error!
            raise InvaliFileContentError('{} "CG-5 SURVEY" blocks found in observation file, '
                                         'but only one expected.'.format(survey_count))

        # Error Msg, wenn der Block mehr als einmal gefunden wird.


class CG5SetupParameters:
    """CG-5 Survey parameters.

    Scintrex CG-5 Setup parameters from the 'CG-5 SETUP PARAMETERS'
    block in the observation file (text format).
    The class is initialized either by the :py:meth:`.__init__`,
    method by passing the attributes as keyword arguments directly,
    or by the class method :py:meth:`.populate_from_obs_file_string`
    that parses the content of a CG-5 observation file (txt).

    Attributes
    ----------
    gcal1 : float
        Calibration factor GCAL1.
    tiltxs : float
        XXXXXXXX
    tiltys : float
       XXXXXXXX
    tiltxo : float
        XXXXXXXX
    tiltyo : float
        XXXXXXXX
    tempco : float
        XXXXXXXX
    drift : float
        Linear drift factor (long-term drift).
    drift_date_time_start : datetime object (TZ aware)
        Start epoch for the determination of th linear long-term drift.
    """

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        *args
            Variable length argument list.
        **kwargs : dict
            Keyword arguments that are parsed to class attributes.
        """
        self.gcal1 = kwargs.get('gcal1', np.nan)  # float
        self.tiltxs = kwargs.get('tiltxs', np.nan)  # float
        self.tiltys = kwargs.get('tiltys', np.nan)  # string
        self.tiltxo = kwargs.get('tiltxo', np.nan)  # float
        self.tiltyo = kwargs.get('tiltyo', np.nan)  # float
        self.tempco = kwargs.get('tempco', np.nan)  # float
        self.drift = kwargs.get('drift', np.nan)  # float
        self.drift_date_time_start = kwargs.get('drift_date_time_start', None)  # datetime object (timezone aware)

    @classmethod
    def create_from_obs_file_string(cls, str_obs_file):
        """ Create instance of class :py:class:`.CG5SetupParameters` by
        parsing a CG-5 observation file string.

        Parameters
        ----------
        str_obs_file : str
            Content of the CG-5 observation file (.txt).

        Returns
        -------
        Object: :py:class:`.CG5SetupParameters`
            Initializes class instance.
        """

        # Parse observation file string:
        expr = r"\/\tCG-5 SETUP PARAMETERS\s*\n\/\s+Gref:\s*(?P<gref>\S+)\s*\n\/\s+Gcal1:\s*(" \
               r"?P<gcal1>\S+)\s*\n\/\s+TiltxS:\s*(?P<tiltxs>\S+)\s*\n\/\s+TiltyS:\s*(" \
               r"?P<tiltys>\S+)\s*\n\/\s+TiltxO:\s*(?P<tiltxo>\S+)\s*\n\/\s+TiltyO:\s*(" \
               r"?P<tiltyo>\S+)\s*\n\/\s+Tempco:\s*(?P<tempco>\S+)\s*\n\/\s+Drift:\s*(" \
               r"?P<drift>\S+)\s*\n\/\s+DriftTime Start:\s*(?P<drift_time_start>\S+)\s*\n\/\s+DriftDate Start:\s*(" \
               r"?P<drift_date_start>\S+)\s*\n"

        # Read setup parameters blocks from obs file string (only one block allowed!):
        block_count = 0  # number of survey blocks in obs file string
        for survey_block in re.finditer(expr, str_obs_file):
            setup_dict = survey_block.groupdict()
            block_count += 1

        if block_count == 1:  # OK => Parse data in string:

            return cls(gref=float(setup_dict['gref']),
                       gcal1=float(setup_dict['gcal1']),
                       tiltxs=float(setup_dict['tiltxs']),
                       tiltys=float(setup_dict['tiltys']),
                       tiltxo=float(setup_dict['tiltxo']),
                       tiltyo=float(setup_dict['tiltyo']),
                       tempco=float(setup_dict['tempco']),
                       drift=float(setup_dict['drift']),
                       drift_date_time_start=dt.datetime.strptime(
                           setup_dict["drift_date_start"] + setup_dict["drift_time_start"],
                           "%Y/%m/%d%H:%M:%S")
                       )
        elif block_count == 0:  # Not available
            return cls()  # Initialize with default values
        else:  # More than 1 block found => Error!
            raise InvaliFileContentError('{} "CG-5 SETUP PARAMETERS" in observation file found, but maximum one '
                                         'expected.'.format(block_count))

    # def is_valid(self) -> bool:
    #     """???"""
    #     # TODO
    #     return True


class CG5OptionsParameters:
    """CG-5 instrumental options (filters, corrections, output).

    Scintrex CG-5 options from the 'CG-5 OPTIONS' block in
    the observation file (text format).
    The class is initialized either by the :py:meth:`.__init__`,
    method by passing the attributes as keyword arguments directly,
    or by the class method :py:meth:`.populate_from_obs_file_string`
    that parses the content of a CG-5 observation file (txt).

    Attributes
    ----------
    tide_correction : bool
        Tide correction (on/off).
    cont_tilt : bool
        Continuous tilt correction (on/off).
    auto_rejection : bool
       Auto rejection of outliers (on/off).
    terrain_correction : bool
        Terrain correction (on/off).
    seismic_filter : bool
        Seismic filter (on/off).
    raw_data : bool
        Raw data output (on/off).
    """

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        *args
            Variable length argument list.
        **kwargs : dict
            Keyword arguments that are parsed to class attributes.
        """
        self.tide_correction = kwargs.get('tide_correction', None)  # bool
        self.cont_tilt = kwargs.get('cont_tilt', None)  # bool
        self.auto_rejection = kwargs.get('auto_rejection', None)  # bool
        self.terrain_correction = kwargs.get('terrain_correction', None)  # bool
        self.seismic_filter = kwargs.get('seismic_filter', None)  # bool
        self.raw_data = kwargs.get('raw_data', None)  # bool
        #
        # self.tide_correction = kwargs.get('tide_correction', False)  # bool
        # self.cont_tilt = kwargs.get('cont_tilt', False)  # bool
        # self.auto_rejection = kwargs.get('auto_rejection', False)  # bool
        # self.terrain_correction = kwargs.get('terrain_correction', False)  # bool
        # self.seismic_filter = kwargs.get('seismic_filter', False)  # bool
        # self.raw_data = kwargs.get('raw_data', False)  # bool

    @classmethod
    def create_from_obs_file_string(cls, str_obs_file):
        """ Create instance of class :py:class:`.CG5OptionsParameters` by
        parsing a CG-5 observation file string.

        Parameters
        ----------
        str_obs_file : str
            Content of the CG-5 observation file (.txt).

        Returns
        -------
        Object: :py:class:`.CG5OptionsParameters`
            Initializes class instance.
        """

        # Parse observation file string:
        expr = r"\/\tCG-5 OPTIONS\s*\n\/\s+Tide Correction:\s*(?P<tide_correction>\S+)\s*\n\/\s+Cont. Tilt:\s*(" \
               r"?P<cont_tilt>\S+)\s*\n\/\s+Auto Rejection:\s*(?P<auto_rejection>\S+)\s*\n\/\s+Terrain Corr.:\s*(" \
               r"?P<terrain_correction>\S+)\s*\n\/\s+Seismic Filter:\s*(?P<seismic_filter>\S+)\s*\n\/\s+Raw Data:\s*(" \
               r"?P<raw_data>\S+)\s*\n"

        # Read setup parameters blocks from obs file string (only one block allowed!):
        block_count = 0  # number of survey blocks in obs file string
        for options_block in re.finditer(expr, str_obs_file):
            options_dict = options_block.groupdict()
            block_count += 1

        if block_count == 1:  # OK => Parse data in string:

            return cls(tide_correction=options_dict['tide_correction'] == "YES",  # bool
                       cont_tilt=options_dict['cont_tilt'] == "YES",  # bool
                       auto_rejection=options_dict['auto_rejection'] == "YES",  # bool
                       terrain_correction=options_dict['terrain_correction'] == "YES",  # bool
                       seismic_filter=options_dict['seismic_filter'] == "YES",  # bool
                       raw_data=options_dict['raw_data'] == "YES",
                       )  # bool

        elif block_count == 0:  # Not available
            return cls()  # Initialize with default values
        else:  # More than 1 block found => Error!
            raise InvaliFileContentError('{} "CG-5 OPTIONS" in observation file found, but maximum one '
                                         'expected.'.format(block_count))

    # def is_valid(self) -> bool:
    #     """???"""
    #     # TODO
    #     return True


class CG5Survey:
    """CG-5 survey data.

    Class instances may contain the information available in the
    following sections of Scintrex CG-5 observation files (txt format):

    - Survey Parameter block (as instance of :py:obj:`.CG5SurveyParameters`)
    - Setup block (as instance of :py:obj:`.CG5SetupParameters`)
    - Options block (as instance of :py:obj:`.CG5OptionsParameters`)
    - Observations (as pandas dataframe)

    If the class is initialized without setting the observations file
    attribute (`obs_filename`), no observation file is load and the
    object is initialized empty.

    Attributes
    ----------
    obs_filename : str
        Name (and path) to CG-5 observation file (txt format)
    survey_parameters : :py:class:`.CG5SurveyParameters`
        Survey Parameter of the Parameters block in the observation file.
    setup_parameters :  :py:class:`.CG5SetupParameters`
        Setup Parameter of the Setup block in the observation file.
    options : :py:class:`.CG5OptionsParameters`
        Instrumental options from the Options block in the observation file.
    obs_df : pandas data frame
       Contains the actual observation data records with the colums defined in `self._OBS_DF_COLUMN_NAMES`.
    """

    # Column names of the dataframe containing tha actual observation data:
    _OBS_DF_COLUMN_NAMES = ('lat_deg',  # Latitude [deg] :
                            'lon_deg',  # Longitude [deg]
                            'alt_m',  # Altitude [m]
                            'g_mgal',  # Determined gravity (corrected) [mGal]
                            'sd_mgal',  # Standard deviation of determined gravity [mGal]
                            'tiltx',
                            'tilty',
                            'temp',
                            'tide',  # Tidal correction determined by the CG-5 [mGal]
                            'duration_sec',  # Duration of the current setup [sec]
                            'rej',  # Number of rejected single measurements
                            'time_str',  # Reference time = mid of setup with duration `duration_sec`) (dropped later)
                            'dec_time_date',
                            'terrain',  # Terrain correction [??]
                            'date',  # Date (dropped later)
                            'station_name',  # Station name : str
                            'dhf_m',  # Distance between instrument top and physical reference point [m]
                            'dhb_m',  # Distance between instrument top and ground [m]
                            'atm_pres_hpa',  # Measured atmospheric pressure [hPa]
                            'setup_id',  # Unique ID of this observation (=setup)
                            )
                            # obs_epoch : datetime object (added to df later)

    # Rename columns: df.rename(columns = {'$b':'B'}, inplace = True)

    # Non-numeric columns in the observation dataframe:
    _OBS_DF_NON_NUMERIC_COLUMNS = ['station_name', 'date', 'time_str']

    def __init__(self,
                 obs_filename='',
                 survey_parameters=CG5SurveyParameters(),
                 setup_parameters=CG5SetupParameters(),
                 options=CG5OptionsParameters()
                 ):
        """
        Parameters
        ----------
        obs_filename : str, optional
            Name (and path) to CG-5 observation file (txt format).
        survey_parameters : :py:class:`.CG5SurveyParameters`, optional
            Survey Parameter of the Parameters block in the observation file.
        setup_parameters :  :py:class:`.CG5SetupParameters`, optional
            Setup Parameter of the Setup block in the observation file.
        options : :py:class:`.CG5OptionsParameters`, optional
            Instrumental options from the Options block in the observation file.
        """
        self.obs_filename = obs_filename
        assert isinstance(survey_parameters, CG5SurveyParameters), \
            "survey_parameters is not an instance of CG5SurveyParameters"
        self.survey_parameters = survey_parameters
        assert isinstance(setup_parameters, CG5SetupParameters), \
            "setup_parameters is not an instance of CG5SetupParameters"
        self.setup_parameters = setup_parameters
        assert isinstance(options, CG5OptionsParameters), \
            "options is not an instance of CG5OptionsParameters"
        self.options = options

        # Read observation file, if a valid filename is available and valid. Otherwise, initialize obs_df as None.
        if self.obs_filename:
            self.read_obs_file(obs_filename)
        else:
            self.obs_df = None  # Initialize as None

    def __str__(self):
        if self.obs_df is None:
            return 'Empty CG-5 Survey.'
        else:
            if not self.survey_parameters.survey_name:
                return 'Unnamed CG-5 Survey with {} observations (file: {}).'.format(len(self.obs_df),
                                                                                     self.obs_filename.split('/')[-1])
            else:
                return 'CG-5 Survey "{}" with {} observations (file: {}).'.format(self.survey_parameters.survey_name,
                                                                                  len(self.obs_df),
                                                                                  self.obs_filename.split('/')[-1])

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

    def read_obs_file(self, obs_filename):
        """Read CG-5 observation file (txt) and populate the object.

        Notes
        -----
        Ignore comment lines that start with "#".

        Parameters
        ----------
        obs_filename : str
            Name (and path) to CG-5 observation file (txt format).
        """
        COMMENT_MARKER = '#'

        self.obs_filename = obs_filename
        # with open(self.obs_filename, 'r') as content_file:
        #     str_obs_file = content_file.read()

        # Read in file and ignore comment lines:
        file_handle = open(self.obs_filename, 'r')
        lines = []
        for line in file_handle:
            line_tmp = line.strip()
            if not line_tmp.startswith(COMMENT_MARKER):
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
            str_obs_file = str_obs_file[0:str_idx+2]

        # ### Match blocks with regex ###
        # CG-5 SURVEY block:
        self.survey_parameters = CG5SurveyParameters.create_from_obs_file_string(str_obs_file)
        self.setup_parameters = CG5SetupParameters.create_from_obs_file_string(str_obs_file)
        self.options = CG5OptionsParameters.create_from_obs_file_string(str_obs_file)

        # Get Observations
        # ### 3 possibilities: ###
        # Initialize empty dataframe and append observation blocks at each station
        # Warning: Better performance, when preparing the data as list (appending) and then converting to df at once.
        #  See: https://stackoverflow.com/questions/13784192/creating-an-empty-pandas-dataframe-then-filling-it

        obs_list = []  # Collect all obs data in this list and then convert to pd dataframe.

        # 1.) Station name & dbh=dhf
        expr = '\/\tNote:   \t(?P<station_name>\S+)\s+(?P<dh_cm>-?[.0-9]+)\s*[\r?\n](?P<obs_data>(?:\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*[\r?\n])+)'
        for obs_block in re.finditer(expr, str_obs_file):
            obs_dict = obs_block.groupdict()
            station_name = self.resolve_station_name(obs_dict['station_name'])
            dhf_m = float(obs_dict['dh_cm']) * 1e-2
            dhb_m = float(obs_dict['dh_cm']) * 1e-2
            atm_pres_hpa = None
            lines = obs_dict['obs_data'].splitlines()
            # Create unique ID (= UNIX timestamp of first observation) for each setup on a station:
            #  - To distinguish multiple setups (with multiple observations each) on multiple stations
            time_str = lines[0].split()[-1] + ' ' + lines[0].split()[11]
            setup_id = int(dt.datetime.timestamp(dt.datetime.strptime(time_str, "%Y/%m/%d %H:%M:%S")))

            for line in lines:
                line_items = line.split()
                line_items.append(station_name)
                line_items.append(dhf_m)
                line_items.append(dhb_m)
                line_items.append(atm_pres_hpa)
                line_items.append(setup_id)
                obs_list.append(line_items)

        # 2.) Station name & dbh=dhf & pressure
        expr = '\/\tNote:   \t(?P<station_name>\S+)\s+(?P<dh_cm>-?[.0-9]+)\s*[\r?\n](?P<obs_data>(?:\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*[\r?\n])+)\/\tNote:   \t(?P<pres>[0-9]{3,4}[.]{0,1}[0-9]*)'
        for obs_block in re.finditer(expr, str_obs_file):
            obs_dict = obs_block.groupdict()
            station_name = self.resolve_station_name(obs_dict['station_name'])
            dhf_m = float(obs_dict['dh_cm']) * 1e-2
            dhb_m = float(obs_dict['dh_cm']) * 1e-2
            atm_pres_hpa = float(obs_dict['pres'])
            lines = obs_dict['obs_data'].splitlines()
            # Create unique ID (= UNIX timestamp of first observation) for each setup on a station:
            #  - To distinguish multiple setups (with multiple observations each) on multiple stations
            time_str = lines[0].split()[-1] + ' ' + lines[0].split()[11]
            setup_id = int(dt.datetime.timestamp(dt.datetime.strptime(time_str, "%Y/%m/%d %H:%M:%S")))

            for line in lines:
                line_items = line.split()
                line_items.append(station_name)
                line_items.append(dhf_m)
                line_items.append(dhb_m)
                line_items.append(atm_pres_hpa)
                line_items.append(setup_id)
                obs_list.append(line_items)

        # 3.) Station name & dhb & dhf
        expr = '\/\tNote:   \t(?P<station_name>\S+)\s+(?P<dhb_cm>-?[.0-9]+)\s+(?P<dhf_cm>-?[.0-9]+)\s*[\r?\n](?P<obs_data>(?:\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*[\r?\n])+)'
        for obs_block in re.finditer(expr, str_obs_file):
            obs_dict = obs_block.groupdict()

            station_name = self.resolve_station_name(obs_dict['station_name'])
            dhf_m = float(obs_dict['dhf_cm']) * 1e-2
            dhb_m = float(obs_dict['dhb_cm']) * 1e-2
            atm_pres_hpa = None
            lines = obs_dict['obs_data'].splitlines()
            # Create unique ID (= UNIX timestamp of first observation) for each setup on a station:
            #  - To distinguish multiple setups (with multiple observations each) on multiple stations
            time_str = lines[0].split()[-1] + ' ' + lines[0].split()[11]
            setup_id = int(dt.datetime.timestamp(dt.datetime.strptime(time_str, "%Y/%m/%d %H:%M:%S")))

            for line in lines:
                line_items = line.split()
                line_items.append(station_name)
                line_items.append(dhf_m)
                line_items.append(dhb_m)
                line_items.append(atm_pres_hpa)
                line_items.append(setup_id)
                obs_list.append(line_items)

        # 4.) Station name & dhb & dhf & pressure
        # expr = "\/\s+Note:\s+(?P<station_name>\S+)\s+(?P<dhb_cm>\S+)\s+(?P<dhf_cm>\S+)\s*\n(?P<obs_data>(?:\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*[\r\n])+)"
        expr = '\/\tNote:   \t(?P<station_name>\S+)\s+(?P<dhb_cm>-?[.0-9]+)\s+(?P<dhf_cm>-?[.0-9]+)\s*[\r?\n](?P<obs_data>(?:\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*[\r?\n])+)\/\tNote:   \t(?P<pres>[0-9]{3,4}[.]{0,1}[0-9]*)'
        for obs_block in re.finditer(expr, str_obs_file):
            obs_dict = obs_block.groupdict()

            station_name = self.resolve_station_name(obs_dict['station_name'])
            dhf_m = float(obs_dict['dhf_cm']) * 1e-2
            dhb_m = float(obs_dict['dhb_cm']) * 1e-2
            atm_pres_hpa = float(obs_dict['pres'])
            lines = obs_dict['obs_data'].splitlines()
            # Create unique ID (= UNIX timestamp of first observation) for each setup on a station:
            #  - To distinguish multiple setups (with multiple observations each) on multiple stations
            time_str = lines[0].split()[-1] + ' ' + lines[0].split()[11]
            setup_id = int(dt.datetime.timestamp(dt.datetime.strptime(time_str, "%Y/%m/%d %H:%M:%S")))

            for line in lines:
                line_items = line.split()
                line_items.append(station_name)
                line_items.append(dhf_m)
                line_items.append(dhb_m)
                line_items.append(atm_pres_hpa)
                line_items.append(setup_id)
                obs_list.append(line_items)

        # Create pandas dataframe of prepared list:
        self.obs_df = pd.DataFrame(obs_list, columns=self._OBS_DF_COLUMN_NAMES)

        # Remove duplicates entries that were matched with and without pressure:
        setup_ids = self.obs_df['setup_id'].unique().tolist()
        setup_ids_diplicates = []
        for setup_id in setup_ids:
            tmp_filter = self.obs_df['setup_id'] == setup_id
            if (~self.obs_df.loc[tmp_filter, 'atm_pres_hpa'].isna()).any():  # Entries with pressure found
                setup_ids_diplicates.append(setup_id)

        if setup_ids_diplicates:
            tmp_filter = ~(self.obs_df['atm_pres_hpa'].isna() & self.obs_df['setup_id'].isin(setup_ids_diplicates))
            self.obs_df = self.obs_df.loc[tmp_filter].copy(deep=True)

        # Convert numeric columns to numeric dtypes:
        cols = self.obs_df.columns.drop(self._OBS_DF_NON_NUMERIC_COLUMNS)
        self.obs_df[cols] = self.obs_df[cols].apply(pd.to_numeric, errors='raise')
        # Sort observations by time and date and reset index:
        self.obs_df.sort_values(by='dec_time_date', inplace=True, ignore_index=True)
        # Convert date and time to datetime objects (aware, if UTC offset is available):
        if self.survey_parameters.date_time is not None:
            self.obs_df['obs_epoch'] = pd.to_datetime(
                self.obs_df['date'] + ' ' + self.obs_df['time_str'], format='%Y/%m/%d %H:%M:%S')
            self.obs_df['obs_epoch'] = self.obs_df['obs_epoch'].dt.tz_localize('UTC')  # Set timezone = UTC
            self.obs_df['obs_epoch'] = self.obs_df['obs_epoch'] + pd.Timedelta(self.survey_parameters.date_time.utcoffset())
        else:  # tz unaware time
            self.obs_df['obs_epoch'] = pd.to_datetime(
                self.obs_df['date'] + ' ' + self.obs_df['time_str'], format='%Y/%m/%d %H:%M:%S')

        self.obs_df.drop(columns=['time_str', 'date'], inplace=True)  # Drop columns that are not required any more

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
        # Check if obs data is available first:
        if self.obs_df is not None:
            # Get list of station names and loop over them:
            if station_names is None:
                station_names = self.obs_df['station_name'].unique()
            fig, ax = plt.subplots()
            for station_name in station_names:
                x = self.obs_df[self.obs_df['station_name'] == station_name].obs_epoch
                y = self.obs_df[self.obs_df['station_name'] == station_name].g_mgal * 1e3  # µGal
                ax.scatter(x, y, label=station_name)
            ax.legend()
            ax.grid(True)
            ax.set_xlabel('Time')
            ax.set_ylabel('g [µGal]')
            ax.set_title('survey: ' + self.survey_parameters.survey_name)

            # Data cursor:
            fig.canvas.mpl_connect('pick_event', DataCursor(plt.gca()))

            plt.show()


# Run as standalone program:
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    path = settings.PATH_OBS_FILE_CG5 + settings.NAME_OBS_FILE_CG5
    s1 = CG5Survey()
    s1.read_obs_file(path)
    # s1.plot_g_values(['1-164-04', '1-164-12', '1-164-11'])
    # s1.plot_g_values(['TEST'])
    s1.plot_g_values()
else:
    pass
