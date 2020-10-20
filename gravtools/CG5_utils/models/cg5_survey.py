import pandas as pd
import numpy as np
import re
import sys
import datetime as dt
import matplotlib.pyplot as plt
from io import StringIO  # Python 3.x required

from gravtools.models.exceptions import InvaliFileContentError
from gravtools import settings


class CG5SurveyParameters:
    """CG-5 CG5Survey parameters from the 'CG-5 SURVEY' block in the observation file."""

    def __init__(self, *args, **kwargs):
        """Initialize object."""
        self.survey_name = kwargs.get('survey_name', '')  # string
        self.client = kwargs.get('client', '')  # string
        self.operator = kwargs.get('operator', '')  # string
        self.long_deg = kwargs.get('long_deg', np.nan)  # np.float
        self.lat_deg = kwargs.get('lat_deg', np.nan)  # np.float
        self.zone = kwargs.get('zone', '')  # string
        self.date_time = kwargs.get('date_time', None)  # datetime object (timezone aware)

    @classmethod
    def populate_from_obs_file_string(cls, str_obs_file):
        """ Create class instance by parsing a CG-5 observation file string."""

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
            longitude_deg = np.float(survey_dict['long_num'])
            if survey_dict['long_dir'] == "W":
                longitude_deg = -longitude_deg
            latitude_deg = np.float(survey_dict['lat_num'])
            if survey_dict['lat_dir'] == "S":
                latitude_deg = -latitude_deg

            return cls(survey_name=survey_dict['survey_name'],
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
                                             ))
        elif survey_count == 0:  # Not available
            return cls()  # Initialize with default values
        else:  # More than 1 block found => Error!
            raise InvaliFileContentError('{} "CG-5 SURVEY" in observation file found, '
                                         'but maximum one expected.'.format(survey_count))

        # Error Msg, wenn der Block mehr als einmal gefunden wird.


class CG5SetupParameters:
    """CG-5 Setup parameters from the 'CG-5 SETUP PARAMETERS' block in the observation file."""

    def __init__(self, *args, **kwargs):
        """Initialize object."""
        self.gcal1 = kwargs.get('gcal1', np.nan)  # np.float
        self.tiltxs = kwargs.get('tiltxs', np.nan)  # np.float
        self.tiltys = kwargs.get('tiltys', np.nan)  # string
        self.tiltxo = kwargs.get('tiltxo', np.nan)  # np.float
        self.tiltyo = kwargs.get('tiltyo', np.nan)  # np.float
        self.tempco = kwargs.get('tempco', np.nan)  # np.float
        self.drift = kwargs.get('drift', np.nan)  # np.float
        self.drift_date_time_start = kwargs.get('drift_date_time_start', None)  # datetime object (timezone aware)

    @classmethod
    def populate_from_obs_file_string(cls, str_obs_file):
        """ Populates class instance by parsing a CG-5 observation file string."""

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

            return cls(gref=np.float(setup_dict['gref']),
                       gcal1=np.float(setup_dict['gcal1']),
                       tiltxs=np.float(setup_dict['tiltxs']),
                       tiltys=np.float(setup_dict['tiltys']),
                       tiltxo=np.float(setup_dict['tiltxo']),
                       tiltyo=np.float(setup_dict['tiltyo']),
                       tempco=np.float(setup_dict['tempco']),
                       drift=np.float(setup_dict['drift']),
                       drift_date_time_start=dt.datetime.strptime(
                           setup_dict["drift_date_start"] + setup_dict["drift_time_start"],
                           "%Y/%m/%d%H:%M:%S"))
        elif block_count == 0:  # Not available
            return cls()  # Initialize with default values
        else:  # More than 1 block found => Error!
            raise InvaliFileContentError('{} "CG-5 SETUP PARAMETERS" in observation file found, but maximum one '
                                         'expected.'.format(block_count))

    def is_valid(self) -> bool:
        """???"""
        # TODO
        return True


class CG5OptionsParameters:
    """CG-5 options from the 'CG-5 OPTIONS' block in the observation file."""

    def __init__(self, *args, **kwargs):
        """Initialize object."""
        self.tide_correction = kwargs.get('tide_correction', False)  # bool
        self.cont_tilt = kwargs.get('cont_tilt', False)  # bool
        self.auto_rejection = kwargs.get('auto_rejection', False)  # bool
        self.terrain_correction = kwargs.get('terrain_correction', False)  # bool
        self.seismic_filter = kwargs.get('seismic_filter', False)  # bool
        self.raw_data = kwargs.get('raw_data', False)  # bool

    @classmethod
    def populate_from_obs_file_string(cls, str_obs_file):
        """ Create class instance by parsing a CG-5 observation file string."""

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
                       raw_data=options_dict['raw_data'] == "YES")  # bool

        elif block_count == 0:  # Not available
            return cls()  # Initialize with default values
        else:  # More than 1 block found => Error!
            raise InvaliFileContentError('{} "CG-5 OPTIONS" in observation file found, but maximum one '
                                         'expected.'.format(block_count))

    def is_valid(self) -> bool:
        """???"""
        # TODO
        return True


class CG5Survey:
    """CG5 CG5Survey data object."""

    # PARAM_ATTRIBUTES = ["survey_name", ]

    def __init__(self, survey_parameters=CG5SurveyParameters(),
                 setup_parameters=CG5SetupParameters(),
                 options=CG5OptionsParameters()):
        assert isinstance(survey_parameters, CG5SurveyParameters), \
            "survey_parameters is not an instance of CG5SurveyParameters"
        self.survey_parameters = survey_parameters
        assert isinstance(setup_parameters, CG5SetupParameters), \
            "setup_parameters is not an instance of CG5SetupParameters"
        self.setup_parameters = setup_parameters
        assert isinstance(options, CG5OptionsParameters), \
            "options is not an instance of CG5OptionsParameters"
        self.options = options

        self.obs_df = None  # Initilaize as None

        # parameters = kwargs.get('parameters', None)
        # assert isinstance(self.parameters, CG5SurveyParameters)
        # assert parameters.is_valid()
        # self.parameters = parameters

    # def set_params(self, param_dict):
    #     for param in param_dict:
    #         if param in self.PARAM_ATTRIBUTES:
    #             self.__setattr__("param_" + param, param_dict[param])

    def __str__(self):
        if self.obs_df is None:
            return 'Empty CG-5 Survey.'
        else:
            if not self.survey_parameters.survey_name:
                return 'Unnamed CG-5 Survey with {} observations.'.format(len(self.obs_df))
            else:
                return 'CG-5 Survey "{}" with {} observations.'.format(self.survey_parameters.survey_name,
                                                                       len(self.obs_df))

    @staticmethod
    def resolve_station_name(station_name_in):
        """Convert station name from Scintrex observation file (as Note) to the naming convention used in the output
        file."""
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
        """Convert dhb and dhf from Notes in the CG-5 observation file to an actual number."""
        if dh_str.startswith('.'):
            dh_str = '-' + dh_str[1:]
        return np.float(dh_str)

    def read_obs_file(self, file_path):
        """Read CG5 observation file an populate object."""
        self.file_path = file_path
        with open(self.file_path, 'r') as content_file:
            str_obs_file = content_file.read()

        str_obs_file += '\n'  # Last character of string has to be a \n so that regex works correctly!

        # ### Match blocks with regex ###
        # CG-5 SURVEY block:
        self.survey_parameters = CG5SurveyParameters.populate_from_obs_file_string(str_obs_file)
        self.setup_parameters = CG5SetupParameters.populate_from_obs_file_string(str_obs_file)
        self.options = CG5OptionsParameters.populate_from_obs_file_string(str_obs_file)

        # Get Observations
        # ### 3 possibilities: ###
        # Initialize empty dataframe and append observation blocks at each station
        # Warning: Better performance, when preparing the data as list (appending) and then converting to df at once.
        #  See: https://stackoverflow.com/questions/13784192/creating-an-empty-pandas-dataframe-then-filling-it

        obs_list = []  # Collect all obs data in this list and then convert to pd dataframe.
        column_names = ['lat_deg', 'lon_deg', 'alt_m', 'g_mgal', 'sd_mgal', 'tiltx', 'tilty', 'temp',
                        'tide', 'duration_sec', 'rej', 'time_str', 'dec_time_date', 'terrain', 'date',
                        'station_name', 'dhf_m', 'dhb_m', 'setup_id']
        non_numeric_columns = ['station_name', 'date', 'time_str']

        # 1.) Only Station name
        expr = "\/\s+Note:\s+(?P<station_name>\S+)\s*\n(?P<obs_data>(?:\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*[\r\n])+)"
        dhf_m = 0.0
        dhb_m = 0.0
        for obs_block in re.finditer(expr, str_obs_file):
            obs_dict = obs_block.groupdict()
            station_name = self.resolve_station_name(obs_dict['station_name'])
            lines = obs_dict['obs_data'].splitlines()
            # Create unique ID (= UNIX timestamp of first observation) for each setup on a station:
            #  - To distinguish multiple setups (with multiple observations each) on multiple stations
            time_str = lines[0].split()[-1] + ' ' + lines[0].split()[11]
            setup_id = dt.datetime.timestamp(dt.datetime.strptime(time_str, "%Y/%m/%d %H:%M:%S"))

            for line in lines:
                line_items = line.split()
                line_items.append(station_name)
                line_items.append(dhf_m)
                line_items.append(dhb_m)
                line_items.append(setup_id)
                obs_list.append(line_items)

        # 2.) Station name & dbh=dhf
        expr = "\/\s+Note:\s+(?P<station_name>\S+)\s+(?P<dh_cm>\S+)\s*\n(?P<obs_data>(?:\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*[\r\n])+)"
        for obs_block in re.finditer(expr, str_obs_file):
            obs_dict = obs_block.groupdict()
            station_name = self.resolve_station_name(obs_dict['station_name'])
            dhf_m = np.float(obs_dict['dh_cm'])
            dhb_m = np.float(obs_dict['dh_cm'])
            lines = obs_dict['obs_data'].splitlines()
            # Create unique ID (= UNIX timestamp of first observation) for each setup on a station:
            #  - To distinguish multiple setups (with multiple observations each) on multiple stations
            time_str = lines[0].split()[-1]+' '+lines[0].split()[11]
            setup_id = dt.datetime.timestamp(dt.datetime.strptime(time_str, "%Y/%m/%d %H:%M:%S"))

            for line in lines:
                line_items = line.split()
                line_items.append(station_name)
                line_items.append(dhf_m)
                line_items.append(dhb_m)
                line_items.append(setup_id)
                obs_list.append(line_items)

        # 3.) Station name & dhb & dhf
        expr = "\/\s+Note:\s+(?P<station_name>\S+)\s+(?P<dhb_cm>\S+)\s+(?P<dhf_cm>\S+)\s*\n(?P<obs_data>(?:\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*[\r\n])+)"
        for obs_block in re.finditer(expr, str_obs_file):
            obs_dict = obs_block.groupdict()

            station_name = self.resolve_station_name(obs_dict['station_name'])
            dhf_m = np.float(obs_dict['dhf_cm'])
            dhb_m = np.float(obs_dict['dhb_cm'])
            lines = obs_dict['obs_data'].splitlines()
            # Create unique ID (= UNIX timestamp of first observation) for each setup on a station:
            #  - To distinguish multiple setups (with multiple observations each) on multiple stations
            time_str = lines[0].split()[-1] + ' ' + lines[0].split()[11]
            setup_id = dt.datetime.timestamp(dt.datetime.strptime(time_str, "%Y/%m/%d %H:%M:%S"))

            for line in lines:
                line_items = line.split()
                line_items.append(station_name)
                line_items.append(dhf_m)
                line_items.append(dhb_m)
                line_items.append(setup_id)
                obs_list.append(line_items)

        # Create pandas dataframe of prepared list:
        self.obs_df = pd.DataFrame(obs_list, columns=column_names)

        # Convert numeric columns to numeric dtypes:
        cols = self.obs_df.columns.drop(non_numeric_columns)
        self.obs_df[cols] = self.obs_df[cols].apply(pd.to_numeric, errors='raise')
        # Sort observations by time and date:
        self.obs_df.sort_values(by='dec_time_date')
        # Reset index after sorting by observation time:
        self.obs_df.reset_index(drop=True)
        # Convert date and time to datetime objects (aware, if UTC offset is available):

        if self.survey_parameters.date_time is not None:
            self.obs_df['obs_epoch'] = pd.to_datetime(
                self.obs_df['date'] + ' ' + self.obs_df['time_str'] + ' ' + self.survey_parameters.date_time.tzname(),
                format='%Y/%m/%d %H:%M:%S %Z')
        else:  # tz unaware time
            self.obs_df['obs_epoch'] = pd.to_datetime(
                self.obs_df['date'] + ' ' + self.obs_df['time_str'], format='%Y/%m/%d %H:%M:%S')
        pass
        self.obs_df.drop(columns=['time_str', 'date'], inplace=True)  # Drop columns that are not required any more

    def plot_g_values(self, station_names=None):
        """Plot g-values of selected or all stations in the df.

        Attributes:
            station_names - List of names of the stations that should be plotted.
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
            plt.show()


# Run as standalone program:
if __name__ == "__main__":

    # path = '/home/heller/pyProjects/gravtools/data/20200907_test.TXT'
    path = settings.PATH_OBS_FILE_CG5 + settings.NAME_OBS_FILE_CG5

    s1 = CG5Survey()
    s1.read_obs_file(path)

    # s1.plot_g_values(['1-164-04', '1-164-12', '1-164-11'])
    # s1.plot_g_values(['TEST'])
    s1.plot_g_values()

else:
    # not run as standalone program, but as module
    pass
