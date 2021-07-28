"""Model classes for pyQt5's model view architecture."""

from PyQt5.QtCore import QAbstractTableModel, Qt
from PyQt5 import QtGui
from PyQt5.QtWidgets import QMessageBox
import datetime as dt
import pandas as pd

from gravtools.models.survey import Survey
from gravtools.models.lsm_diff import LSMDiff
from gravtools.models.lsm_nondiff import LSMNonDiff
from gravtools.models.mlr_bev_legacy import BEVLegacyProcessing

NONE_REPRESENTATION_IN_TABLE_VIEW = ''  # Representation of None values in table views in the GUI


class StationTableModel(QAbstractTableModel):
    """Model for displaying the station data in a table view (QTableView)."""

    def __init__(self, stat_df):
        QAbstractTableModel.__init__(self)
        self.load_stat_df(stat_df)

    def load_stat_df(self, stat_df):
        """Load station data from pandas dataframe to table model."""
        # self._data = stat_df.copy()
        self._data = stat_df

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:

                value = self._data.iloc[index.row(), index.column()]
                # Custom formatter (string is expected as return type):
                if value is None:  #
                    return NONE_REPRESENTATION_IN_TABLE_VIEW
                elif isinstance(value, float):
                    if value != value:  # True, if value is "NaN"
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    else:
                        number_of_decimal_places_per_column = {
                            1: 3,
                            2: 3,
                            3: 3,
                            4: 1,
                            5: 1,
                            6: 1,
                        }
                        if index.column() in number_of_decimal_places_per_column.keys():
                            num_dec_places = number_of_decimal_places_per_column[index.column()]
                            return '{1:.{0}f}'.format(num_dec_places, value)
                else:  # all other
                    return str(value)

            if role == Qt.TextAlignmentRole:
                value = self._data.iloc[index.row(), index.column()]
                if isinstance(value, int) or isinstance(value, float):
                    # Align right, vertical middle.
                    return Qt.AlignVCenter + Qt.AlignRight

            if role == Qt.BackgroundRole:
                is_observed_flag = self._data.iloc[index.row(), 7]  # is_observed
                if is_observed_flag:
                    return QtGui.QColor('cyan')
        return None

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data.columns[section])
            if orientation == Qt.Vertical:
                return str(self._data.index[section])

    def setData(self, index, value, role):
        """Example: https://www.semicolonworld.com/question/58510/how-to-display-a-pandas-data-frame-with-pyqt5"""
        if not index.isValid():
            return False
        if role == Qt.EditRole:
            # Get column and row indices for dataframe:
            row = self._data.index[index.row()]
            col = self._data.columns[index.column()]

            # Handle data type in a generic way, if required:
            # - list of available dtypes: https://stackoverflow.com/questions/37561991/what-is-dtypeo-in-pandas
            # dtype_col = self._data[col].dtype  # dtype of the dataframe columns
            # if dtype != object:  # "string"
            #     value = None if value == '' else dtype.type(value)  # changes dtype of the dataframe...

            if col == 'is_datum':
                # convert "value" (str) to bool and set itm in dataframe:
                if value == 'True':
                    self._data.at[row, col] = True
                    self.dataChanged.emit(index, index)  # Is it necessary?
                elif value == 'False':
                    self._data.at[row, col] = False
                    self.dataChanged.emit(index, index)  # Is it necessary?
                else:
                    # print(f'Input not valid: {value}')
                    QMessageBox.warning(self.parent(), 'Warning!',
                                        f'Input "{value}" not valid! Only "True" or "False" allowed.')
                    return False
                return True  # Data successfully set

        return False

    def flags(self, index):
        """Enable editing of table items."""
        flags = super(self.__class__, self).flags(index)
        flags |= Qt.ItemIsSelectable
        flags |= Qt.ItemIsEnabled
        flags |= Qt.ItemIsDragEnabled
        flags |= Qt.ItemIsDropEnabled
        # flags |= Qt.ItemIsUserCheckable
        if index.column() == 10:  # Column: "is_datum"
            flags |= Qt.ItemIsEditable
        return flags

    @property
    def get_data(self):
        return self._data


class SetupTableModel(QAbstractTableModel):
    """Model for displaying the station data.

    Attributes
    ----------
    _surveys : dict of py.obj:`gravtools.survey.Survey`
        Dict of survey object that are assigne by reference to the campaing data object
        (py.obj:`gravtools.survey.Campaign`).
    _data : pandas dataframe
        Data frame that holds the observation data of a specific setup or survey to be displayed in the table view.
    _data_survey_name : str
        Name of the currently selected survey.
    """

    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {  # key = column name; value = number of decimal places
        'g_mugal': 1,
        'sd_g_mugal': 1,
        'epoch_unix': 1,
        'delta_t_h': 3,
    }

    def __init__(self, surveys):
        """Initialize the observation table view model.

        Parameters
        ----------
        surveys : dict of py.obj:`gravtools.survey.Survey` objects
        """
        QAbstractTableModel.__init__(self)
        self._surveys = {}
        self._data = None  # Observations (or at subset of them) of the survey with the name `self._data_survey_name`
        self.load_surveys(surveys)
        self._data_survey_name = ''  # Name of the Survey that is currently represented by `self._data`

    def load_surveys(self, surveys):
        """Load observation data (dict of survey objects in the campaign object) to the observation model.

        Notes
        -----
        The data is assigned by reference, i.e. all changes in `_surveys` will propagate to the data origin.
        """
        self._surveys = surveys

    def update_view_model(self, survey_name, setup_id):
        """Update the `_data` DataFrame that hold the actual data that is displayed."""

        try:
            setup_df = self._surveys[survey_name].setup_df
        except KeyError:
            QMessageBox.critical(self.parent(), 'Error!', f'Survey "{survey_name}" is not available in this campaign.')
        else:
            self._data_survey_name = survey_name
            if setup_id is None:  # No setup ID provided => Take all observations in survey
                if setup_df is None:
                    self._data = None
                else:
                    self._data = setup_df.copy(deep=True)
            else:  # Only take observations of the specified setup
                if setup_df is None:
                    self._data = None
                else:
                    self._data = setup_df[setup_df['setup_id'] == setup_id].copy(deep=True)

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if self._data is not None:  # No data available yet
            if role == Qt.DisplayRole:
                if orientation == Qt.Horizontal:
                    return str(self._data.columns[section])
                if orientation == Qt.Vertical:
                    return str(self._data.index[section])

    def rowCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[0]
        else:
            return 0

    def columnCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[1]
        else:
            return 0

    def data(self, index, role=Qt.DisplayRole):
        if self._data is not None:
            if index.isValid():
                if role == Qt.DisplayRole:
                    value = self._data.iloc[index.row(), index.column()]
                    # Custom formatter (string is expected as return type):
                    if value is None:  #
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    elif isinstance(value, float):
                        if value != value:  # True, if value is "NaN"
                            return NONE_REPRESENTATION_IN_TABLE_VIEW
                        else:
                            try:
                                col_name_str = Survey.get_setup_df_column_name(index.column())
                                num_dec_places = self._DECIMAL_PLACES_PER_FLOAT_COLUMN[col_name_str]
                                return '{1:.{0}f}'.format(num_dec_places, value)
                            except KeyError:
                                return str(value)
                    elif isinstance(value, dt.datetime):
                        return value.strftime("%Y-%m-%d, %H:%M:%S")
                    else:  # all other
                        return str(value)

                if role == Qt.TextAlignmentRole:
                    value = self._data.iloc[index.row(), index.column()]
                    if isinstance(value, int) or isinstance(value, float):
                        # Align right, vertical middle.
                        return Qt.AlignVCenter + Qt.AlignRight


class ObservationTableModel(QAbstractTableModel):
    """Model for displaying the observation data.

    Attributes
    ----------
    _surveys : dict of py.obj:`gravtools.survey.Survey`
        Dict of survey object that are assigne by reference to the campaing data object
        (py.obj:`gravtools.survey.Campaign`).
    _data : pandas dataframe
        Data frame that holds the observation data of a specific setup or survey to be displayed in the table view.
    _data_survey_name : str
        Name of the currently selected survey.
    _setup_data : pandas DataFrame
        Contains one observation per setup.

    """

    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {
        'lon_deg': 3,
        'lat_deg': 3,
        'alt_m': 3,
        'g_obs_mugal': 1,
        'sd_g_obs_mugal': 1,
        'g_red_mugal': 1,
        'sd_g_red_mugal': 1,
        'corr_terrain': 1,
        'corr_tide_mugal': 1,
        'corr_tide_red_mugal': 1,
        'temp': 1,
        'tiltx': 1,
        'tilty': 1,
        'dhf_m': 3,
        'dhb_m': 3,
        'vg_mugalm': 1,
    }

    def __init__(self, surveys):
        """Initialize the observation table view model.

        Parameters
        ----------
        surveys : dict of py.obj:`gravtools.survey.Survey` objects
        """
        self._surveys = {}
        self._data = None  # Observations (or at subset of them) of the survey with the name `self._data_survey_name`
        QAbstractTableModel.__init__(self)
        self.load_surveys(surveys)
        self._data_survey_name = ''  # Name of the Survey that is currently represented by `self._data`
        self._setup_data = None  # Setup data (or at subset) of the survey with the name `self._data_survey_name`

    def load_surveys(self, surveys):
        """Load observation data (dict of survey objects in the campaign object) to the observation model.

        Notes
        -----
        The data is assigned by reference, i.e. all changes in `_surveys` will propagate to the data origin.
        """
        self._surveys = surveys

    def update_view_model(self, survey_name, setup_id):
        """Update the `_data` data frame that hold the actual data that is viewed."""

        try:
            obs_df = self._surveys[survey_name].obs_df  # Select the survey
            setup_df = self._surveys[survey_name].setup_df
        except KeyError:
            QMessageBox.critical(self.parent(), 'Error!', f'Survey "{survey_name}" is not available in this campaign.')
        else:
            self._data_survey_name = survey_name
            if setup_id is None:  # No setup ID provided => Take all observations in survey
                self._data = obs_df.copy(deep=True)
                if setup_df is None:
                    self._setup_data = None
                else:
                    self._setup_data = setup_df.copy(deep=True)
            else:  # Only take observations of the specified setup
                self._data = obs_df[obs_df['setup_id'] == setup_id].copy(deep=True)
                if setup_df is None:
                    self._setup_data = None
                else:
                    self._setup_data = setup_df[setup_df['setup_id'] == setup_id].copy(deep=True)

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if self._data is not None:  # No data available yet
            if role == Qt.DisplayRole:
                if orientation == Qt.Horizontal:
                    return str(self._data.columns[section])
                if orientation == Qt.Vertical:
                    return str(self._data.index[section])

    def rowCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[0]
        else:
            return 0

    def columnCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[1]
        else:
            return 0

    def data(self, index, role=Qt.DisplayRole):
        if self._data is not None:
            if index.isValid():
                if role == Qt.DisplayRole:
                    value = self._data.iloc[index.row(), index.column()]
                    # Custom formatter (string is expected as return type):
                    if value is None:  #
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    elif isinstance(value, float):
                        if value != value:  # True, if value is "NaN"
                            return NONE_REPRESENTATION_IN_TABLE_VIEW
                        else:
                            try:
                                col_name_str = Survey.get_obs_df_column_name(index.column())
                                num_dec_places = self._DECIMAL_PLACES_PER_FLOAT_COLUMN[col_name_str]
                                return '{1:.{0}f}'.format(num_dec_places, value)
                            except:
                                return str(value)
                    elif isinstance(value, dt.datetime):
                        return value.strftime("%Y-%m-%d, %H:%M:%S")
                    else:  # all other
                        return str(value)

                if role == Qt.TextAlignmentRole:
                    value = self._data.iloc[index.row(), index.column()]
                    if isinstance(value, int) or isinstance(value, float):
                        # Align right, vertical middle.
                        return Qt.AlignVCenter + Qt.AlignRight

                if role == Qt.BackgroundRole:
                    # keep_obs_index = self.observation_model.get_data.columns.to_list().index('keep_obs')
                    # keep_obs_flag = self._data.iloc[index.row(), keep_obs_index]
                    keep_obs_flag = self._data.iloc[index.row(), 18]  # keep_obs
                    if not keep_obs_flag:
                        return QtGui.QColor('red')

    def flags(self, index):
        """Enable editing of table items."""
        flags = super(self.__class__, self).flags(index)
        flags |= Qt.ItemIsSelectable
        flags |= Qt.ItemIsEnabled
        flags |= Qt.ItemIsDragEnabled
        flags |= Qt.ItemIsDropEnabled
        col_name_str = Survey.get_obs_df_column_name(index.column())
        if col_name_str == 'keep_obs':
            flags |= Qt.ItemIsEditable
        return flags

    def setData(self, index, value, role):
        """Example: https://www.semicolonworld.com/question/58510/how-to-display-a-pandas-data-frame-with-pyqt5"""
        if not index.isValid():
            return False
        if role == Qt.EditRole:
            # Get column and row indices for dataframe:
            row = self._data.index[index.row()]
            col = self._data.columns[index.column()]

            if col == 'keep_obs':
                # convert "value" (str) to bool and set itm in dataframe:
                if value == 'True':
                    self._data.at[row, col] = True
                    self.dataChanged.emit(index, index)  # Is it necessary?
                    # # Change data in `obs_df`:
                    # obs_df_row_index_int = self._data.index[index.row()]
                    # self._surveys[self._data_survey_name].obs_df.iat[
                    #     obs_df_row_index_int, Survey.get_obs_df_column_index('keep_obs')] = True
                elif value == 'False':
                    self._data.at[row, col] = False
                    self.dataChanged.emit(index, index)  # Is it necessary?
                    # # Change data in `obs_df`:
                    # obs_df_row_index_int = self._data.index[index.row()]
                    # self._surveys[self._data_survey_name].obs_df.iat[
                    #     obs_df_row_index_int, Survey.get_obs_df_column_index('keep_obs')] = False
                else:
                    QMessageBox.warning(self.parent(), 'Warning!',
                                        f'Input "{value}" not valid! Only "True" or "False" allowed.')
                    return False
                return True  # Data successfully set
        return False

    @property
    def data_survey_name(self):
        return self._data_survey_name

    @property
    def get_data(self):
        return self._data

    @property
    def get_setup_data(self):
        return self._setup_data


class ResultsObservationModel(QAbstractTableModel):
    """Model for displaying the observations-related results."""

    # Number of decimal pla
    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {
        'g_diff_mugal': 1,
        'sd_g_diff_mugal': 1,
        'sd_g_diff_est_mugal': 1,
        'v_diff_mugal': 1,
        'sd_g_mugal': 1,
        'g_mugal': 1,
        'abw_mugal': 1,
        'delta_t_h': 3,
        'corr_drift_mugal': 1,
        'w_diff_mugal': 3,
        'r_diff_obs': 3,
        'ref_epoch_dt': 3,
        'w_obs_est_mugal': 3,  # LSM_non_diff
        'r_obs_est': 3,  # LSM_non_diff
        'v_obs_est_mugal': 1,  # LSM_non_diff
        'sd_g_obs_est_mugal': 1,  # LSM_non_diff
        'g_obs_mugal': 1,  # LSM_non_diff
        'sd_g_obs_mugal': 1,  # LSM_non_diff
    }

    # Columns that will be shown in the table view, if available in the data (Also defines the order of columns):
    # - This is also e pre-selection for data columns that can be plotted in the observation plots in the results tab
    # - keys: Actual names of the dataframe columns
    # - items: Header names for the Table View Widget
    _SHOW_COLUMNS_IN_TABLE_DICT = {
        'survey_name': 'Survey',  # LSM_non_diff, LSM_diff
        'station_name_from': 'Station from',  # LSM_diff
        'station_name_to': 'Station to',  # LSM_diff
        'station_name': 'Station',  # MLR, LSM_non_diff
        'ref_epoch_dt': 'Epoch',  # LSM_non_diff, LSM_diff
        'delta_t_h': 'd_t [h]',  # LSM_non_diff, LSM_diff
        'g_mugal': 'g [µGal]',  # MLR
        'g_diff_mugal': 'g [µGal]',  # LSM_diff
        'g_obs_mugal': 'g [µGal]',  # LSM_non_diff
        'sd_g_mugal': 'SD [µGal]',  # MLR
        'sd_g_diff_mugal': 'SD [µGal]',  # LSM_diff
        'sd_g_obs_mugal': 'SD [µGal]',  # LSM_non_diff
        'sd_g_diff_est_mugal': 'SD_est [µGal]',  # LSM_diff
        'sd_g_obs_est_mugal': 'SD_est [µGal]',  # LSM_non_diff
        'v_diff_mugal': 'v [µGal]',  # LSM_diff
        'v_obs_est_mugal': 'v [µGal]',  # LSM_non_diff
        'w_diff_mugal': 'w [µGal]',  # LSM_diff
        'w_obs_est_mugal': 'w [µGal]',  # LSM_non_diff
        'r_diff_obs': 'r [µGal]',  # LSM_diff
        'r_obs_est': 'r [µGal]',  # LSM_non_diff
        'tau_test_result': 'Outlier test',  # LSM_non_diff, LSM_diff
        'corr_drift_mugal': 'Drift corr. [µGal]',  # MLR: Estimated drift correction
        'abw_mugal': 'abw [µGal]',  # MLR: Difference between drift-corrected instrument reading and the estimated station gravity
    }
    _SHOW_COLUMNS_IN_TABLE = list(_SHOW_COLUMNS_IN_TABLE_DICT.keys())  # Actual list of columns to be shown

    # List of
    _SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI = [
        'survey_name',
        'station_name_from',
        'station_name_to',
        'station_name',
        'ref_epoch_dt',
        'delta_t_h',
        'g_mugal',
        'g_diff_mugal',
        'g_obs_mugal',
        'sd_g_mugal',
        'sd_g_diff_mugal',
        'sd_g_obs_mugal',
        'sd_g_diff_est_mugal',
        'sd_g_obs_est_mugal',
        'v_diff_mugal',
        'v_obs_est_mugal',
        'tau_test_result',
    ]

    # Columns that can be plotted. Column names (keys) and description (values):
    # - numerical data only!
    _PLOTABLE_DATA_COLUMNS = {
        # LSM_diff:
        'v_diff_mugal': 'Post-fit residuals [µGal]',  # LSM_diff
        'sd_g_diff_est_mugal': 'A posteriori SD of diff. obs. [µGal]',  # LSM_diff
        'g_diff_mugal': 'Differential observation [µGal]',  # LSM_diff
        'sd_g_diff_mugal': 'SD of differential observation [µGal]',  # LSM_diff
        # MLR:
        'g_mugal': 'Instrument reading [µGal]',  # MLR
        'sd_g_mugal': 'Standard deviation of the instrument reading [µGal]',  # MLR
        'corr_drift_mugal': 'Estimated drift correction [µGal]',  # MLR
        'abw_mugal': 'Drift-corrected reading minus estimated station gravity (not absolute) [µGal]',  # MLR
        'r_diff_obs': 'Redundancy components []',  # LSM_diff
        'w_diff_mugal': 'Standardized post-fit residuals []',  # LSM_diff
        # LSM_non_diff:
        'w_obs_est_mugal': 'Standardized post-fit residuals []',  # LSM_non_diff
        'r_obs_est': 'Redundancy components []',  # LSM_non_diff
        'v_obs_est_mugal': 'Post-fit residuals [µGal]',  # LSM_non_diff
        'sd_g_obs_est_mugal': 'A posteriori SD of observation [µGal]',  # LSM_non_diff
        'g_obs_mugal': 'Non-differential observation [µGal]',  # LSM_non_diff
        'sd_g_obs_mugal': 'SD of non-differential observation [µGal]',  # LSM_non_diff
    }

    def __init__(self, lsm_runs):
        """Initialize the observation-results table view model.

        Parameters
        ----------
        lsm_runs : list of py.obj:`gravtools.lsm.LSM` objects
        """
        QAbstractTableModel.__init__(self)
        self._lsm_runs = []
        self._data = None  # Observations (or at subset of them) of the survey with the name `self._data_survey_name`
        self.load_lsm_runs(lsm_runs)
        self._lsm_run_index = None  # Name of the Survey that is currently represented by `self._data`
        self._data_column_names = None
        self.flag_gui_simple_mode = False

    def load_lsm_runs(self, lsm_runs: list):
        """Load adjustment results.

        Notes
        -----
        The data is assigned by reference, i.e. all changes in `_surveys` will propagate to the data origin.
        """
        self._lsm_runs = lsm_runs

    def update_view_model(self, lsm_run_index: int, station_name=None, survey_name=None, gui_simple_mode=False):
        """Update the `_data` DataFrame that hold the actual data that is displayed."""
        self.flag_gui_simple_mode = gui_simple_mode
        if lsm_run_index == -1:  # No data available => Invalid index => Reset model data
            self._data = None
            self._lsm_run_index = None
            self._data_column_names = None
        else:
            try:
                results_obs_df = self._lsm_runs[lsm_run_index].get_results_obs_df
            except KeyError:
                QMessageBox.critical(self.parent(), 'Error!', f'LSM run with index "{lsm_run_index}" not found!')
            else:
                self._lsm_run_index = lsm_run_index

                # Get list of columns to be depicted via the table model:
                # - Keep order of itemms in `self._SHOW_COLUMNS_IN_TABLE`
                results_obs_df_columns_set = frozenset(results_obs_df.columns)
                table_model_columns = [x for x in self.get_table_columns if x in results_obs_df_columns_set]

                # Apply filter:
                if ((station_name is not None) or (survey_name is not None)) and (results_obs_df is not None):
                    tmp_filter = pd.Series([True] * len(results_obs_df))  # Init boolean filter series
                    column_names = results_obs_df.columns

                    if station_name is not None:  # Filer data for station names
                        flag_is_diff_obs = ('station_name_from' in column_names) and ('station_name_to' in column_names)

                        # Differential observations => Check columns 'station_name_from' and 'station_name_to'
                        if flag_is_diff_obs:
                            tmp_filter = tmp_filter & ((results_obs_df['station_name_from'] == station_name) |
                                                       (results_obs_df['station_name_to'] == station_name))
                        else:  # No differential observations => Check column 'station_name'
                            tmp_filter = tmp_filter & (results_obs_df['station_name'] == station_name)

                    if survey_name is not None:  # Filer data for survey name
                        tmp_filter = tmp_filter & (results_obs_df['survey_name'] == survey_name)
                    try:
                        self._data = results_obs_df.loc[tmp_filter, table_model_columns].copy(deep=True)
                        # self._data = self._data.reindex(columns=table_model_columns)  # Reoorder columns
                    except:  # Just in case any problem occurs...
                        self._data = results_obs_df.loc[tmp_filter, :].copy(deep=True)  # Show all columns
                else:  # No filter
                    try:
                        self._data = results_obs_df.loc[:, table_model_columns].copy(deep=True)
                        # self._data = self._data.reindex(columns=table_model_columns)  # Reoorder columns
                    except:  # Just in case any problem occurs...
                        self._data = results_obs_df.copy(deep=True)
                self._data_column_names = self._data.columns.to_list()  # Show all columns

    def rowCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[0]
        else:
            return 0

    def columnCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[1]
        else:
            return 0

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                value = self._data.iloc[index.row(), index.column()]
                column_name = self._data_column_names[index.column()]
                # Custom formatter (string is expected as return type):
                if value is None:  #
                    return NONE_REPRESENTATION_IN_TABLE_VIEW
                elif isinstance(value, float):
                    if value != value:  # True, if value is "NaN"
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    else:
                        if column_name in self._DECIMAL_PLACES_PER_FLOAT_COLUMN.keys():
                            num_dec_places = self._DECIMAL_PLACES_PER_FLOAT_COLUMN[column_name]
                            return '{1:.{0}f}'.format(num_dec_places, value)
                        else:
                            return str(value)
                elif isinstance(value, dt.datetime):
                    return value.strftime("%Y-%m-%d, %H:%M:%S")
                else:  # all other
                    return str(value)

            if role == Qt.TextAlignmentRole:
                value = self._data.iloc[index.row(), index.column()]
                if isinstance(value, int) or isinstance(value, float):
                    # Align right, vertical middle.
                    return Qt.AlignVCenter + Qt.AlignRight

            if role == Qt.BackgroundRole:
                if 'tau_test_result' in self._data_column_names:
                    if index.column() == self._data_column_names.index('tau_test_result'):
                        value = self._data.iloc[index.row(), index.column()]
                        if value == 'passed':
                            return QtGui.QColor('green')
                        elif value == 'failed':
                            return QtGui.QColor('red')
        return None

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._SHOW_COLUMNS_IN_TABLE_DICT[str(self._data.columns[section])]
            if orientation == Qt.Vertical:
                return str(self._data.index[section])

    def get_plotable_columns(self) -> dict:
        """Returns a dict with names and description of the plotable numerical columns.

        Notes
        -----
        If the model data is empty return an empty dict.
        """
        plotable_columns_dict = {}
        if self._data_column_names is not None:  # Data model is not empty
            for key, item in self._PLOTABLE_DATA_COLUMNS.items():
                # Check if key is a column-name in the current data
                if key in self._data_column_names:
                    plotable_columns_dict[key] = item
        return plotable_columns_dict

    @property
    def get_model_data_df(self):
        """Returns the model data dataframe."""
        return self._data

    @property
    def get_model_data_df_for_plotting(self):
        """Return all plotable rows on the model data dataframe.

        Notes
        -----
        Data is plotable, if the data can be referred station(s) and to a specific reference epoch! E.g. pseudo
        observations introduced for datum constraints are not plotable!
        """
        filter_tmp = self._data['is_constraint'] == False
        # Add further options to thze filter series, is needed!
        return self._data.loc[filter_tmp, :]

    @property
    def get_table_columns(self):
        """Returns a list with all columns names of the dataframe that should be copied to the view model."""
        if self.flag_gui_simple_mode:
            return [value for value in self._SHOW_COLUMNS_IN_TABLE if value in self._SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI]
        else:
            return self._SHOW_COLUMNS_IN_TABLE

    def get_short_column_description(self, column_name: str) -> str:
        """Returns the short description of the model column."""
        try:
            return self._SHOW_COLUMNS_IN_TABLE_DICT[column_name]
        except AttributeError:
            return ''


class ResultsStationModel(QAbstractTableModel):
    """Model for displaying the stations-related results."""

    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {
        'lon_deg': 3,
        'lat_deg': 3,
        'height_m': 3,
        'g_mugal': 1,
        'sd_g_mugal': 1,
        'vg_mugalm': 1,
        'g_est_mugal': 1,
        'sd_g_est_mugal': 1,
        'diff_g_est_mugal': 1,
        'diff_sd_g_est_mugal': 1,
        'g0_mugal': 1,
        'sig_g0_mugal': 1,
    }

    _SHOW_COLUMNS_IN_TABLE = [
        'station_name',
        'g_mugal',
        'sd_g_mugal',
        'is_datum',
        'g_est_mugal',
        'sd_g_est_mugal',
        'diff_g_est_mugal',
        'diff_sd_g_est_mugal',
        'g0_mugal',
        'sig_g0_mugal',
    ]

    def __init__(self, lsm_runs):
        """Initialize the station-results table view model.

        Parameters
        ----------
        lsm_runs : list of py.obj:`gravtools.lsm.LSM` objects
        """
        QAbstractTableModel.__init__(self)
        self._lsm_runs = []
        self._data = None  # Observations (or at subset of them) of the survey with the name `self._data_survey_name`
        self.load_lsm_runs(lsm_runs)
        self._lsm_run_index = None  # Name of the Survey that is currently represented by `self._data`
        self._data_column_names = None

    def load_lsm_runs(self, lsm_runs: list):
        """Load adjustment results.

        Notes
        -----
        The data is assigned by reference, i.e. all changes in `_surveys` will propagate to the data origin.
        """
        self._lsm_runs = lsm_runs

    def update_view_model(self, lsm_run_index: int, station_name=None, survey_name=None):
        """Update the `_data` DataFrame that hold the actual data that is displayed.

        Notes
        -----
        Data selection based on survey name (only display stations that were observed in the selected survey) is not
        implemented yet.
        """
        if lsm_run_index == -1:  # No data available => Invalid index => Reset model data
            self._data = None
            self._lsm_run_index = None
            self._data_column_names = None
        else:
            try:
                results_stat_df = self._lsm_runs[lsm_run_index].get_results_stat_df
            except KeyError:
                QMessageBox.critical(self.parent(), 'Error!', f'LSM run with index "{lsm_run_index}" not found!')
            else:
                self._lsm_run_index = lsm_run_index

                # Get list of columns to be depicted via the table model:
                # - Keep order of itemms in `self._SHOW_COLUMNS_IN_TABLE`
                results_stat_df_columns_set = frozenset(results_stat_df.columns)
                table_model_columns = [x for x in self._SHOW_COLUMNS_IN_TABLE if x in results_stat_df_columns_set]

                if (station_name is not None) and (results_stat_df is not None):
                    tmp_filter = results_stat_df['station_name'] == station_name
                    self._data = results_stat_df.loc[tmp_filter, table_model_columns].copy(deep=True)
                else:  # No filter
                    self._data = results_stat_df.loc[:, table_model_columns].copy(deep=True)
                self._data_column_names = self._data.columns.to_list()

    def rowCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[0]
        else:
            return 0

    def columnCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[1]
        else:
            return 0

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                value = self._data.iloc[index.row(), index.column()]
                column_name = self._data_column_names[index.column()]
                # Custom formatter (string is expected as return type):
                if value is None:  #
                    return NONE_REPRESENTATION_IN_TABLE_VIEW
                elif isinstance(value, float):
                    if value != value:  # True, if value is "NaN"
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    else:
                        if column_name in self._DECIMAL_PLACES_PER_FLOAT_COLUMN.keys():
                            num_dec_places = self._DECIMAL_PLACES_PER_FLOAT_COLUMN[column_name]
                            return '{1:.{0}f}'.format(num_dec_places, value)
                        else:
                            return str(value)
                else:  # all other
                    return str(value)

            if role == Qt.TextAlignmentRole:
                value = self._data.iloc[index.row(), index.column()]
                if isinstance(value, int) or isinstance(value, float):
                    # Align right, vertical middle.
                    return Qt.AlignVCenter + Qt.AlignRight
        return None

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data.columns[section])
            if orientation == Qt.Vertical:
                return str(self._data.index[section])


class ResultsDriftModel(QAbstractTableModel):
    """Model for displaying the drift-related results."""

    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {
        'coefficient': 3,
        'sd_coeff': 3,
    }

    # Column names (keys) and short description (values):
    _PLOT_COLUMNS_DICT = {
        'survey_name': 'Survey',
        'degree': 'Degree',
        'coefficient': 'Coefficient',
        'sd_coeff': 'SD',
        'coeff_unit': 'Unit',
    }
    _PLOT_COLUMNS = list(_PLOT_COLUMNS_DICT.keys())  # Actual list of columns to be shown

    def __init__(self, lsm_runs):
        """Initialize the drift-results table view model.

        Parameters
        ----------
        lsm_runs : list of py.obj:`gravtools.lsm.LSM` objects
        """
        QAbstractTableModel.__init__(self)
        self._lsm_runs = []
        self._data = None  # Observations (or at subset of them) of the survey with the name `self._data_survey_name`
        self.load_lsm_runs(lsm_runs)
        self._lsm_run_index = None  # Name of the Survey that is currently represented by `self._data`
        self._data_column_names = None

    def load_lsm_runs(self, lsm_runs: list):
        """Load adjustment results.

        Notes
        -----
        The data is assigned by reference, i.e. all changes in `_surveys` will propagate to the data origin.
        """
        self._lsm_runs = lsm_runs

    def update_view_model(self, lsm_run_index: int, survey_name=None):
        """Update the `_data` DataFrame that hold the actual data that is displayed."""
        if lsm_run_index == -1:  # No data available => Invalid index => Reset model data
            self._data = None
            self._lsm_run_index = None
            self._data_column_names = None
        else:
            try:
                results_drift_df = self._lsm_runs[lsm_run_index].get_results_drift_df
            except KeyError:
                QMessageBox.critical(self.parent(), 'Error!', f'LSM run with index "{lsm_run_index}" not found!')
            else:
                self._lsm_run_index = lsm_run_index

                if survey_name is None:  # No filter
                    self._data = results_drift_df.loc[:, self._PLOT_COLUMNS].copy(deep=True)
                else:  # Filer data for survey name
                    if results_drift_df['survey_name'].isna().all():
                        # On drift polynomial estimated for all surveys in campaign
                        self._data = results_drift_df.loc[:, self._PLOT_COLUMNS].copy(deep=True)
                    else:
                        tmp_filter = results_drift_df['survey_name'] == survey_name
                        self._data = results_drift_df.loc[tmp_filter, self._PLOT_COLUMNS].copy(deep=True)
                self._data_column_names = self._data.columns.to_list()

    def rowCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[0]
        else:
            return 0

    def columnCount(self, parent=None):
        if self._data is not None:
            return self._data.shape[1]
        else:
            return 0

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                value = self._data.iloc[index.row(), index.column()]
                column_name = self._data_column_names[index.column()]
                # Custom formatter (string is expected as return type):
                if value is None:  #
                    return NONE_REPRESENTATION_IN_TABLE_VIEW
                elif isinstance(value, float):
                    if value != value:  # True, if value is "NaN"
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    else:
                        if column_name in self._DECIMAL_PLACES_PER_FLOAT_COLUMN.keys():
                            num_dec_places = self._DECIMAL_PLACES_PER_FLOAT_COLUMN[column_name]
                            return '{1:.{0}f}'.format(num_dec_places, value)
                        else:
                            return str(value)
                else:  # all other
                    return str(value)

            if role == Qt.TextAlignmentRole:
                value = self._data.iloc[index.row(), index.column()]
                if isinstance(value, int) or isinstance(value, float):
                    # Align right, vertical middle.
                    return Qt.AlignVCenter + Qt.AlignRight
        return None

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._PLOT_COLUMNS_DICT[str(self._data.columns[section])]
            if orientation == Qt.Vertical:
                return str(self._data.index[section])

    def get_short_column_description(self, column_name: str) -> str:
        """Returns the short description of the model column."""
        try:
            return self._PLOT_COLUMNS_DICT[column_name]
        except AttributeError:
            return ''
