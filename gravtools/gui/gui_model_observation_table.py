"""Model classes for PyQt's model view architecture.

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

from PyQt5.QtCore import QAbstractTableModel, Qt, QPersistentModelIndex
from PyQt5 import QtGui
from PyQt5.QtWidgets import QMessageBox, QStyledItemDelegate
import datetime as dt
import numpy as np

from gravtools.models.survey import Survey

NONE_REPRESENTATION_IN_TABLE_VIEW = ''  # Representation of None values in table views in the GUI


class ObservationTableModel(QAbstractTableModel):
    """Model for displaying the observation data.

    Attributes
    ----------
    _surveys : dict of py.obj:`gravtools.survey.Survey`
        Dict of survey object that are assigned by reference to the campaign data object
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
        'atm_pres_hpa': 1,
        'corr_atm_pres_red_mugal': 1,
        'norm_atm_pres_hpa': 1,
        'linear_scale': 5,
        'corr_oceanload_instrument_mugal': 1,
        'corr_oceanload_red_mugal': 1,
        'duration_sec': 0,
    }

    # Columns that will be shown in the table view, if available in the data (Also defines the order of columns):
    # - keys: Actual names of the dataframe columns
    # - items: Header names for the Table View Widget
    _SHOW_COLUMNS_IN_TABLE_DICT = {
        'keep_obs': 'Keep obs.',
        'station_name': 'Station',
        'setup_id': 'Setup ID',
        'loop_id': 'Loop ID',
        'lon_deg': 'Lon [°]',
        'lat_deg': 'Lat [°]',
        'alt_m': 'h [m]',
        'obs_epoch': 'Epoch UTC',
        'g_obs_mugal': 'g_obs [µGal]',
        'sd_g_obs_mugal': 'SD_obs [µGal]',
        'g_red_mugal': 'g_red [µGal]',
        'sd_g_red_mugal': 'SD_red [µGal]',
        'corr_terrain': 'Terrain corr.',
        'corr_tide_mugal': 'Instr. tide corr. [µGal]',
        'corr_tide_red_mugal': 'Tide corr. [µGal]',
        'temp': 'Temp. corr. [?]',
        'tiltx': 'Tilt_x [asec]',
        'tilty': 'Tilt_y [asec]',
        'dhb_m': 'dhb [m]',
        'dhf_m': 'dhf [m]',
        'vg_mugal': 'VG [µGal]',
        'duration_sec': 'Duration [sec]',
        'atm_pres_hpa': 'p [hPa]',
        'norm_atm_pres_hpa': 'pn [hPa]',
        'corr_atm_pres_red_mugal': 'p corr [µGal]',
        'linear_scale': 'lin. scale',
        'corr_oceanload_instrument_mugal': 'Instr. oceanload corr. [µGal]',
        'corr_oceanload_red_mugal': 'Oceanload corr. [µGal]',
    }
    _SHOW_COLUMNS_IN_TABLE = list(_SHOW_COLUMNS_IN_TABLE_DICT.keys())  # Actual list of columns to be shown

    # List of columns that are shown in the simple mode of the GUI
    _SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI = [
        'station_name',
        'setup_id',
        'obs_epoch',
        'g_obs_mugal',
        'sd_g_obs_mugal',
        'g_red_mugal',
        'sd_g_red_mugal',
        'dhb_m',
        'dhf_m',
        'tiltx',
        'tilty',
        'corr_tide_red_mugal',
        'corr_tide_mugal',
        'vg_mugal',
        'duration_sec',
        'keep_obs',
        'atm_pres_hpa',
        'corr_atm_pres_red_mugal',
        'corr_oceanload_instrument_mugal',
        'corr_oceanload_red_mugal',
    ]

    def __init__(self, surveys):
        """Initialize the observation table view model.

        Parameters
        ----------
        surveys : dict of py.obj:`gravtools.survey.Survey` objects
        """
        self._surveys = {}
        self._data = None  # Observations (or at subset of them) of the survey with the name `self._data_survey_name`
        self._data_column_names = None
        QAbstractTableModel.__init__(self)
        self.load_surveys(surveys)
        self._data_survey_name = ''  # Name of the Survey that is currently represented by `self._data`
        self._setup_data = None  # Setup data (or at subset) of the survey with the name `self._data_survey_name`
        self.flag_gui_simple_mode = False
        # self._keep_obs_check_states = dict()  # To keep track of the checkbox states in the keep_obs column

    def reset_model(self):
        """Reset the model to a state without selected data to display."""
        self._data = None  # Observations (or at subset of them) of the survey with the name `self._data_survey_name`
        self._data_column_names = None
        self._data_survey_name = ''  # Name of the Survey that is currently represented by `self._data`
        self._setup_data = None  # Setup data (or at subset) of the survey with the name `self._data_survey_name`

    def load_surveys(self, surveys):
        """Load observation data (dict of survey objects in the campaign object) to the observation model.

        Notes
        -----
        The data is assigned by reference, i.e. all changes in `_surveys` will propagate to the data origin.
        """
        if surveys is None:
            self._surveys = {}
            return
        self._surveys = surveys

    def update_view_model(self, survey_name, setup_id, gui_simple_mode=False):
        """Update the `_data` data frame that hold the actual data that is viewed."""
        self.flag_gui_simple_mode = gui_simple_mode
        if survey_name is not None:
            try:
                obs_df = self._surveys[survey_name].obs_df  # Select the survey
                setup_df = self._surveys[survey_name].setup_df
            except KeyError:
                QMessageBox.critical(self.parent(), 'Error!', f'Survey "{survey_name}" is not available in this campaign.')
            else:
                self._data_survey_name = survey_name

                # Get list of columns to be depicted via the table model:
                # - Keep order of items in `self._SHOW_COLUMNS_IN_TABLE`
                obs_df_columns_set = frozenset(obs_df.columns)
                table_model_columns = [x for x in self.get_table_columns if x in obs_df_columns_set]

                if setup_id is None:  # No setup ID provided => Take all observations in survey
                    # self._data = obs_df.copy(deep=True)
                    self._data = obs_df.loc[:, table_model_columns].copy(deep=True)
                    if setup_df is None:
                        self._setup_data = None
                    else:
                        self._setup_data = setup_df.copy(deep=True)
                else:  # Only take observations of the specified setup
                    # self._data = obs_df[obs_df['setup_id'] == setup_id].copy(deep=True)
                    self._data = obs_df.loc[obs_df['setup_id'] == setup_id, table_model_columns].copy(deep=True)
                    if setup_df is None:
                        self._setup_data = None
                    else:
                        self._setup_data = setup_df[setup_df['setup_id'] == setup_id].copy(deep=True)
                # Column names of the actual table model:
                self._data_column_names = self._data.columns.to_list()
        else:
            self.reset_model()

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if self._data is not None:  # No data available yet
            if role == Qt.DisplayRole:
                if orientation == Qt.Horizontal:
                    return self._SHOW_COLUMNS_IN_TABLE_DICT[str(self._data.columns[section])]
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
                if role == Qt.DisplayRole or role == Qt.UserRole:
                    value = self._data.iloc[index.row(), index.column()]

                    # Only checkboxes in the "keep_obs" flag column:
                    if role == Qt.DisplayRole:
                        if index.column() == self._data_column_names.index('keep_obs'):
                            return ''

                    # Custom formatter (string is expected as return type):
                    if value is None:  #
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    elif isinstance(value, float):
                        if value != value:  # True, if value is "NaN"
                            return NONE_REPRESENTATION_IN_TABLE_VIEW
                        else:
                            try:
                                col_name_str = self._data_column_names[index.column()]
                                num_dec_places = self._DECIMAL_PLACES_PER_FLOAT_COLUMN[col_name_str]
                                return '{1:.{0}f}'.format(num_dec_places, value)
                            except Exception:
                                return str(value)
                    elif isinstance(value, np.int64):
                        return str(value)
                    elif isinstance(value, dt.datetime):
                        return value.strftime("%Y-%m-%d, %H:%M:%S")
                    else:  # all other
                        return str(value)

                if role == Qt.TextAlignmentRole:
                    value = self._data.iloc[index.row(), index.column()]
                    if isinstance(value, float) or isinstance(value, int) or isinstance(value, np.int64):
                        # Align right, vertical middle.
                        return Qt.AlignVCenter + Qt.AlignRight

                if role == Qt.BackgroundRole:
                    try:
                        if not self._data.iloc[index.row(), self._data_column_names.index('keep_obs')]:
                            return QtGui.QColor('red')
                    except Exception:
                        pass

                if role == Qt.CheckStateRole:
                    try:
                        if index.column() == self._data_column_names.index('keep_obs'):
                            keep_obs_flag = self._data.iloc[index.row(), self._data_column_names.index('keep_obs')]
                            if keep_obs_flag:
                                return Qt.Checked
                            else:
                                return Qt.Unchecked
                    except Exception:
                        print(f'ERROR: row: {index.row()}, col: {index.column()}')

    def flags(self, index):
        """Enable editing of table items."""
        flags = super(self.__class__, self).flags(index)
        flags |= Qt.ItemIsSelectable
        flags |= Qt.ItemIsEnabled
        if index.column() == self._data_column_names.index('keep_obs'):
            # flags |= Qt.ItemIsEditable  # No longer needed => Only checkbox in "keep_obs" column!
            flags |= Qt.ItemIsUserCheckable
        return flags

    def setData(self, index, value, role):
        """Example: https://www.semicolonworld.com/question/58510/how-to-display-a-pandas-data-frame-with-pyqt5"""
        if index.isValid():
            if role == Qt.EditRole:
                # Get column and row indices for dataframe:
                row = self._data.index[index.row()]
                col = self._data.columns[index.column()]
                if col == 'keep_obs':
                    # convert "value" (str) to bool and set item in dataframe:
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

            if role == Qt.CheckStateRole:
                row = self._data.index[index.row()]
                col = self._data.columns[index.column()]
                idx_min = self.index(index.row(), 0)
                idx_max = self.index(index.row(), len(self._data_column_names) - 1)
                if col == 'keep_obs':
                    if value == Qt.Unchecked:
                        # print(f'row {index.row()}: Unchecked!')
                        self._data.at[row, col] = False
                        # Change data in `obs_df`:
                        self._surveys[self._data_survey_name].obs_df.iloc[
                            self._data.index[index.row()], Survey.get_obs_df_column_index('keep_obs')] = False
                        self.dataChanged.emit(idx_min, idx_max)  # Change color of obs. table row
                        self.dataChanged.emit(index, index, [9999])  # Change obs. plot and tree
                    elif value == Qt.Checked:
                        # print(f'row {index.row()}: Checked!')
                        self._data.at[row, col] = True
                        # Change data in `obs_df`:
                        self._surveys[self._data_survey_name].obs_df.iloc[
                            self._data.index[index.row()], Survey.get_obs_df_column_index('keep_obs')] = True
                        self.dataChanged.emit(idx_min, idx_max)  # Change color of obs. table row
                        self.dataChanged.emit(index, index, [9999])  # Change obs. plot and tree
                    else:
                        QMessageBox.warning(self.parent(), 'Warning!',
                                            f'Invalid value fpr keep observation flag: "{value}"')
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

    @property
    def get_table_columns(self):
        """Returns a list with all columns names of the dataframe that should be copied to the view model."""
        if self.flag_gui_simple_mode:
            return [value for value in self._SHOW_COLUMNS_IN_TABLE if value in self._SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI]
        else:
            return self._SHOW_COLUMNS_IN_TABLE

    def get_data_df(self):
        """Returns the model data as pandas DataFrame."""
        print('test')
