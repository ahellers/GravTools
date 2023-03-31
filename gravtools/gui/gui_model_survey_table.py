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

from PyQt5.QtCore import QAbstractTableModel, Qt
from PyQt5 import QtGui
from PyQt5.QtWidgets import QMessageBox
import datetime as dt

NONE_REPRESENTATION_IN_TABLE_VIEW = ''  # Representation of None values in table views in the GUI


class SurveyTableModel(QAbstractTableModel):
    """Model for displaying the surveys in a table view (QTableView)."""

    # Columns that will be shown in the table view:
    # - keys: Actual attribute names of the Survey object
    # - items: Header names for the Table View Widget
    _SHOW_ATTRIBUTES_IN_TABLE_DICT = {
        'name': 'Name',
        'keep_survey': 'Use data',
        'date': 'Date',
        'start_time_hhmmss_str': 'Start time',
        'end_time_hhmmss_str': 'End time',
        'duration_hhmmss_str': 'Duration',
        'operator': 'Operator',
        'institution': 'Institution',
        'gravimeter_type': 'Gravimeter',
        'gravimeter_serial_number': 'S/N',
        'number_of_setups': 'n/o setups',
        'number_of_observations': 'n/o obs.',
    }
    # _ATTRIBUTE_NAMES = list(_SHOW_ATTRIBUTES_IN_TABLE_DICT.keys())  # Actual list of columns to be shown

    # List of columns that are shown in the simple mode of the GUI
    _SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI = [
        # 'name',
        'keep_survey',
        'date',
        'operator',
        'institution',
        'gravimeter_type',
        'gravimeter_serial_number',
        'number_of_setups',
        'number_of_observations',
        'start_time_hhmmss_str',
        'end_time_hhmmss_str',
        'duration_hhmmss_str',
    ]

    _ATTRIBUTES_TOOLTIPS = {
        'name': 'Survey name',
        'keep_survey': '`True` implies that the observations of this survey are used to calculate setup observations.',
        'date': 'Date',
        'operator': 'Operator',
        'institution': 'Institution',
        'gravimeter_type': 'Gravimeter type',
        'gravimeter_serial_number': 'Serial number of the gravimeter',
        'number_of_setups': 'Number of setups in this survey',
        'number_of_observations': 'Number of observations in this setup',
        'start_time_hhmmss_str': 'Start time',
        'end_time_hhmmss_str': 'End time',
        'duration_hhmmss_str': 'Duration',
    }

    # Number of decimal pla
    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {
        'number_of_setups': 0,
        'number_of_observations': 0,
    }

    def __init__(self, surveys: dict, gui_simple_mode=False):
        QAbstractTableModel.__init__(self)
        # self._data = None
        # self._attribute_names = []
        self._show_attributes_in_table_dict = {}
        self._survey_name_list = []
        self.flag_gui_simple_mode = gui_simple_mode
        # self.load_survey(surveys)
        self._data = surveys
        self._survey_name_list = list(surveys.keys())
        if not gui_simple_mode:
            self._show_attributes_in_table_dict = self._SHOW_ATTRIBUTES_IN_TABLE_DICT
        else:
            for key, value in self._SHOW_ATTRIBUTES_IN_TABLE_DICT.items():
                if key in self._SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI:
                    self._show_attributes_in_table_dict[key] = value
        self._attribute_names = list(self._show_attributes_in_table_dict.keys())  # Actual list of columns to be shown

    # def load_surveys(self, surveys):
    #     """Load Survey object."""
    #     self._data = surveys
    #     self._survey_name_list = list(surveys.keys())
    #
    #     # self._data = stat_df.copy()
    #     # Get list of columns to be depicted via the table model:
    #     # - Keep order of items in `self._SHOW_COLUMNS_IN_TABLE`
    #     stat_df_columns_set = frozenset(stat_df.columns)
    #     table_model_columns = [x for x in self.get_table_columns if x in stat_df_columns_set]
    #     # if self.flag_gui_simple_mode:
    #     #     self._data = stat_df.loc[:, table_model_columns]
    #     # else:
    #     #     self._data = stat_df.loc[:, table_model_columns]
    #     self._data = stat_df.loc[:, table_model_columns]
    #     self._data_column_names = self._data.columns.to_list()

    def rowCount(self, parent=None):
        return len(self._data)  # Number of surveys

    def columnCount(self, parent=None):
        return len(self._attribute_names)

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            survey_name = self._survey_name_list[index.row()]
            attribute_name = self._attribute_names[index.column()]
            try:
                value = getattr(self._data[survey_name], attribute_name)
            except AttributeError:
                value = 'n.a.'  # Not available
            # attribute_name = self._attribute_names[index.column()]
            # value = self._data[survey_name]   .iloc[index.row(), index.column()]

            if role == Qt.DisplayRole:

                # Draw checkboxes only in the "is_datum" column:
                if attribute_name == 'keep_survey':
                    return ''

                # Custom formatter (string is expected as return type):
                if value is None:  #
                    return NONE_REPRESENTATION_IN_TABLE_VIEW
                elif isinstance(value, float):
                    if value != value:  # True, if value is "NaN"
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    else:
                        if attribute_name in self._DECIMAL_PLACES_PER_FLOAT_COLUMN.keys():
                            num_dec_places = self._DECIMAL_PLACES_PER_FLOAT_COLUMN[attribute_name]
                            return '{1:.{0}f}'.format(num_dec_places, value)
                        else:
                            return str(value)
                elif isinstance(value, dt.date):
                    return value.strftime(value.strftime('%Y-%m-%d'))

                else:  # all other
                    return str(value)

            if role == Qt.TextAlignmentRole:
                # value = self._data.iloc[index.row(), index.column()]
                if isinstance(value, int) or isinstance(value, float):
                    # Align right, vertical middle.
                    return Qt.AlignVCenter + Qt.AlignRight

            if role == Qt.BackgroundRole:
                if not self._data[survey_name].keep_survey:
                    return QtGui.QColor('red')

            if role == Qt.CheckStateRole:
                try:
                    if attribute_name == 'keep_survey':
                        keep_obs_flag = self._data[survey_name].keep_survey
                        if keep_obs_flag:
                            return Qt.Checked
                        else:
                            return Qt.Unchecked
                except Exception:
                    return None
        return None

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if self._data is not None:
                if orientation == Qt.Horizontal:
                    attribute_name = self._attribute_names[section]
                    return self._show_attributes_in_table_dict[attribute_name]
                if orientation == Qt.Vertical:
                    return self._survey_name_list[section]

        if role == Qt.ToolTipRole:
            if orientation == Qt.Horizontal:
                attribute_name = self._attribute_names[section]
                try:
                    return self._ATTRIBUTES_TOOLTIPS[attribute_name]
                except KeyError:
                    return ''

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
                    QMessageBox.warning(self.parent(), 'Warning!',
                                        f'Input "{value}" not valid! Only "True" or "False" allowed.')
                    return False
                return True  # Data successfully set

        if role == Qt.CheckStateRole:
            row = self._data.index[index.row()]
            col = self._data.columns[index.column()]
            if col == 'is_datum':
                if value == Qt.Unchecked:
                    self._data.at[row, col] = False
                    self.dataChanged.emit(index, index)  # Update only one item
                elif value == Qt.Checked:
                    self._data.at[row, col] = True
                    self.dataChanged.emit(index, index)  # Update only one item
                else:
                    QMessageBox.warning(self.parent(), 'Warning!',
                                        f'Invalid value fpr keep observation flag: "{value}"')
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
        if index.column() == self._attribute_names.index('keep_survey'):
            # flags |= Qt.ItemIsEditable  # Use checkbox only!
            flags |= Qt.ItemIsUserCheckable
        return flags

    @property
    def get_data(self):
        return self._data

    # @property
    # def get_table_columns(self):
    #     """Returns a list with all columns names of the dataframe that should be copied to the view model."""
    #     if self.flag_gui_simple_mode:
    #         return [value for value in self._SHOW_COLUMNS_IN_TABLE if value in self._SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI]
    #     else:
    #         return self._SHOW_COLUMNS_IN_TABLE
