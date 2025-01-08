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

from gravtools import settings

NONE_REPRESENTATION_IN_TABLE_VIEW = ''  # Representation of None values in table views in the GUI


class StationTableModel(QAbstractTableModel):
    """Model for displaying the station data in a table view (QTableView)."""

    # Columns that will be shown in the table view, if available in the data (Also defines the order of columns):
    # - keys: Actual names of the dataframe columns
    # - items: Header names for the Table View Widget
    _SHOW_COLUMNS_IN_TABLE_DICT = {
        'station_name': 'Station',
        'long_deg': 'Lon [°]',
        'lat_deg': 'Lat [°]',
        'height_m': 'h [m]',
        'g_mugal': 'g [µGal]',
        'sd_g_mugal': 'SD [µGal]',
        'vg_mugalm': 'VG [µGal/m]',
        'is_observed': 'Observed',
        'source_type': 'Source type',
        'source_name': ' Source name',
        'is_datum': 'Datum',
        'in_survey': 'Surveys',
    }
    _SHOW_COLUMNS_IN_TABLE = list(_SHOW_COLUMNS_IN_TABLE_DICT.keys())  # Actual list of columns to be shown

    # List of columns that are shown in the simple mode of the GUI
    _SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI = [
        'station_name',
        'long_deg',
        'lat_deg',
        'height_m',
        'g_mugal',
        'sd_g_mugal',
        'vg_mugalm',
        'is_datum',
        'is_observed',
    ]

    # Number of decimal pla
    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {
        'long_deg': 3,
        'lat_deg': 3,
        'height_m': 3,
        'g_mugal': 1,
        'sd_g_mugal': 1,
        'vg_mugalm': 1,
    }

    def __init__(self, stat_df, gui_simple_mode=False):
        QAbstractTableModel.__init__(self)
        self._data_column_names = None
        self._data = None
        self.flag_gui_simple_mode = gui_simple_mode
        self.load_stat_df(stat_df)

    def load_stat_df(self, stat_df):
        """Load station data from pandas dataframe to table model."""
        # self._data = stat_df.copy()
        # Get list of columns to be depicted via the table model:
        # - Keep order of items in `self._SHOW_COLUMNS_IN_TABLE`
        stat_df_columns_set = frozenset(stat_df.columns)
        table_model_columns = [x for x in self.get_table_columns if x in stat_df_columns_set]
        # if self.flag_gui_simple_mode:
        #     self._data = stat_df.loc[:, table_model_columns]
        # else:
        #     self._data = stat_df.loc[:, table_model_columns]
        self._data = stat_df.loc[:, table_model_columns]
        self._data_column_names = self._data.columns.to_list()

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():

            value = self._data.iloc[index.row(), index.column()]
            column_name = self._data_column_names[index.column()]

            if role == Qt.DisplayRole or role == Qt.UserRole:

                # Draw checkboxes only in the "is_datum" column:
                if role == Qt.DisplayRole:
                    if column_name == 'is_datum':
                        return ''

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
                # value = self._data.iloc[index.row(), index.column()]
                if isinstance(value, int) or isinstance(value, float):
                    # Align right, vertical middle.
                    return Qt.AlignVCenter + Qt.AlignRight

            if role == Qt.BackgroundRole:
                if column_name == 'is_datum':
                    if value:
                        return QtGui.QColor(settings.DATUM_STATION_COLOR[0], settings.DATUM_STATION_COLOR[1],
                                            settings.DATUM_STATION_COLOR[2])

                if self._data.iloc[index.row(), self._data_column_names.index('is_observed')]:  # is_observed
                    return QtGui.QColor('cyan')

            if role == Qt.CheckStateRole:
                try:
                    if column_name == 'is_datum':
                        keep_obs_flag = self._data.iloc[index.row(), self._data_column_names.index('is_datum')]
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
                    return self._SHOW_COLUMNS_IN_TABLE_DICT[str(self._data.columns[section])]
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
        if index.column() == self._data_column_names.index('is_datum'):  # Column: "is_datum"
            # flags |= Qt.ItemIsEditable  # Use checkbox only!
            flags |= Qt.ItemIsUserCheckable
        return flags

    @property
    def get_data(self):
        return self._data

    @property
    def get_table_columns(self):
        """Returns a list with all columns names of the dataframe that should be copied to the view model."""
        if self.flag_gui_simple_mode:
            return [value for value in self._SHOW_COLUMNS_IN_TABLE if value in self._SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI]
        else:
            return self._SHOW_COLUMNS_IN_TABLE
