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
from PyQt5.QtWidgets import QMessageBox
import datetime as dt

NONE_REPRESENTATION_IN_TABLE_VIEW = ''  # Representation of None values in table views in the GUI


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
        'ref_epoch_t0_dt': 'Ref. epoch',
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
        self._data = None
        self.load_lsm_runs(lsm_runs)
        self._lsm_run_index = None
        self._data_column_names = None

    def load_lsm_runs(self, lsm_runs: list):
        """Load adjustment results.

        Notes
        -----
        The data is assigned by reference.
        """
        self._lsm_runs = lsm_runs

    def update_view_model(self, lsm_run_index: int, survey_name=None):
        """Update the `_data` DataFrame that hold the actual data that is displayed."""
        flag_error_init = False
        if lsm_run_index == -1:  # No data available => Invalid index => Reset model data
            flag_error_init = True
        else:
            try:
                results_drift_df = self._lsm_runs[lsm_run_index].get_results_drift_df
            except KeyError:
                QMessageBox.critical(self.parent(), 'Error!', f'LSM run with index "{lsm_run_index}" not found!')
            else:
                if results_drift_df is None:  # E.g. no results in campaign yet
                    flag_error_init = True
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
        if flag_error_init:
            self._data = None
            self._lsm_run_index = None
            self._data_column_names = None


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
            if role == Qt.DisplayRole or role == Qt.UserRole:
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
        return None

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if self._data is not None:
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
