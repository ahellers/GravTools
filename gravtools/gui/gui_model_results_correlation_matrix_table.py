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

from gravtools import settings
from gravtools.models.misc import time_it

NONE_REPRESENTATION_IN_TABLE_VIEW = ''  # Representation of None values in table views in the GUI


class ResultsCorrelationMatrixModel(QAbstractTableModel):
    """Model for displaying the correlation matrix table."""

    # Number of decimal places for displaying the correlation coefficients
    _DECIMAL_PLACES_CORR_COEF = 3

    def __init__(self, lsm_runs):
        """Initialize the correlation matrix table view model.

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

    # @time_it
    def load_lsm_runs(self, lsm_runs: list):
        """Load adjustment results.

        Notes
        -----
        The data is assigned by reference, i.e. all changes in `_surveys` will propagate to the data origin.
        """
        self._lsm_runs = lsm_runs

    # @time_it
    def update_view_model(self, lsm_run_index: int, station_name=None, survey_name=None):
        """Update the `_data` DataFrame that hold the actual data that is displayed.

        Notes
        -----
        Data selection based on a station name or on a survey name is not implemented yet.
        """

        if lsm_run_index == -1:  # No data available => Invalid index => Reset model data
            self._data = None
            self._lsm_run_index = None
            self._data_column_names = None
        else:
            try:
                self._data = self._lsm_runs[lsm_run_index].get_correlation_matrix  # Rxx matrix (np.array)
                if hasattr(self._lsm_runs[lsm_run_index], 'x_estimate_names'):
                    self._data_column_names = self._lsm_runs[
                        lsm_run_index].x_estimate_names  # Names of estimates in same order as in Rxx
                else:
                    self._data_column_names = None
            except KeyError:
                self._data = None
                self._data_column_names = None
                QMessageBox.critical(self.parent(), 'Error!', f'LSM run with index "{lsm_run_index}" not found!')
            except Exception as e:
                QMessageBox.critical(self.parent(), 'Error!', str(e))
                self._data = None
                self._data_column_names = None
            else:
                self._lsm_run_index = lsm_run_index

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
                value = self._data[index.row(), index.column()]  # np.array
                # column_name = self._data_column_names[index.column()]
                # Custom formatter (string is expected as return type):
                if value is None:  #
                    return NONE_REPRESENTATION_IN_TABLE_VIEW
                elif isinstance(value, float):
                    if value != value:  # True, if value is "NaN"
                        return NONE_REPRESENTATION_IN_TABLE_VIEW
                    else:
                        return '{1:.{0}f}'.format(self._DECIMAL_PLACES_CORR_COEF, value)
                else:  # all other
                    return str(value)

            if role == Qt.TextAlignmentRole:
                value = self._data[index.row(), index.column()]
                if isinstance(value, int) or isinstance(value, float):
                    # Align right, vertical middle.
                    return Qt.AlignVCenter + Qt.AlignRight

            if role == Qt.BackgroundRole:
                if index.row() == index.column():
                    return QtGui.QColor(settings.CORRELATION_COEF_DIAG_ELEMENTS)
                else:
                    value = self._data[index.row(), index.column()]
                    if isinstance(value, int) or isinstance(value, float):
                        # Get absolute value:
                        value = abs(value)
                        color_idx = int(value * len(settings.CORRELATION_COEF_COLORS))
                        color_idx = max(0, color_idx)  # color_idx < 0 become 0
                        color_idx = min(len(settings.CORRELATION_COEF_COLORS)-1, color_idx)
                        return QtGui.QColor(settings.CORRELATION_COEF_COLORS[color_idx])
        return None

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if self._data_column_names is None:
                return str(section)
            else:
                if orientation == Qt.Horizontal:
                    # return self._SHOW_COLUMNS_IN_TABLE_DICT[str(self._data.columns[section])]
                    return str(self._data_column_names[section])
                if orientation == Qt.Vertical:
                    # return str(self._data.index[section])
                    return str(self._data_column_names[section])

