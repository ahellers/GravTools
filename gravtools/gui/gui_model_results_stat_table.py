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

NONE_REPRESENTATION_IN_TABLE_VIEW = ''  # Representation of None values in table views in the GUI


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
        'se_g_est_mugal': 1,
        'diff_g_est_mugal': 1,
        'diff_sd_g_est_mugal': 1,
        'diff_se_g_est_mugal': 1,
        'g0_mugal': 1,
        'sig_g0_mugal': 1,
        'dhf_sensor_mean_m': 4,
        'dhf_sensor_std_m': 4,
        'dhf_sensor_min_m': 4,
        'dhf_sensor_max_m': 4,
    }

    # Columns that will be shown in the table view, if available in the data (also defines the order of columns):
    # - keys: Actual names of the dataframe columns
    # - items: Header names for the Table View Widget
    _SHOW_COLUMNS_IN_TABLE_DICT = {
        'station_name': 'Station',
        'g_mugal': 'g [µGal]',
        'sd_g_mugal': 'SD [µGal]',
        'is_datum': 'Datum',
        'g_est_mugal': 'g_est [µGal]',
        'sd_g_est_mugal': 'SD_est [µGal]',
        'se_g_est_mugal': 'SE_est [µGal]',
        'diff_g_est_mugal': 'delta_g [µGal]',
        'diff_sd_g_est_mugal': 'delta_SD [µGal]',
        'diff_se_g_est_mugal': 'delta_SE [µGal]',
        'g0_mugal': 'g0 [µGal]',
        'dhf_sensor_mean_m': 'dhf mean [m]',
        'dhf_sensor_std_m': 'dhf SD [m]',
        'dhf_sensor_min_m': 'dhf min [m]',
        'dhf_sensor_max_m': 'dhf max [m]',
    }
    _SHOW_COLUMNS_IN_TABLE = list(_SHOW_COLUMNS_IN_TABLE_DICT.keys())  # Actual list of columns to be shown

    # For the following columns descriptive statistics, such as SD and mean, will be provided in the GUI:
    _COLUMNS_STATISTICS = (
        'sd_g_est_mugal',
        'se_g_est_mugal',
        'g_est_mugal',
    )

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

    def update_view_model(self, lsm_run_index: int, station_name=None, survey_name=None, surveys=None):
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
                # - Keep order of items in `self._SHOW_COLUMNS_IN_TABLE`
                results_stat_df_columns_set = frozenset(results_stat_df.columns)
                table_model_columns = [x for x in self._SHOW_COLUMNS_IN_TABLE if x in results_stat_df_columns_set]

                if results_stat_df is None:
                    self._data = None
                    self._data_column_names = None
                    return

                tmp_filter = results_stat_df['station_name'] != True  # All items are True

                if station_name is not None:
                    tmp_filter = tmp_filter & (results_stat_df['station_name'] == station_name)

                if (survey_name is not None) and (surveys is not None):
                    try:
                        observed_stations_list = surveys[survey_name].observed_stations
                    except KeyError:  # Survey name not available
                        QMessageBox.critical(self.parent(), 'Error!', f'Survey "{survey_name}" not found in campaign!')
                    else:
                        tmp_filter = tmp_filter & results_stat_df['station_name'].isin(observed_stations_list)

                self._data = results_stat_df.loc[tmp_filter, table_model_columns].copy(deep=True)  # Apply filter
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

            value = self._data.iloc[index.row(), index.column()]
            column_name = self._data_column_names[index.column()]

            if role == Qt.DisplayRole or role == Qt.UserRole:
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
                if isinstance(value, int) or isinstance(value, float):
                    # Align right, vertical middle.
                    return Qt.AlignVCenter + Qt.AlignRight

            if role == Qt.BackgroundRole:
                if 'is_datum' in self._data_column_names:
                    if self._data.iloc[index.row(), self._data_column_names.index('is_datum')]:
                        return QtGui.QColor(settings.DATUM_STATION_COLOR[0], settings.DATUM_STATION_COLOR[1], settings.DATUM_STATION_COLOR[2])

        return None

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole and self._data is not None:
            if orientation == Qt.Horizontal:
                return self._SHOW_COLUMNS_IN_TABLE_DICT[str(self._data.columns[section])]
            if orientation == Qt.Vertical:
                return str(self._data.index[section])

    @property
    def mean_g_est(self) -> (float, float):
        """Return the mean of all estimated gravity values in the current model in µGal."""
        if self._data is None:
            return None
        return self._data['g_est_mugal'].mean()

    @property
    def median_g_est(self) -> (float, float):
        """Return the median of all estimated gravity values in the current model in µGal."""
        if self._data is None:
            return None
        return self._data['g_est_mugal'].median()

    @property
    def std_g_est(self) -> (float, float):
        """Return the standard deviation of all estimated gravity values in the current model in µGal."""
        if self._data is None:
            return None
        return self._data['g_est_mugal'].median()

    @property
    def iqr_g_est(self) -> (float, float):
        """Return the IQR of all estimated gravity values in the current model in µGal."""
        if self._data is None:
            return None
        q25 = self._data['g_est_mugal'].quantile(0.25)
        q75 = self._data['g_est_mugal'].quantile(0.75)
        return q75 - q25

    @property
    def mean_sd_est(self) -> (float, float):
        """Return the mean of the standard deviations of all estimated gravity values in the current model in µGal."""
        if self._data is None:
            return None
        return self._data['sd_est_mugal'].mean()

    @property
    def median_sd_est(self) -> (float, float):
        """Return the median of the standard deviations of all estimated gravity values in the current model in µGal."""
        if self._data is None:
            return None
        return self._data['sd_est_mugal'].median()

    @property
    def std_sd_est(self) -> (float, float):
        """Return the standard deviation of the SD's of all estimated gravity values in the current model in µGal."""
        if self._data is None:
            return None
        return self._data['sd_est_mugal'].std()

    @property
    def iqr_sd_est(self) -> (float, float):
        """Return IQR of the standard deviations of all estimated gravity values in the current model in µGal."""
        if self._data is None:
            return None
        q25 = self._data['sd_est_mugal'].quantile(0.25)
        q75 = self._data['sd_est_mugal'].quantile(0.75)
        return q75 - q25

    def get_columns_for_descriptive_statistics(self) -> dict:
        """Returns a dict with names and description of the numerical columns suitable for descriptive statistics.

        Notes
        -----
        If the model data is empty return an empty dict.
        """
        columns_dict = {}
        if self._data_column_names is not None:  # Data model is not empty
            for col_name in self._COLUMNS_STATISTICS:
                # Check if key is a column-name in the current data
                if col_name in self._data_column_names:
                    columns_dict[col_name] = self._SHOW_COLUMNS_IN_TABLE_DICT[col_name]
        return columns_dict

    def get_short_column_description(self, column_name: str) -> str:
        """Returns the short description of the model column."""
        try:
            return self._SHOW_COLUMNS_IN_TABLE_DICT[column_name]
        except AttributeError:
            return ''

    def get_model_data_df(self):
        """Returns the model data dataframe."""
        return self._data