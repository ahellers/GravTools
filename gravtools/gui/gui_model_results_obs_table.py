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
import pandas as pd


NONE_REPRESENTATION_IN_TABLE_VIEW = ''  # Representation of None values in table views in the GUI


class ResultsObservationModel(QAbstractTableModel):
    """Model for displaying the observations-related results."""

    # Number of decimal pla
    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {
        'g_diff_mugal': 1,
        'sd_g_diff_mugal': 1,
        'sd_g_diff_est_mugal': 1,
        'v_diff_mugal': 1,
        'sd_v_diff_mugal': 1,
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
        'sd_v_obs_est_mugal': 1,  # LSM_non_diff
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
        'ref_epoch_dt': 'Epoch UTC',  # LSM_non_diff, LSM_diff
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
        'sd_v_obs_est_mugal': 'SD_v [µGal]',  # LSM_non_diff
        'sd_v_diff_mugal': 'SD_v [µGal]',  # LSM_non_diff
        'w_diff_mugal': 'w',  # LSM_diff
        'w_obs_est_mugal': 'w',  # LSM_non_diff
        'r_diff_obs': 'r',  # LSM_diff
        'r_obs_est': 'r',  # LSM_non_diff
        'tau_test_result': 'Outlier test',  # LSM_non_diff, LSM_diff
        'corr_drift_mugal': 'Drift corr. [µGal]',  # MLR: Estimated drift correction
        'abw_mugal': 'abw [µGal]',  # MLR: Difference between drift-corrected instrument reading and the estimated station gravity
    }
    _SHOW_COLUMNS_IN_TABLE = list(_SHOW_COLUMNS_IN_TABLE_DICT.keys())  # Actual list of columns to be shown

    # List of columns that are shown in the simple mode of the GUI
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
        'w_diff_mugal': 'w',  # LSM_diff
        'r_diff_obs': 'Redundancy components []',  # LSM_diff
        'sd_v_diff_mugal': 'SD of post-fit residuals [µGal]',  # LSM_diff
        'sd_g_diff_est_mugal': 'A posteriori SD of diff. obs. [µGal]',  # LSM_diff
        'g_diff_mugal': 'Differential observation [µGal]',  # LSM_diff
        'sd_g_diff_mugal': 'SD of differential observation [µGal]',  # LSM_diff
        # MLR:
        'g_mugal': 'Instrument reading [µGal]',  # MLR
        'sd_g_mugal': 'Standard deviation of the instrument reading [µGal]',  # MLR
        'corr_drift_mugal': 'Estimated drift correction [µGal]',  # MLR
        'abw_mugal': 'Drift-corrected reading minus estimated station gravity (not absolute) [µGal]',  # MLR
        # LSM_non_diff:
        'v_obs_est_mugal': 'Post-fit residuals [µGal]',  # LSM_non_diff
        'w_obs_est_mugal': 'Standardized post-fit residuals []',  # LSM_non_diff
        'r_obs_est': 'Redundancy components []',  # LSM_non_diff
        'sd_v_obs_est_mugal': 'SD of post-fit residuals [µGal]',  # LSM_non_diff
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
        flag_error_init = False
        if lsm_run_index == -1:  # No data available => Invalid index => Reset model data
            flag_error_init = True
        else:
            try:
                results_obs_df = self._lsm_runs[lsm_run_index].get_results_obs_df
            except KeyError:
                QMessageBox.critical(self.parent(), 'Error!', f'LSM run with index "{lsm_run_index}" not found!')
            else:
                if results_obs_df is None:  # E.g. no results in campaign yet
                    flag_error_init = True
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
                        except:  # Just in case any problem occurs...
                            self._data = results_obs_df.loc[tmp_filter, :].copy(deep=True)  # Show all columns
                    else:  # No filter
                        try:
                            self._data = results_obs_df.loc[:, table_model_columns].copy(deep=True)
                        except:  # Just in case any problem occurs...
                            self._data = results_obs_df.copy(deep=True)
                    self._data_column_names = self._data.columns.to_list()  # Show all columns
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
