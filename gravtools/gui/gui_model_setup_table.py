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

    # Columns that will be shown in the table view, if available in the data (Also defines the order of columns):
    # - keys: Actual names of the dataframe columns
    # - items: Header names for the Table View Widget
    _SHOW_COLUMNS_IN_TABLE_DICT = {
        'station_name': 'Station',
        'setup_id': 'Setup ID',
        'epoch_dt': 'Epoch UTC',
        'g_mugal': 'g [µGal]',
        'sd_g_mugal': 'SD [µGal]',
        'epoch_unix': 'Epoch Unix',
        'delta_t_h': 'd_t [h]',
        'delta_t_campaign_h': 'd_t camp [h]',
        'sd_setup_mugal': 'SD of obs [µGal]',
        'number_obs': 'Number of obs.',
        'dhf_sensor_m': 'dhf_sensor [m]',
        'linear_scale': 'linear scale',
    }
    _SHOW_COLUMNS_IN_TABLE = list(_SHOW_COLUMNS_IN_TABLE_DICT.keys())  # Actual list of columns to be shown

    # List of columns that are shown in the simple mode of the GUI
    _SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI = [
        'station_name',
        'setup_id',
        'g_mugal',
        'sd_g_mugal',
        'epoch_dt',
    ]

    _DECIMAL_PLACES_PER_FLOAT_COLUMN = {  # key = column name; value = number of decimal places
        'g_mugal': 1,
        'sd_g_mugal': 1,
        'epoch_unix': 1,
        'delta_t_h': 3,
        'delta_t_campaign_h': 3,
        'sd_setup_mugal': 1,
        'dhf_sensor_m': 3,
        'linear_scale': 5,
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
        self._data_column_names = None
        self.flag_gui_simple_mode = False
        self._reference_height_type = ''
        self._tide_correction_type = ''
        self._atm_pres_correction_type = ''
        self._scale_correction_type = ''
        self._oceanload_correction_type = ''

    def reset_model(self):
        """Reset the model to a state without selected data to display."""
        self._data = None  # Observations (or at subset of them) of the survey with the name `self._data_survey_name`
        self._data_survey_name = ''  # Name of the Survey that is currently represented by `self._data`
        self._data_column_names = None
        self._reference_height_type = ''
        self._tide_correction_type = ''
        self._atm_pres_correction_type = ''
        self._scale_correction_type = ''
        self._oceanload_correction_type = ''

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
        """Update the `_data` DataFrame that hold the actual data that is displayed."""
        if survey_name is not None:
            try:
                setup_df = self._surveys[survey_name].setup_df
            except KeyError:
                QMessageBox.critical(self.parent(), 'Error!', f'Survey "{survey_name}" is not available in this campaign.')
            else:
                self._data_survey_name = survey_name
                try:
                    self._reference_height_type = self._surveys[survey_name].setup_reference_height_type
                    self._tide_correction_type = self._surveys[survey_name].setup_tide_correction_type
                except AttributeError:
                    self._reference_height_type = ''
                    self._tide_correction_type = ''
                try:
                    self._atm_pres_correction_type = self._surveys[survey_name].setup_atm_pres_correction_type
                except AttributeError:
                    self._atm_pres_correction_type = ''
                try:
                    self._scale_correction_type = self._surveys[survey_name].setup_scale_correction_type
                except AttributeError:
                    self._scale_correction_type = ''
                try:
                    self._oceanload_correction_type = self._surveys[survey_name].setup_oceanload_correction_type
                except AttributeError:
                    self._oceanload_correction_type = ''
                self.flag_gui_simple_mode = gui_simple_mode

                # Get list of columns to be depicted via the table model:
                # - Keep order of items in `self._SHOW_COLUMNS_IN_TABLE`
                if setup_df is not None:
                    setup_df_columns_set = frozenset(setup_df.columns)
                    table_model_columns = [x for x in self.get_table_columns if x in setup_df_columns_set]

                if setup_id is None:  # No setup ID provided => Take all observations in survey
                    if setup_df is None:
                        self._data = None
                        self._data_column_names = None
                    else:
                        # self._data = setup_df.copy(deep=True)
                        self._data = setup_df.loc[:, table_model_columns].copy(deep=True)
                        self._data_column_names = self._data.columns.to_list()
                else:  # Only take observations of the specified setup
                    if setup_df is None:
                        self._data = None
                        self._data_column_names = None
                    else:
                        # self._data = setup_df[setup_df['setup_id'] == setup_id].copy(deep=True)
                        self._data = setup_df.loc[setup_df['setup_id'] == setup_id, table_model_columns].copy(deep=True)
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

    @property
    def get_table_columns(self):
        """Returns a list with all columns names of the dataframe that should be copied to the view model."""
        if self.flag_gui_simple_mode:
            return [value for value in self._SHOW_COLUMNS_IN_TABLE if value in self._SHOW_COLUMNS_IN_TABLE_SIMPLE_GUI]
        else:
            return self._SHOW_COLUMNS_IN_TABLE

    @property
    def get_ref_heigth_type(self):
        """Returns the reference height type of the setup observations."""
        return self._reference_height_type

    @property
    def get_tidal_corr_type(self):
        """Returns the tidal corrections applied on the setup observations."""
        return self._tide_correction_type

    @property
    def get_atm_pres_corr_type(self):
        """Returns the atmospheric pressure corrections applied on the setup observations."""
        return self._atm_pres_correction_type

    @property
    def get_scale_corr_type(self):
        """Returns the scale corrections applied on the setup observations."""
        return self._scale_correction_type

    @property
    def get_oceanload_corr_type(self):
        """Returns the ocean-loading correction applied on the setup observations."""
        return self._oceanload_correction_type
