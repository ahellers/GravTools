"""Dialog for loading CG6 observation files.

Copyright (C) 2024  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

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
from PyQt5.QtWidgets import QDialog, QFileDialog, QMessageBox
from PyQt5.QtCore import Qt


from gravtools.gui.dialog_load_cg6_obs_files import Ui_DialogLoadCg6ObservationFiles
from gravtools import settings


class DialogLoadCg6ObservationFiles(QDialog, Ui_DialogLoadCg6ObservationFiles):
    """Dialog for loading CG6 observation files in different formats."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)

        # Fill comboboxes:
        for key, item in settings.CG6_SURVEY_DATA_SOURCE_TYPES_SHORT.items():
            self.comboBox_format_selection.addItem(item, userData=key)
            idx = self.comboBox_format_selection.findData(key)
            tooltip = settings.CG6_SURVEY_DATA_SOURCE_TYPES[key]
            self.comboBox_format_selection.setItemData(idx, tooltip, Qt.ToolTipRole)

        for key, item in settings.CG6_OBS_FILE_LOCATION_TYPES.items():
            self.comboBox_location_type.addItem(item['label'], userData=key)
            idx = self.comboBox_location_type.findData(key)
            self.comboBox_location_type.setItemData(idx, item['tooltip'], Qt.ToolTipRole)

        for key, item in settings.CG6_OBS_FILE_ERROR_TYPES.items():
            self.comboBox_uncertainty_type.addItem(item['label'], userData=key)
            idx = self.comboBox_uncertainty_type.findData(key)
            self.comboBox_uncertainty_type.setItemData(idx, item['tooltip'], Qt.ToolTipRole)

        for key, item in settings.CG6_PRESSURE_IN_COLUMN.items():
            self.comboBox_pres_in_column.addItem(item['label'], userData=key)
            idx = self.comboBox_pres_in_column.findData(key)
            self.comboBox_pres_in_column.setItemData(idx, item['tooltip'], Qt.ToolTipRole)

        for key, item in settings.CG6_DHB_IN_COLUMN.items():
            self.comboBox_dhb_in_column.addItem(item['label'], userData=key)
            idx = self.comboBox_dhb_in_column.findData(key)
            self.comboBox_dhb_in_column.setItemData(idx, item['tooltip'], Qt.ToolTipRole)

        # Connect signals and slots:
        self.comboBox_format_selection.currentIndexChanged.connect(self.on_combobox_current_index_changed)
        self.pushButton_select_files.pressed.connect(self.select_files)
        self.pushButton_clear_list.pressed.connect(self.listWidget_open_files.clear)
        self.comboBox_pres_in_column.currentIndexChanged.connect(self.check_alternative_column_data_selection)
        self.comboBox_dhb_in_column.currentIndexChanged.connect(self.check_alternative_column_data_selection)

        # Set default values for GUI items (AFTER connecting signals and slots):
        idx = self.comboBox_format_selection.findData(settings.CG6_SURVEY_DATA_SOURCE_TYPE_DEFAULT)
        self.comboBox_format_selection.setCurrentIndex(idx)
        idx = self.comboBox_pres_in_column.findData(settings.CG6_PRESSURE_IN_COLUMN_DEFAULT)
        self.comboBox_pres_in_column.setCurrentIndex(idx)
        idx = self.comboBox_dhb_in_column.findData(settings.CG6_DHB_IN_COLUMN_DEFAULT)
        self.comboBox_dhb_in_column.setCurrentIndex(idx)
        idx = self.comboBox_uncertainty_type.findData(settings.CG6_OBS_FILE_ERROR_TYPES_DEFAULT)
        self.comboBox_uncertainty_type.setCurrentIndex(idx)
        idx = self.comboBox_location_type.findData(settings.CG6_OBS_FILE_LOCATION_TYPES_DEFAULT)
        self.comboBox_location_type.setCurrentIndex(idx)

    def check_alternative_column_data_selection(self):
        """Check selection in comboboxes for alternative data in columns and raise a warning."""
        dhb_data = self.comboBox_dhb_in_column.currentData()
        pres_data = self.comboBox_pres_in_column.currentData()
        if (dhb_data == pres_data) and not (dhb_data == 'no_data' and pres_data == 'no_data'):
            QMessageBox.warning(self, 'Warning!', f'Non-unique assignment of alternative data in observation file '
                                                  f'columns: {settings.CG6_DHB_IN_COLUMN[dhb_data]["label"]}.\n\n'
                                                  f'Please change the column selection!')

    def on_combobox_current_index_changed(self, idx):
        """On combobox index change."""
        if idx == -1:  # Invalid index
            return
        format_key = self.comboBox_format_selection.itemData(idx)
        # Update format description label
        description = settings.CG6_SURVEY_DATA_SOURCE_TYPES[format_key]
        self.label_format_description.setText(description)

    def select_files(self):
        """Select observation files using a file selection dialog"""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        obs_file_filenames, _ = QFileDialog.getOpenFileNames(self,
                                                             'Select CG6 observation file(s)',
                                                             self.parent().campaign.output_directory,
                                                             "CG6 observation file (*.DAT *.TXT)",
                                                             options=options)
        if not obs_file_filenames:
            return
        self.listWidget_open_files.clear()
        self.listWidget_open_files.addItems(obs_file_filenames)

    @property
    def file_list(self) -> list:
        """Return a list of all files in the list widget."""
        # self.listWidget_open_file.
        return [self.listWidget_open_files.item(idx).text() for idx in range(self.listWidget_open_files.count())]

    @property
    def file_format(self) -> str:
        """Returns the file format identifier string."""
        return self.comboBox_format_selection.currentData()

    @property
    def location_type(self) -> str:
        """Returns the location type identifier."""
        return self.comboBox_location_type.currentData()

    @property
    def error_type(self) -> str:
        """Returns the error type identifier."""
        return self.comboBox_uncertainty_type.currentData()

    @property
    def pres_in_column(self) -> str:
        """Returns the column identifier for obs_df for pressure data."""
        return self.comboBox_pres_in_column.currentData()

    @property
    def dhb_in_column(self) -> str:
        """Returns the column identifier for height differences dhb."""
        return self.comboBox_dhb_in_column.currentData()
