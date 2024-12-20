"""Dialog for exporting results.

Copyright (C) 2023  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

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

from PyQt5.QtWidgets import QDialog
from PyQt5.QtCore import Qt, pyqtSlot

from gravtools.gui.dialog_export_results import Ui_Dialog_export_results

class DialogExportResults(QDialog, Ui_Dialog_export_results):
    """Dialog to define the estimation settings."""

    def __init__(self, campaign, _has_geopandas=False, parent=None):
        super().__init__(parent)
        # Get data:
        self._lsm_runs = campaign.lsm_runs
        self.flag_observation_data_available = len(campaign.surveys) > 0
        self.flag_lsm_runs_available = len(campaign.lsm_run_times) > 0
        # Run the .setupUi() method to show the GUI
        self.setupUi(self)
        # Populate the combo box to select an LSM run (enable/disable groupBoxes accordingly):
        self.comboBox_select_lsm_run.clear()
        self.comboBox_select_lsm_run.addItems(['Current observation selection (no LSM run)'] + campaign.lsm_run_times)

        # Set the lineEdit with the export path:
        self.label_export_path_show.setText(campaign.output_directory)

        # Set initial item in comboBox (last entry/LSM run):
        idx = self.comboBox_select_lsm_run.count() - 1  # Index of last lsm run
        self.comboBox_select_lsm_run.setCurrentIndex(idx)

        # Select LSM run, etc.
        if self.flag_lsm_runs_available:
            lsm_run_idx = idx - 1
        else:
            lsm_run_idx = -1

        # Enable/disable GUI items and add lsm run comment:
        self.write_lsm_run_comment_to_gui(lsm_run_idx)
        self.enable_gui_widgets_based_on_lsm_run_selection(lsm_run_idx)

        # connect signals and slots:
        self.comboBox_select_lsm_run.currentIndexChanged.connect(self.on_comboBox_select_lsm_run_current_index_changed)

        # Optional dependency for GIS data export:
        if _has_geopandas:
            self.groupBox_gis.setEnabled(True)
        else:
            self.groupBox_gis.setEnabled(False)
            self.checkBox_gis_write_shapefile.setChecked(Qt.Unchecked)

    @pyqtSlot(int)
    def on_comboBox_select_lsm_run_current_index_changed(self, index: int):
        """Invoked whenever the index of the selected item in the combobox changed."""
        lsm_run_idx = index - 1
        self.write_lsm_run_comment_to_gui(lsm_run_idx)
        self.enable_gui_widgets_based_on_lsm_run_selection(lsm_run_idx)

    def write_lsm_run_comment_to_gui(self, lsm_run_idx):
        """Writes the lsm run comment of the selected lsm run to the GUI line edit."""
        if lsm_run_idx > -1:
            self.label_export_comment_show.setText(self.get_lsm_run_comment(lsm_run_idx))
        else:
            self.label_export_comment_show.setText('')

    def get_lsm_run_comment(self, lsm_run_idx: int):
        """Returns the lsm run comment of the run with the specified index."""
        try:
            return self._lsm_runs[lsm_run_idx].comment
        except:
            return ''

    def enable_gui_widgets_based_on_lsm_run_selection(self, lsm_run_idx: int):
        """Enable or disable GUI elements based on the lsm run selection."""

        flag_lsm_run_selected = self.flag_lsm_runs_available and (lsm_run_idx > -1)

        if flag_lsm_run_selected:
            self.groupBox_other_files.setEnabled(True)
            self.groupBox_nsb_file.setEnabled(True)
            self.groupBox_observation_list.setEnabled(True)
            if _has_geopandas:
                self.groupBox_gis.setEnabled(True)
        else:
            self.groupBox_gis.setEnabled(False)
            self.groupBox_other_files.setEnabled(False)
            self.groupBox_nsb_file.setEnabled(False)
            self.groupBox_observation_list.setEnabled(False)
            if self.flag_observation_data_available:
                self.groupBox_observation_list.setEnabled(True)
        if flag_lsm_run_selected or self.flag_observation_data_available:
            self.buttonBox.buttons()[0].setEnabled(True)  # OK button in buttonBox
        else:
            self.buttonBox.buttons()[0].setEnabled(False)  # OK button in buttonBox

        if flag_lsm_run_selected:  # => lsm_run_idx >= 0
            lsm_method = self._lsm_runs[lsm_run_idx].lsm_method
            if lsm_method == 'LSM_diff' or lsm_method == 'LSM_non_diff' or lsm_method == 'MLR_BEV':  # network adjust.
                self.groupBox_nsb_file.setEnabled(True)
                self.checkBox_save_vg_plot_png.setEnabled(False)
            elif lsm_method == 'VG_LSM_nondiff':  # VG estimation
                self.groupBox_nsb_file.setEnabled(False)
                self.checkBox_save_vg_plot_png.setEnabled(True)
            else:
                raise AssertionError(f'Invalid LSM method: {lsm_method}!')