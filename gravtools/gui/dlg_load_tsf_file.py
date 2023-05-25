"""Dialog for loading time series data from TSF files.

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
import os

from PyQt5.QtWidgets import QDialog, QFileDialog

from gravtools.gui.dialog_load_tsf_file import Ui_DialogLoadTsfFile

class DialogLoadTsfFile(QDialog, Ui_DialogLoadTsfFile):
    """Dialog for loading time series data from TSF files."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)

        self.pushButton_select_file.pressed.connect(self.select_tsf_file)

    def select_tsf_file(self):
        """Launches dialog to select TSF file and get filename."""
        if self.lineEdit_filename:
            initial_path = self.lineEdit_filename.text()
        else:
            try:
                initial_path = self.parent().parent().campaign.output_directory
            except:
                initial_path = os.getcwd()

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        tsf_filename, _ = QFileDialog.getOpenFileName(self,
                                                           'Select TSF file with campaign data',
                                                           initial_path,
                                                           "TSF file (*.TSF)",
                                                           options=options)
        if not tsf_filename:
            return

        # Valid filename:
        self.lineEdit_filename.setText(tsf_filename)
        try:
            survey_name = os.path.basename(tsf_filename).split('.')[0]
            self.lineEdit_survey_name.setText(survey_name)
        except:
            pass

