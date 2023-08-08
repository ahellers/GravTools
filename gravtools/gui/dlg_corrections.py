"""Dialog for managing observation corrections.

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

from PyQt5.QtWidgets import QDialog, QMessageBox
from PyQt5 import QtCore

from gravtools.gui.dialog_corrections import Ui_Dialog_corrections
from gravtools import settings

class DialogCorrections(QDialog, Ui_Dialog_corrections):
    """Dialog to select and apply observation corrections."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)
        # connect signals and slots:

        # Other inits:

        # Set up combo box for interpolation method selection and select default method:
        for idx, (method_name, tooltip) in enumerate(settings.SCIPY_INTERP1_INTERPOLATION_METHODS.items()):
            self.comboBox_tides_interpolation_method.addItem(method_name)
            self.comboBox_tides_interpolation_method.setItemData(idx, tooltip, QtCore.Qt.ToolTipRole)
            if method_name == settings.SCIPY_INTERP1_INTERPOLATION_DEFAULT_METHOD:
                self.comboBox_tides_interpolation_method.setCurrentIndex(idx)

        self.doubleSpinBox_atm_pres_admittance.setValue(settings.ATM_PRES_CORRECTION_ADMITTANCE_DEFAULT)
