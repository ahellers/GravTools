"""Customized PyQt GUI elements for GravTools.

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

from PyQt5.QtWidgets import QMessageBox, QScrollArea, QLabel, QGridLayout


class ScrollMessageBox(QMessageBox):
   """Custom PyQt message box with scroll area for large texts.

   Source:
   https://stackoverflow.com/questions/47345776/pyqt5-how-to-add-a-scrollbar-to-a-qmessagebox
   """
   def __init__(self, *args, **kwargs):
      QMessageBox.__init__(self, *args, **kwargs)
      chldn = self.children()
      scrll = QScrollArea(self)
      scrll.setWidgetResizable(True)
      grd = self.findChild(QGridLayout)
      lbl = QLabel(chldn[1].text(), self)
      lbl.setWordWrap(True)
      scrll.setWidget(lbl)
      scrll.setMinimumSize (400,200)
      grd.addWidget(scrll,0,1)
      chldn[1].setText('')
      self.exec_()