"""Context menu handler for the survey table.

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

from PyQt5.QtWidgets import QMenu
from PyQt5 import QtGui


class SurveyTableHandler:
    """Context menu handler for the survey table."""

    def __init__(self, parent=None):
        self.parent = parent

    def right_click(self):
        # bar = self.parent.menuBar()
        top_menu = QMenu(self.parent)

        menu = top_menu.addMenu("Menu")
        config = menu.addMenu("Configuration ...")

        _load = config.addAction("&Load ...")
        _save = config.addAction("&Save ...")

        config.addSeparator()

        config1 = config.addAction("Config1")
        config2 = config.addAction("Config2")
        config3 = config.addAction("Config3")

        action = menu.exec_(QtGui.QCursor.pos())

        if action == _load:
            # do this
            pass
        elif action == _save:
            # do this
            pass
        elif action == config1:
            # do this
            pass
        elif action == config2:
            # do this
            pass
        elif action == config3:
            # do this
            pass