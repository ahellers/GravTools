"""Misc. methods for the GUI.

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

import pyqtgraph as pg
from PyQt5.QtCore import Qt


def get_station_color_dict(stations: list) -> dict:
    """Returns a dict with unique pyqtgraph QColor items for each station in the input list."""
    station_color_dict = {}
    number_of_stations = len(stations)
    for stat_id, station in enumerate(stations, start=0):
        station_color_dict[station] = pg.intColor(stat_id, hues=number_of_stations)
    return station_color_dict


def checked_state_to_bool(checked_state) -> bool:
    """Converts Qt checked states to boolean values."""
    if checked_state == Qt.Checked or checked_state == Qt.PartiallyChecked:
        return True
    elif checked_state == Qt.Unchecked:
        return False
    else:
        raise AttributeError('Invalid input argument!')


if __name__ == "__main__":
    """Main Program for debugging."""
    stat_list = ['stat1', 'stat2', 'stat3', 'stat4']
    station_color_dict = get_station_color_dict(stat_list)
    print(station_color_dict)



