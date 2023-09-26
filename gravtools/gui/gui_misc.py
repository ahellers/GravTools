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
import random


def get_station_color_dict(stations: list, randomize=False) -> dict:
    """Returns a dict with unique pyqtgraph QColor items for each station in the input list."""
    station_color_dict = {}
    number_of_stations = len(stations)

    stat_id_list = [i for i in range(0, number_of_stations)]
    if randomize:
        random.seed(1)  # Use always the same seed to get reproducibility
        random.shuffle(stat_id_list)

    for stat_id, station in enumerate(stations, start=0):
        stat_id_altered = stat_id_list[stat_id]
        station_color_dict[station] = pg.intColor(stat_id_altered, hues=number_of_stations)
    return station_color_dict


def checked_state_to_bool(checked_state) -> bool:
    """Converts Qt checked states to boolean values."""
    if checked_state == Qt.Checked or checked_state == Qt.PartiallyChecked:
        return True
    elif checked_state == Qt.Unchecked:
        return False
    else:
        raise AttributeError('Invalid input argument!')


def resize_table_view_columns(table_view, n, add_pixel=0):
    """Resize the width of columns of a QTableView based on the header and the content of the first n rows.

    Parameters
    ----------
    table_view: QTableView object
        The QTableView object that should be concerned.
    n : int
        The column width is adjusted to the width of the column headers and the content of the first n rows. Is the
        model contains less than n items, all available items are considered.
    add_pixel : int >= 0
        Number of pixels to be added to the minimum width to improve readability.
    """
    fm = table_view.fontMetrics()
    for col_num in range(0, table_view.model().columnCount()):
        name= table_view.model().headerData(section=col_num, orientation=Qt.Horizontal, role=Qt.DisplayRole)
        max_width = fm.width(name)
        for row_idx in range(0, n):
            data_index = table_view.model().index(row_idx, col_num)
            if not data_index.isValid():
                break
            data = table_view.model().data(data_index)
            width_data = fm.width(data)
            if width_data > max_width:
                max_width = width_data
        table_view.setColumnWidth(col_num, max_width + add_pixel)


if __name__ == "__main__":
    """Main Program for debugging."""
    stat_list = ['stat1', 'stat2', 'stat3', 'stat4']
    station_color_dict = get_station_color_dict(stat_list)
    print(station_color_dict)



