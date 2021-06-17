"""Support functions for the GUI."""

import pyqtgraph as pg


def get_station_color_dict(stations: list) -> dict:
    """Returns a dict with unique pyqtgraph QColor items for each station in the input list."""
    station_color_dict = {}
    number_of_stations = len(stations)
    for stat_id, station in enumerate(stations, start=0):
        station_color_dict[station] = pg.intColor(stat_id, hues=number_of_stations)
    return station_color_dict


if __name__ == "__main__":
    """Main Program for debugging."""
    stat_list = ['stat1', 'stat2', 'stat3', 'stat4']
    station_color_dict = get_station_color_dict(stat_list)
    print(station_color_dict)
