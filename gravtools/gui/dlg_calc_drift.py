"""Dialog for the calculation of the instrumental drift.

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
import os
import pytz
import datetime as dt


from PyQt5.QtWidgets import QDialog, QMessageBox
from PyQt5 import QtCore
import pyqtgraph as pg
import numpy as np

from gravtools.gui.dialog_calc_drift import Ui_Dialog_calculate_drift
from gravtools import settings
from gravtools.gui.gui_misc import get_station_color_dict
from gravtools.models.misc import unique_ordered_list

class TimeAxisItem(pg.AxisItem):
    """"Needed to handle the x-axes tags representing date and time.
    From: https://stackoverflow.com/questions/49046931/how-can-i-use-dateaxisitem-of-pyqtgraph

    Notes
    -----
    The timestamps need to refer to UTC!"
    """
    def tickStrings(self, values, scale, spacing) -> str:
        """Handles the x-axes tags representing date and time."""
        return [dt.datetime.fromtimestamp(value, tz=pytz.utc).strftime(settings.Y_TICK_DATETIME_FORMAT) for value in values]

class DialogCalcDrift(QDialog, Ui_Dialog_calculate_drift):
    """Dialog for drift calibration tests."""

    def __init__(self, parent=None):
        super().__init__(parent)
        # Run the .setupUi() method to show the GUI
        self.setupUi(self)

        # Other initial settings:
        self.set_up_drift_plot_widget()

        # Connect signals and slots
        self.comboBox_survey.currentIndexChanged.connect(self.update_station_combobox)
        self.comboBox_station.currentIndexChanged.connect(self.calc_and_plot)
        self.spinBox_degree.valueChanged.connect(self.calc_and_plot)
        self.checkBox_multiple_setups_only.stateChanged.connect(self.update_station_combobox)
        self.checkBox_active_obs_only.stateChanged.connect(self.calc_and_plot)
        self.spinBox_min_no_obs.valueChanged.connect(self.calc_and_plot)

        # Set up combo box for selecting the drift calculation method:
        for idx, (method_name, tooltip) in enumerate(settings.DRIFT_CALC_METHODS.items()):
            self.comboBox_method.addItem(method_name)
            self.comboBox_method.setItemData(idx, tooltip, QtCore.Qt.ToolTipRole)
            if method_name == settings.DRIFT_CALC_DEFAULT_METHOD:
                self.comboBox_method.setCurrentIndex(idx)

        self.station_color_dict = {}

    def invoke(self):
        """Invoke dialog."""
        # Get current campaign data:
        self.reset_gui()
        if self.parent().campaign is None:
            pass
        else:
            # Update GUI items based on current campaign data
            self.update_survey_and_station_combobox()
        return self.exec()

    def update_survey_and_station_combobox(self):
        """Update GUI items based on the current campaign data."""
        # Invoke only when executing the dialog.
        # - Silence the related signals meanwhile.

        if self.parent().campaign is None:
            self.reset_gui()
            return
        survey_names = self.parent().campaign.survey_names
        if not survey_names:  # No surveys in campaign
            self.reset_gui()
            return
        selected_survey = self.comboBox_survey.currentText()
        self.comboBox_survey.blockSignals(True)
        for survey_name in survey_names:
            self.comboBox_survey.addItem(survey_name)
        # Select previously selected item, if available!
        if selected_survey:  # Not empty
            idx = self.comboBox_survey.findText(selected_survey)
            if idx != -1:  # Item found!
                self.comboBox_survey.setCurrentIndex(idx)
        self.comboBox_survey.blockSignals(False)
        self.update_station_combobox()

    def update_station_combobox(self):
        """Update the station selection combo box based on the currently selected survey."""
        current_survey_name = self.comboBox_survey.currentText()
        current_survey = self.parent().campaign.surveys[current_survey_name]

        station_names = current_survey.observed_stations
        self.station_color_dict = get_station_color_dict(station_names, randomize=True)
        if self.checkBox_multiple_setups_only.isChecked():
            station_names_filtered = []
            for station in station_names:
                tmp_filter = current_survey.obs_df['station_name'] == station
                if len(current_survey.obs_df.loc[tmp_filter, 'setup_id'].unique()) > 1:
                    station_names_filtered.append(station)
            station_names = station_names_filtered
        station_names = ['All'] + station_names

        selected_station = self.comboBox_station.currentText()
        self.comboBox_station.blockSignals(True)
        self.comboBox_station.clear()
        for station_name in station_names:
            self.comboBox_station.addItem(station_name)
        # Select previously selected item, if available!
        if selected_station:  # Not empty
            idx = self.comboBox_station.findText(selected_station)
            if idx != -1:  # Item found!
                self.comboBox_station.setCurrentIndex(idx)
            else:
                self.comboBox_station.setCurrentIndex(0)  # First item = 'All'
        self.comboBox_station.blockSignals(False)
        self.comboBox_station.currentIndexChanged.emit(self.comboBox_station.currentIndex())

    def reset_gui(self):
        """Reset the GUI dialog and empty all items and plots."""
        # print('Reset')
        self.comboBox_survey.clear()
        self.comboBox_station.clear()
        # Plot:
        self.drift_calib_plot.clear()
        self.drift_calib_plot.legend.clear()
        self.drift_calib_plot.setTitle('')

    def calc_and_plot(self):
        """Calculate drift parameters and plot the current data."""
        # Invoke via signal whenever a relevant selection in the GUI dialog changed (send signal!)
        # - Make sure that this method is not executed several times in a row... block signals if required!

        print('Update drift plot')
        # Reset:
        self.drift_calib_plot.clear()
        self.drift_calib_plot.legend.clear()
        self.drift_calib_plot.setTitle('')
        # Get data:
        current_survey_name = self.comboBox_survey.currentText()
        survey = self.parent().campaign.surveys[current_survey_name]

        # Get list of selected stations:
        selected_station = self.comboBox_station.currentText()
        if selected_station == 'All':
            stations = [self.comboBox_station.itemText(i) for i in range(self.comboBox_station.count())]
            stations = stations[1:]  # Drop the first items ("All")
        else:
            stations = [self.comboBox_station.currentText()]

        if not stations:
            return

        # Calc. drift:
        method = self.comboBox_method.currentText()
        poly_coefficients ={}
        if method == 'numpy.polyfit':
            for station in stations:
                poly_coefficients[station] = survey.calc_drift_at_station_polyfit(station,
                                                                                  degree=self.spinBox_degree.value(),
                                                                                  obs_type='reduced',
                                                                                  active_only=self.checkBox_active_obs_only.isChecked(),
                                                                                  min_number_obs=self.spinBox_min_no_obs.value())
        else:
            raise RuntimeError(f'Unknown drift calculation method: {method}')

        # ### Plot data:
        # - Give stations different colours and show legend accordingly
        # - Mark/highlight non-active observations
        # self.station_color_dict = get_station_color_dict(stations, randomize=True)
        current_survey = self.parent().campaign.surveys[current_survey_name]
        obs_df = current_survey.obs_df.copy(deep=True)
        if obs_df is None:
            QMessageBox.warning(self, 'Warning!', f'No observation data available for survey {current_survey_name}')
            self.reset_gui()
            return

        # Timeline for drift plot
        # - Hours since survey start
        delta_t_min_h = 0.0 # = 0
        delta_t_max_h = (current_survey.end_time - current_survey.start_time) / np.timedelta64(1, 's') / 3600
        delta_t_h = np.linspace(delta_t_min_h, delta_t_max_h, settings.DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION)
        # - Unix time
        t_ref_unix = dt.datetime.fromisoformat('1970-01-01T00:00:00')
        t_ref_unix = pytz.utc.localize(t_ref_unix)
        delta_t_unix_min = (current_survey.start_time - t_ref_unix) / np.timedelta64(1, 's')
        delta_t_unix_max = (current_survey.end_time - t_ref_unix) / np.timedelta64(1, 's')
        delta_t_unix = np.linspace(delta_t_unix_min, delta_t_unix_max,
                                     settings.DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION)

        # Get gravity offset for plotting:
        tmp_filter = obs_df['station_name'].isin(stations)
        g_offset_mugal = round(obs_df.loc[tmp_filter, 'g_red_mugal'].mean() / 1000) * 1000  # integer mGal

        for station in stations:
            brush_color = self.station_color_dict[station]

            # Plot observations:
            tmp_filter = obs_df['station_name'] == station
            if self.checkBox_active_obs_only.isChecked():
                tmp_filter = tmp_filter & obs_df['keep_obs']
            obs_df_short = obs_df.loc[tmp_filter].copy(deep=True)
            obs_df_short['t_ref_unix'] = (obs_df_short['obs_epoch'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
            scatter = pg.ScatterPlotItem()
            spots = []
            # - prep. data for scatterplot:
            for index, row in obs_df_short.iterrows():
                if row['keep_obs']:
                    brush_color_scatter = brush_color
                    symbol_scatter = 'o'
                else:
                    brush_color_scatter = brush_color
                    symbol_scatter = 'x'
                spot_dic = {'pos': (row['t_ref_unix'], row['g_red_mugal'] - g_offset_mugal),
                    'size': 10,
                    'pen': {'color': 'k',
                            'width': 1},
                    'brush': brush_color_scatter,
                    'symbol': symbol_scatter}
                spots.append(spot_dic)
            scatter.addPoints(spots)
            self.drift_calib_plot.addItem(scatter)

            # Plot drift function:
            coeff_list = poly_coefficients[station]
            if coeff_list is None:  # Not enough ob. available to fit a polynomial
                legend_str = station
            else:
                yy_mugal = np.polyval(coeff_list, delta_t_h)
                pen = pg.mkPen(color='k', width=3)
                self.drift_calib_plot.plot(delta_t_unix, yy_mugal - g_offset_mugal,
                                           # name=f'{station}',
                                           pen=pen, symbol='d', symbolSize=4, symbolBrush=brush_color)
                if len(coeff_list) == 2:  # deg=1
                    legend_str = f'{station}: {coeff_list[0]/1000*24:0.4f} mGal/d'
                elif len(coeff_list) == 3:  # deg=2
                    legend_str = f'{station}: {coeff_list[0]/1000*24:0.4f} mGal/d², {coeff_list[1]/1000*24:0.4f} mGal/d'
                elif len(coeff_list) == 4:  # deg=3
                    legend_str = f'{station}: {coeff_list[0]/1000*24:0.4f} mGal/d³, {coeff_list[1]/1000*24:0.4f} mGal/d², {coeff_list[2]/1000*24:0.4f} mGal/d'
                else:
                    legend_str = station


            # Add item to legend
            s_item_tmp = pg.ScatterPlotItem()
            s_item_tmp.setBrush(brush_color)
            s_item_tmp.setPen({'color': 'k', 'width': 1})
            s_item_tmp.setSize(10)
            self.drift_calib_plot.legend.addItem(s_item_tmp, legend_str)



        # Adjust plot window:
        self.drift_calib_plot.setLabel(axis='left', text=f'g [µGal] + {g_offset_mugal / 1000:.1f} mGal')
        # self.drift_calib_plot.setTitle(f'Drift function w.r.t. setup observations')
        self.drift_calib_plot.autoRange()

            # legend
            # # Add residuals to legend:
            # # - https://pyqtgraph.readthedocs.io/en/latest/graphicsItems/legenditem.html
            # for station, color in self.station_colors_dict_results.items():
            #     if (stations is None) or (station in stations):
            #         s_item_tmp = pg.ScatterPlotItem()
            #         s_item_tmp.setBrush(color)
            #         s_item_tmp.setPen({'color': settings.VG_PLOT_SCATTER_PLOT_PEN_COLOR,
            #                            'width': settings.VG_PLOT_SCATTER_PLOT_PEN_WIDTH})
            #         s_item_tmp.setSize(settings.VG_PLOT_SCATTER_PLOT_SYMBOL_SIZE)
            #         self.vg_plot.legend.addItem(s_item_tmp, res_plot_legend_str + f' ({station})')



    def set_up_drift_plot_widget(self):
        """Set up `self.graphicsLayoutWidget_drift_plot` widget."""
        self.glw_drift_plot = self.graphicsLayoutWidget_drift_plot
        self.glw_drift_plot.setBackground('w')  # white background color
        self.drift_calib_plot = self.glw_drift_plot.addPlot(0, 0, name='ts_plot',
                                                axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        # self.drift_calib_plot.setLabel(axis='bottom', text='Time')  # Takes too much space
        self.drift_calib_plot.addLegend()
        self.drift_calib_plot.showGrid(x=True, y=True)
