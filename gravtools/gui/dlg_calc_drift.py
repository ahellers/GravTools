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
import pytz
import datetime as dt
import os


from PyQt5.QtWidgets import QDialog, QMessageBox
from PyQt5 import QtCore
import pyqtgraph as pg
import numpy as np
import pandas as pd

from gravtools.gui.dialog_calc_drift import Ui_Dialog_calculate_drift
from gravtools import settings
from gravtools.gui.gui_misc import get_station_color_dict
from gravtools import __version__

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
        self.checkBox_remove_g_offset.stateChanged.connect(self.calc_and_plot)

        # Set up combo box for selecting the drift calculation method:
        for idx, (method_name, tooltip) in enumerate(settings.DRIFT_CALC_METHODS.items()):
            self.comboBox_method.addItem(method_name)
            self.comboBox_method.setItemData(idx, tooltip, QtCore.Qt.ToolTipRole)
            if method_name == settings.DRIFT_CALC_DEFAULT_METHOD:
                self.comboBox_method.setCurrentIndex(idx)

        # Init. attributes:
        self.station_color_dict = {}
        self.drift_poly_coeff = {}
        self.drift_start_time = None

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
        # if not current_survey_name:
        #     return
        current_survey = self.parent().campaign.surveys[current_survey_name]

        station_names = current_survey.observed_stations
        self.station_color_dict = get_station_color_dict(station_names, randomize=True)
        if self.checkBox_multiple_setups_only.isChecked():
            station_names_filtered = []
            for station in station_names:
                tmp_filter = current_survey.obs_df['station_name'] == station
                if self.checkBox_active_obs_only.isChecked():
                    tmp_filter = tmp_filter & current_survey.obs_df['keep_obs']
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
        # Gui items:
        self.comboBox_survey.blockSignals(True)
        self.comboBox_survey.clear()
        self.comboBox_survey.blockSignals(False)
        self.comboBox_station.blockSignals(True)
        self.comboBox_station.clear()
        self.comboBox_station.blockSignals(False)
        # Plot:
        self.drift_calib_plot.clear()
        self.drift_calib_plot.legend.clear()
        self.drift_calib_plot.setTitle('')
        # Attributes:
        self.drift_poly_coeff = {}
        self.drift_start_time = None

    def calc_and_plot(self):
        """Calculate drift parameters and plot the current data."""
        # Invoke via signal whenever a relevant selection in the GUI dialog changed (send signal!)
        # - Make sure that this method is not executed several times in a row... block signals if required!

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
        self.drift_poly_coeff = {}
        self.drift_start_time = None
        if method == 'numpy.polyfit':
            for station in stations:
                self.drift_poly_coeff[station] = survey.calc_drift_at_station_polyfit(station,
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
        if self.checkBox_remove_g_offset.isChecked():
            tmp_filter = obs_df['station_name'] == stations[0]
        else:
            tmp_filter = obs_df['station_name'].isin(stations)
        g_offset_mugal = round(obs_df.loc[tmp_filter, 'g_red_mugal'].mean() / 1000) * 1000  # integer mGal

        for station in stations:
            brush_color = self.station_color_dict[station]

            # Plot observations:
            tmp_filter = obs_df['station_name'] == station
            if self.checkBox_active_obs_only.isChecked():
                tmp_filter = tmp_filter & obs_df['keep_obs']
            obs_df_short = obs_df.loc[tmp_filter].copy(deep=True)
            # - Remove g offset w.r.t. the first station in the list:
            if self.checkBox_remove_g_offset.isChecked():
                if station != stations[0]:  # Apply offset w.r.t. the first station
                    g_mean_stat = obs_df_short['g_red_mugal'].mean()
                    g_offset_stat_mugal = g_mean_first_stat - g_mean_stat
                    # obs_df_short['g_red_mugal'] = obs_df_short['g_red_mugal'] + g_offset_stat
                elif station == stations[0]:
                    g_mean_first_stat = obs_df_short['g_red_mugal'].mean()
                    g_offset_stat_mugal = 0.0
            else:
                g_offset_stat_mugal = 0.0
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
                spot_dic = {'pos': (row['t_ref_unix'], row['g_red_mugal'] - g_offset_mugal + g_offset_stat_mugal),
                    'size': 10,
                    'pen': {'color': 'k',
                            'width': 1},
                    'brush': brush_color_scatter,
                    'symbol': symbol_scatter}
                spots.append(spot_dic)
            scatter.addPoints(spots)
            self.drift_calib_plot.addItem(scatter)

            # Plot drift function:
            coeff_list = self.drift_poly_coeff[station]
            if coeff_list is None:  # Not enough ob. available to fit a polynomial
                legend_str = station
            else:
                yy_mugal = np.polyval(coeff_list, delta_t_h)
                pen = pg.mkPen(color='k', width=3)
                self.drift_calib_plot.plot(delta_t_unix, yy_mugal - g_offset_mugal + g_offset_stat_mugal,
                                           # name=f'{station}',
                                           pen=pen, symbol='d', symbolSize=4, symbolBrush=brush_color)
                if len(coeff_list) == 2:  # deg=1
                    legend_str = f'{station}: {coeff_list[0]/1000*24:0.4f} mGal/d'
                elif len(coeff_list) == 3:  # deg=2
                    legend_str = f'{station}: {coeff_list[0]/1000*(24**2):0.4f} mGal/d², {coeff_list[1]/1000*24:0.4f} mGal/d'
                elif len(coeff_list) == 4:  # deg=3
                    legend_str = f'{station}: {coeff_list[0]/1000*(24**3):0.4f} mGal/d³, {coeff_list[1]/1000*(24**2):0.4f} mGal/d², {coeff_list[2]/1000*24:0.4f} mGal/d'
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

    def set_up_drift_plot_widget(self):
        """Set up `self.graphicsLayoutWidget_drift_plot` widget."""
        self.glw_drift_plot = self.graphicsLayoutWidget_drift_plot
        self.glw_drift_plot.setBackground('w')  # white background color
        self.drift_calib_plot = self.glw_drift_plot.addPlot(0, 0, name='ts_plot',
                                                axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        # self.drift_calib_plot.setLabel(axis='bottom', text='Time')  # Takes too much space
        self.drift_calib_plot.addLegend()
        self.drift_calib_plot.showGrid(x=True, y=True)

    def save(self):
        """Save the drift results as log file and plot."""
        campaign_dir = self.parent().campaign.output_directory
        filename = self.comboBox_survey.currentText() + '_drift.txt'
        filename = os.path.join(self.parent().campaign.output_directory, filename)
        self.save_drift_file(filename)
        filename = self.comboBox_survey.currentText() + '_drift.png'
        filename = os.path.join(self.parent().campaign.output_directory, filename)
        self.save_drift_plot(filename)

    def save_drift_file(self, filename, verbose=True):
        """Save the drift calculation results to a file."""
        if verbose:
            print(f'Save drift results as {filename}')

        # Remove None entries from dict:
        poly_ceoffs = {}
        for key, value in self.drift_poly_coeff.items():
            if value is not None:
                poly_ceoffs[key] = value
        if len(poly_ceoffs) == 0:
            raise RuntimeError('No valid drift polynomials available for export to file.')

        poly_function_str = f'f(t) = a0(t0)'
        colnames = ['a0 [mGal]']
        for i in range(self.spinBox_degree.value()):
            poly_function_str = poly_function_str + f' + a{i + 1}*(t-t0)^{i + 1}'
            colnames.append(f'a{i + 1} [mGal/d^{i + 1}]')
        poly_df = pd.DataFrame.from_dict(poly_ceoffs, orient='index')
        poly_df = poly_df[poly_df.columns[::-1]]
        for colnum in range(poly_df.shape[1]):  # muGal/h => mGal/day
            poly_df.iloc[:, colnum] = poly_df.iloc[:, colnum] / 1000 * (24**colnum)
        poly_df.columns = colnames
        poly_df.loc[:, 'Station'] = poly_df.index
        colnames = ['Station'] + colnames
        poly_df = poly_df[colnames]
        poly_df_str = poly_df.to_string(index=False, float_format=lambda x: '{:.4f}'.format(x))
        start_time = self.parent().campaign.surveys[self.comboBox_survey.currentText()].start_time
        try:
            start_time_tz = start_time.tzname()
        except:
            start_time_tz = ''

        tmp_str = f'### Drift results for survey {self.comboBox_survey.currentText()} ###\n'
        tmp_str = tmp_str + f'Calculated by GravTools {__version__}\n'
        tmp_str = tmp_str + f'File created: {dt.datetime.now().strftime("%Y-%m-%d, %H:%M:%S")}\n'
        tmp_str = tmp_str + f'Method: {settings.DRIFT_CALC_METHODS[self.comboBox_method.currentText()]}\n'
        tmp_str = tmp_str + f'Polynomial degree: {self.spinBox_degree.value()}\n'
        tmp_str = tmp_str + f'Min. number of obs.: {self.spinBox_min_no_obs.value()}\n'
        tmp_str = tmp_str + f'Active obs. only: {self.checkBox_active_obs_only.isChecked()}\n'
        tmp_str = tmp_str + f'Multiple setups only: {self.checkBox_multiple_setups_only.isChecked()}\n'
        tmp_str = tmp_str + '\n'
        tmp_str = tmp_str + poly_function_str + '\n'
        tmp_str = tmp_str + f'... with [t]=days and [g]=mGal\n'
        if start_time_tz:
            tmp_str = tmp_str + f'... with t0 = {start_time.isoformat()} ({start_time_tz})\n'
        else:
            tmp_str = tmp_str + f'... with t0 = {start_time.isoformat()}\n'
        tmp_str = tmp_str + f'\n'
        tmp_str = tmp_str + poly_df_str + '\n'
        tmp_str = tmp_str + f'\n'

        if verbose:
            print(tmp_str)

        with open(filename, 'w') as f:
            f.write(tmp_str)

    def save_drift_plot(self, filename, verbose=True):
        """Save the current drift plot as png file."""
        # exporter = pg.exporters.ImageExporter(self.graphicsLayoutWidget_results_drift_plot.scene())
        exporter = pg.exporters.ImageExporter(self.graphicsLayoutWidget_drift_plot.scene())
        flag_export_successful = exporter.export(filename)
        if verbose and flag_export_successful:
            print(f'Save drift plot as {filename}')

