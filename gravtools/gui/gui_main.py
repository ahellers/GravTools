import sys
import os
from PyQt5.QtWidgets import QApplication, QDialog, QMainWindow, QFileDialog, QMessageBox, QTreeWidgetItem, QHeaderView
from PyQt5.QtCore import QDir, QAbstractTableModel, Qt, QSortFilterProxyModel, pyqtSlot, QRegExp, QModelIndex
from PyQt5 import QtGui

import datetime as dt
import pyqtgraph as pg
import numpy as np
import pandas as pd
import pytz

from MainWindow import Ui_MainWindow
from dialog_new_campaign import Ui_Dialog_new_Campaign
from dialog_load_stations import Ui_Dialog_load_stations
from dialog_corrections import Ui_Dialog_corrections
from dialog_autoselection_settings import Ui_Dialog_autoselection_settings
from dialog_estimation_settings import Ui_Dialog_estimation_settings

from gravtools.models.survey import Campaign, Survey, Station
from gravtools import settings
from gui_models import StationTableModel, ObservationTableModel, SetupTableModel, ResultsStationModel, \
    ResultsObservationModel, ResultsDriftModel

DEFAULT_OUTPUT_DIR = os.path.abspath(os.getcwd())  # Current working directory
DEFAULT_CG5_OBS_FILE_PATH = os.path.abspath(os.getcwd())  # Current working directory
IS_VERBOSE = True  # Define, whether screen output is enabled.

MARKER_SYMBOL_ORDER = ('o', 't', 'x', 's', 'star', '+', 'd', 't1', 'p', 't2', 'h', 't3')
MARKER_COLOR_ODER = ('b', 'r', 'g', 'c', 'm', 'y')


def checked_state_to_bool(checked_state) -> bool:
    """Converts Qt checked states to boolean values."""
    if checked_state == Qt.Checked or checked_state == Qt.PartiallyChecked:
        return True
    elif checked_state == Qt.Unchecked:
        return False
    else:
        raise AttributeError('Invalid input argument!')


class TimeAxisItem(pg.AxisItem):
    """"Needed to handle the x-axes tags representing date and time.
    From: https://stackoverflow.com/questions/49046931/how-can-i-use-dateaxisitem-of-pyqtgraph

    Notes
    -----
    The timestamps need to refer to UTC!"
    """

    def tickStrings(self, values, scale, spacing) -> str:
        """Handles the x-axes tags representing date and time."""
        return [dt.datetime.fromtimestamp(value, tz=pytz.utc) for value in values]


class MainWindow(QMainWindow, Ui_MainWindow):
    """Main Window of the application."""

    # General options for the main window:
    BRUSH_ACTIVE_OBS = pg.mkBrush('g')
    BRUSH_INACTIVE_OBS = pg.mkBrush('r')

    def __init__(self):
        """Initializer."""

        # Instance Attributes:
        self.campaign = None

        # GUI:
        super().__init__()
        self.setupUi(self)

        # Connect signals and slots
        self.action_Exit.triggered.connect(self.exit_application)
        self.action_New_Campaign.triggered.connect(self.on_menu_file_new_campaign)
        self.action_Add_Stations.triggered.connect(self.on_menu_file_load_stations)
        self.action_Corrections.triggered.connect(self.on_menu_observations_corrections)
        self.action_Autoselection_settings.triggered.connect(self.on_menu_observations_autoselection_settings)
        self.action_Estimation_settings.triggered.connect(self.on_menu_estimation_settings)
        self.pushButton_obs_apply_autoselect_current_data.pressed.connect(self.on_apply_autoselection)
        self.pushButton_obs_comp_setup_data.pressed.connect(self.on_pushbutton_obs_comp_setup_data)
        self.pushButton_obs_run_estimation.pressed.connect(self.on_pushbutton_obs_run_estimation)
        self.pushButton_results_delete_lsm_run.pressed.connect(self.on_pushbutton_results_delete_lsm_run)
        # self.actionShow_Stations.triggered.connect(self.show_station_data)
        self.action_from_CG5_observation_file.triggered.connect(self.on_menu_file_load_survey_from_cg5_observation_file)
        self.lineEdit_filter_stat_name.textChanged.connect(self.on_lineEdit_filter_stat_name_textChanged)
        self.checkBox_filter_observed_stat_only.stateChanged.connect(self.on_checkBox_filter_observed_stat_only_toggled)
        self.checkBox_obs_plot_setup_data.stateChanged.connect(self.on_checkBox_obs_plot_setup_data_state_changed)
        # Observations tree widget:
        self.treeWidget_observations.itemSelectionChanged.connect(self.on_obs_tree_widget_item_selected)
        self.treeWidget_observations.itemChanged.connect(self.on_tree_widget_item_changed)
        self.checkBox_obs_plot_reduced_observations.clicked.connect(self.on_obs_tree_widget_item_selected)
        self.comboBox_results_lsm_run_selection.currentIndexChanged.connect(
            self.on_comboBox_results_lsm_run_selection_current_index_changed)
        self.comboBox_results_selection_station.currentIndexChanged.connect(
            self.on_comboBox_results_selection_station_current_index_changed)
        self.comboBox_results_selection_survey.currentIndexChanged.connect(
            self.on_comboBox_results_selection_survey_current_index_changed)
        self.comboBox_results_obs_plot_select_data_column.currentIndexChanged.connect(
            self.on_comboBox_results_obs_plot_select_data_column_current_index_changed)

        # Set up GUI items and widgets:
        self.set_up_survey_tree_widget()
        self.set_up_obseration_plots_widget()
        self.set_up_obseration_results_plots_widget()
        self.set_up_drift_plot_widget()
        # self.observations_splitter.setSizes([1000, 10])

        # Initialize dialogs if necessary at the start of the application:
        self.dlg_corrections = DialogCorrections()
        self.dlg_autoselect_settings = DialogAutoselectSettings()
        self.dlg_estimation_settings = DialogEstimationSettings()

        # Overwrite/change setting from ui file, if necessary:
        self.dlg_estimation_settings.comboBox_adjustment_method.addItems(settings.ADJUSTMENT_METHODS.values())

        # Init models:
        self.observation_model = None
        self.setup_model = None
        self.results_station_model = None
        self.results_observation_model = None
        self.results_drift_model = None

        # Get system fonts:
        self.system_default_fixed_width_font = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont)

        # Set fonts:
        self.plainTextEdit_results_log.setFont(self.system_default_fixed_width_font)  # Monospace font

    def set_up_drift_plot_widget(self):
        """Set up `self.graphicsLayoutWidget_results_drift_plot`."""
        self.glw_drift_plot = self.graphicsLayoutWidget_results_drift_plot
        self.glw_drift_plot.setBackground('w')  # white background color
        # Create sub-plots:
        self.drift_plot = self.glw_drift_plot.addPlot(0, 0, name='drift_plot',
                                                             axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.drift_plot.setLabel(axis='left', text='')
        self.drift_plot.addLegend()

    def update_drift_plot(self):
        """Update the drift plot in the results tab.

        This method is used as slot. Hence, it will be invoked by signals from various GUI widgets that change the
        drift plot in the according plotting widget.
        """
        pass
        # Get GUI parameters:
        # - Selected LSM run:
        lsm_run_idx, lsm_run_time_str = self.get_selected_lsm_run()
        # - Selected station:
        idx_selected_station, selected_station_name = self.get_selected_station()
        if selected_station_name == 'All stations':
            selected_station_names = None
        else:
            selected_station_names = [selected_station_name]
        # - Selected Survey:
        idx_selected_survey, selected_survey_name = self.get_selected_survey()
        if selected_survey_name == 'All surveys':
            selected_survey_names = None
        else:
            selected_survey_names = [selected_survey_name]

        # Select lsm run:
        lsm_run = self.campaign.lsm_runs[lsm_run_idx]

        # - Enable/disable plot settings groupBox:
        if lsm_run.lsm_method == 'LSM_diff':
            self.groupBox_results_drift_plot.setEnabled(True)
        else:
            self.groupBox_results_drift_plot.setEnabled(False)

        if lsm_run.lsm_method == 'LSM_diff':
            offset_mugal = self.spinBox_results_drift_plot_v_offset.value()
            self.plot_drift_lsm_diff(lsm_run, surveys=selected_survey_names, stations=selected_station_names,
                                     offset=offset_mugal)
        elif lsm_run.lsm_method  == 'MLR_BEV':
            self.plot_drift_mlr_bev_legacy(lsm_run, surveys=selected_survey_names, stations=selected_station_names)
        else:
            self.drift_plot.clear()  # Clear drift plot

        # ---- invoke lsm-method specific plotting method ---

        # Get data from LSM object (selected LSM run):
        # - LSM method: LSMdiff/bev_mlr/...
        # - setup_df: pre-fit observations
        # - stat_obs_df: estimated station gravity
        # - drift_pol_df: drift polynomial coefficients

        # Merge observation and station data:

        # Filter data by surveys & stations (from GUI)

    def plot_drift_lsm_diff(self, lsm_run, surveys=None, stations=None, offset=0):
        """Create a drift plot for LSM runs based on differential observations (method: LSMdiff)

        Parameters:
        -----------
        surveys : `None` (default) or list of survey names (str)
            To filter for surveys that will be displayed.
        stations : `None` (default) or list of station names (str)
            To filter for stations that will be displayed.
        """
        NUM_ITEMS_IN_DRIFT_FUNCTION = 100

        self.drift_plot.clear()


        stat_obs_df = lsm_run.stat_obs_df
        drift_pol_df = lsm_run.drift_pol_df

        # Loop over surveys (setup data) in the selectd lsm run object and plot data:
        for survey_name, setup_df_orig in lsm_run.setups.items():
            # print(survey_name, setup)
            drift_pol_df_short = drift_pol_df.loc[drift_pol_df['survey_name'] == survey_name]

            # Prep data:
            setup_df = setup_df_orig.copy(deep=True)  # Make hard copy to protect original data!
            stat_obs_df_short = stat_obs_df.loc[:, ['station_name', 'g_est_mugal', 'sd_g_est_mugal']]
            setup_df = pd.merge(setup_df, stat_obs_df_short, on='station_name')
            setup_df['g_plot_mugal'] = setup_df['g_mugal'] - setup_df['g_est_mugal']
            setup_df.sort_values(by='delta_t_h', inplace=True)

            # Evaluate drift polynomial:
            coeff_list = drift_pol_df_short['coefficient'].to_list()
            coeff_list.reverse()
            coeff_list.append(0)
            delta_t_min_h = setup_df['delta_t_h'].min()  # = 0
            delta_t_max_h = setup_df['delta_t_h'].max()
            delta_t_h = np.linspace(delta_t_min_h, delta_t_max_h, NUM_ITEMS_IN_DRIFT_FUNCTION)
            drift_polynomial_mugal = np.polyval(coeff_list, delta_t_h)

            # Drift function time reference as UNIX time (needed for plots):
            epoch_unix_min = setup_df['epoch_unix'].min()
            epoch_unix_max = setup_df['epoch_unix'].max()
            delta_t_epoch_unix = np.linspace(epoch_unix_min, epoch_unix_max, NUM_ITEMS_IN_DRIFT_FUNCTION)

            # !!! Due to the differential observations, the constant bias (N0) of the gravity reading cannot be estimated!
            # In order to draw the drift polynomial function w.r.t. the gravity meter observations (for the sake of visual
            # assessment of the drift function), the const. bias N0 is approximated, see below.
            offset_mugal = setup_df['g_plot_mugal'].mean() - drift_polynomial_mugal.mean()
            yy_mugal = drift_polynomial_mugal + offset_mugal

            # Plot drift function:
            pen = pg.mkPen(color='b')
            self.drift_plot.plot(delta_t_epoch_unix, yy_mugal, name=f'drift: {survey_name}', pen=pen, symbol='o',
                                 symbolSize=10, symbolBrush=('b'))


        # Adjust plot window:
        self.drift_plot.showGrid(x=True, y=True)
        self.drift_plot.setLabel(axis='left', text='g [µGal]')
        self.drift_plot.setTitle(f'Drift function w.r.t. setup observations (arbitrary offset = {offset_mugal:0.1f} µGal)')
        self.drift_plot.autoRange()



        # EXAMPLE:
        # data = results_obs_df[column_name].values
        # obs_epoch_timestamps = (results_obs_df['ref_epoch_dt'].values - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1,'s')
        #
        # self.plot_xy_data(self.plot_obs_results, obs_epoch_timestamps, data, plot_name=column_name, color='b',
        #                   symbol='o', symbol_size=10)


    def plot_drift_mlr_bev_legacy(self, lsm_run, surveys=None, stations=None):
        """Create a drift plot for LSM runs using multiple linear regression (method: MLR BEV legacy)

        Parameters:
        -----------
        surveys : `None` (default) or list of survey names (str)
            To filter for surveys that will be displayed.
        stations : `None` (default) or list of station names (str)
            To filter for stations that will be displayed.
        """
        self.drift_plot.clear()
        pass

    def set_up_obseration_results_plots_widget(self):
        """Set up `self.graphicsLayoutWidget_results_observations_plots`."""
        self.glw_obs_results = self.graphicsLayoutWidget_results_observations_plots
        self.glw_obs_results.setBackground('w')  # white background color

        # Create sub-plots:
        self.plot_obs_results = self.glw_obs_results.addPlot(0, 0, name='obs_results',
                                                             axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_results.setLabel(axis='left', text='')
        self.plot_obs_results.addLegend()

    def plot_observation_results(self, results_obs_df=None, column_name=''):
        """Plots observation data to the GraphicsLayoutWidget.

        Notes
        -----
        If input parameters `` or/and `` is/are `None` or '', the plot content is deleted and the plot is resetted.
        """
        # Clear plot in any case:
        self.plot_obs_results.clear()

        if results_obs_df is not None:  # Data available for plotting
            # Get data:
            data = results_obs_df[column_name].values
            obs_epoch_timestamps = (results_obs_df['ref_epoch_dt'].values - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1,'s')

            self.plot_xy_data(self.plot_obs_results, obs_epoch_timestamps, data, plot_name=column_name, color='b',
                              symbol='o', symbol_size=10)
            self.plot_obs_results.showGrid(x=True, y=True)
            column_description = self.results_observation_model.get_plotable_columns()[column_name]
            self.plot_obs_results.setLabel(axis='left', text=column_description)
            self.plot_obs_results.autoRange()

    def update_results_obs_plots(self):
        """Update the observation results plots in the results tab."""
        # Get data from selected column
        col_idx, column_name = self.get_selected_obs_data_column()
        if col_idx != -1:
            filtered_results_obs_df = self.results_observation_model.get_model_data_df
        else:  # Invalid indices
            filtered_results_obs_df = None
        self.plot_observation_results(filtered_results_obs_df, column_name)

    def set_up_results_drift_view_model(self):
        """Set up the view model for the drift results table view."""
        try:
            self.results_drift_model = ResultsDriftModel(self.campaign.lsm_runs)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No LSM-adjustment results available!')
            self.statusBar().showMessage(f"No LSM-adjustment results available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.tableView_results_drift.setModel(self.results_drift_model)
            self.tableView_results_drift.resizeColumnsToContents()

    def update_results_drift_table_view(self, lsm_run_index: int, survey_name=None):
        """Update the drift results table view after changing the table model."""
        self.results_drift_model.update_view_model(lsm_run_index, survey_name=survey_name)
        self.results_drift_model.layoutChanged.emit()  # Show changes in table view
        self.tableView_results_drift.resizeColumnsToContents()

    def set_up_results_observations_view_model(self):
        """Set up the view model for the observation results table view."""
        try:
            self.results_observation_model = ResultsObservationModel(self.campaign.lsm_runs)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No LSM-adjustment results available!')
            self.statusBar().showMessage(f"No LSM-adjustment results available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.tableView_results_observations.setModel(self.results_observation_model)
            self.tableView_results_observations.resizeColumnsToContents()

    def update_results_observation_table_view(self, lsm_run_index: int, station_name=None, survey_name=None):
        """Update the observation results table view after changing the table model."""
        self.results_observation_model.update_view_model(lsm_run_index,
                                                         station_name=station_name,
                                                         survey_name=survey_name)
        self.results_observation_model.layoutChanged.emit()  # Show changes in table view
        self.tableView_results_observations.resizeColumnsToContents()
        self.update_comboBox_results_obs_plot_select_data_column_based_on_table_view()

    def set_up_results_stations_view_model(self):
        """Set up the view model for the station results table view."""
        try:
            self.results_station_model = ResultsStationModel(self.campaign.lsm_runs)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No LSM-adjustment results available!')
            self.statusBar().showMessage(f"No LSM-adjustment results available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.tableView_results_stations.setModel(self.results_station_model)
            self.tableView_results_stations.resizeColumnsToContents()

    def update_results_station_table_view(self, lsm_run_index: int, station_name=None, survey_name=None):
        """Update the station results table view after changing the table model."""
        self.results_station_model.update_view_model(lsm_run_index,
                                                     station_name=station_name,
                                                     survey_name=survey_name)
        self.results_station_model.layoutChanged.emit()  # Show changes in table view
        self.tableView_results_stations.resizeColumnsToContents()

    @pyqtSlot(int)
    def on_comboBox_results_lsm_run_selection_current_index_changed(self, index: int):
        """Invoked whenever the index of the selected item in the combobox changed."""
        # print(index)
        self.update_results_tab()

    @pyqtSlot(int)
    def on_comboBox_results_selection_station_current_index_changed(self, index: int):
        """Invoked whenever the index of the selected item in the combobox changed."""
        self.update_results_tab()

    @pyqtSlot(int)
    def on_comboBox_results_selection_survey_current_index_changed(self, index: int):
        """Invoked whenever the index of the selected item in the combobox changed."""
        self.update_results_tab()

    @pyqtSlot(int)
    def on_comboBox_results_obs_plot_select_data_column_current_index_changed(self, index: int):
        """Invoked whenever the index of the selected item in the combobox changed."""
        self.update_results_obs_plots()

    def update_results_tab(self, select_latest_item=False):
        """Update the results tab of the main tab widget with the content of the selected lsm run."""
        # Update combo box content for lsm run selection:
        self.update_comboBox_lsm_run_selection(select_latest_item=select_latest_item)

        # Get the currently selected lsm run object:
        idx, time_str = self.get_selected_lsm_run()
        if idx != -1:  # Valid index
            lsm_run = self.campaign.lsm_runs[idx]

            # Update the list in the combobox:
            self.update_comboBox_results_selection_station(observed_stations=lsm_run.observed_stations)
            self.update_comboBox_results_selection_surrvey(survey_names=list(lsm_run.setups.keys()))

            # Get data from LSM object an populate GUI widgets:
            # - Info tab:
            self.label_results_comment.setText(lsm_run.comment)
            self.label_results_adjustment_method.setText(settings.ADJUSTMENT_METHODS[lsm_run.lsm_method])
            self.label_results_time_and_date.setText(lsm_run.init_time.strftime("%Y-%m-%d, %H:%M:%S"))
            if lsm_run.s02_a_posteriori is not None:
                self.label_results_sig0.setText(f'{lsm_run.s02_a_posteriori:1.3f}')
            else:
                self.label_results_sig0.clear()
            if lsm_run.write_log:
                self.plainTextEdit_results_log.setPlainText(lsm_run.log_str)

            # TODO: Add further assignments for displaying data here!

            # Get station and/or survey names for filtering the displayed data:
            stat_idx, current_station_name = self.get_selected_station()
            if current_station_name == 'All stations':
                station_name = None
            else:
                station_name = current_station_name
            survey_idx, current_survey_name = self.get_selected_survey()
            if current_survey_name == 'All surveys':
                survey_name = None
            else:
                survey_name = current_survey_name

            # Update widgets:
            self.update_results_station_table_view(idx, station_name=station_name, survey_name=survey_name)
            self.update_results_observation_table_view(idx, station_name=station_name, survey_name=survey_name)
            self.update_results_drift_table_view(idx, survey_name=survey_name)
            self.update_results_obs_plots()
            self.update_drift_plot()
        else:  # invalid index => Reset results views
            self.label_results_comment.clear()
            self.label_results_adjustment_method.clear()
            self.label_results_time_and_date.clear()
            self.label_results_sig0.clear()
            self.plainTextEdit_results_log.clear()
            self.update_results_station_table_view(idx, station_name=None, survey_name=None)  # Can handle idx=-1
            self.update_results_observation_table_view(idx, station_name=None, survey_name=None)  # Can handle idx=-1
            self.update_results_drift_table_view(idx, survey_name=None)
            self.update_comboBox_results_selection_station(observed_stations=[])
            self.update_comboBox_results_selection_surrvey(survey_names=[])
            self.update_results_obs_plots()
            self.update_drift_plot()

    def update_comboBox_results_obs_plot_select_data_column_based_on_table_view(self):
        """Update the observaterion results data column selection combo box in the results tab."""
        self.comboBox_results_obs_plot_select_data_column.blockSignals(True)
        # Get data columns with data that is plottable from the observations results table view model:
        data_columns_dict = self.results_observation_model.get_plotable_columns()
        data_columns = list(data_columns_dict.keys())
        # Get current item:
        idx, current_column_name = self.get_selected_obs_data_column()
        self.comboBox_results_obs_plot_select_data_column.clear()
        self.comboBox_results_obs_plot_select_data_column.addItems(data_columns)
        # Try to select the previous item again:
        if idx != -1:  # Previous selection available
            try:
                self.comboBox_results_obs_plot_select_data_column.setCurrentText(current_column_name)
            except:
                pass
        self.comboBox_results_obs_plot_select_data_column.blockSignals(False)

    def update_comboBox_results_selection_station(self, observed_stations: list):
        """Update station selection combo box in the results tab, based on the observed station in the lsm run."""
        self.comboBox_results_selection_station.blockSignals(True)
        # Get current item:
        stat_idx, current_station_name = self.get_selected_station()
        self.comboBox_results_selection_station.clear()
        self.comboBox_results_selection_station.addItems(['All stations'] + observed_stations)
        # Try to select the previous item again:
        if stat_idx != -1:  # Previous selection available
            self.comboBox_results_selection_station.setCurrentText(current_station_name)
        self.comboBox_results_selection_station.blockSignals(False)

    def update_comboBox_results_selection_surrvey(self, survey_names: list):
        """Update the surve seletion combobox in the results tab, based on the current lsm run."""
        self.comboBox_results_selection_survey.blockSignals(True)
        # Get current item:
        survey_idx, current_survey_name = self.get_selected_survey()
        self.comboBox_results_selection_survey.clear()
        self.comboBox_results_selection_survey.addItems(['All surveys'] + survey_names)
        # Try to select the previous item again:
        if survey_idx != -1:  # Previous selection available
            self.comboBox_results_selection_survey.setCurrentText(current_survey_name)
        self.comboBox_results_selection_survey.blockSignals(False)

    def update_comboBox_lsm_run_selection(self, select_latest_item=False):
        """Update the LSM run selection combo box in the results tab, based on the available runs in the campaign."""
        self.comboBox_results_lsm_run_selection.blockSignals(True)
        # Get current item:
        idx, current_time_str = self.get_selected_lsm_run()
        self.comboBox_results_lsm_run_selection.clear()
        self.comboBox_results_lsm_run_selection.addItems(self.campaign.lsm_run_times)
        # Try to select the previous item again:
        if select_latest_item:
            idx = self.comboBox_results_lsm_run_selection.count() - 1
            self.comboBox_results_lsm_run_selection.setCurrentIndex(idx)
            # last index in list
        else:
            if idx != -1:  # Previous selection available
                try:
                    self.comboBox_results_lsm_run_selection.setCurrentText(current_time_str)
                except:
                    self.comboBox_results_lsm_run_selection.setCurrentIndex(0)  # 'All stations'
        self.comboBox_results_lsm_run_selection.blockSignals(False)

    def on_pushbutton_results_delete_lsm_run(self):
        """Delete the lsm object with index `idx` in the campaign."""
        idx, time_str = self.get_selected_lsm_run()
        if idx == -1:  # combo box is empty => No lsm run available!
            QMessageBox.warning(self, 'Warning!', 'No LSM adjustment run to be deleted!')
        else:
            try:
                self.campaign.delete_lsm_run(idx)
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f'LSM run "{time_str}" not deleted.')
            else:
                self.statusBar().showMessage(f'LSM run "{time_str}" deleted.')
                self.update_results_tab()

    def get_selected_lsm_run(self):
        """Get the selected lsm run in the results tab."""
        time_str = self.comboBox_results_lsm_run_selection.currentText()
        idx = self.comboBox_results_lsm_run_selection.currentIndex()
        return idx, time_str

    def get_selected_station(self):
        """Get the selected station in the results tab."""
        station_name = self.comboBox_results_selection_station.currentText()
        idx = self.comboBox_results_selection_station.currentIndex()
        return idx, station_name

    def get_selected_survey(self):
        """Get the selected survey in the results tab."""
        survey_name = self.comboBox_results_selection_survey.currentText()
        idx = self.comboBox_results_selection_survey.currentIndex()
        return idx, survey_name

    def get_selected_obs_data_column(self):
        """Get the selected observation results data column."""
        column_name = self.comboBox_results_obs_plot_select_data_column.currentText()
        idx = self.comboBox_results_obs_plot_select_data_column.currentIndex()
        return idx, column_name

    def on_checkBox_obs_plot_setup_data_state_changed(self):
        """Invoke, whenever the state of the checkbox changes."""
        self.refresh_observation_plot()

    def refresh_observation_plot(self):
        """Refresh the observation plot."""
        survey_name, setup_id = self.get_obs_tree_widget_selected_item()
        self.plot_observations(survey_name, setup_id)

    def on_pushbutton_obs_comp_setup_data(self):
        """Invoked when pushing the button 'pushbutton_obs_comp_setup_data'."""
        self.compute_setup_data_for_campaign()
        survey_name, setup_id = self.get_obs_tree_widget_selected_item()
        self.observation_model.update_view_model(survey_name, setup_id)
        self.refresh_observation_plot()
        self.update_setup_table_view(survey_name, setup_id)

    def update_setup_table_view(self, survey_name, setup_id):
        """Update the setups table view after changing the table model."""
        self.setup_model.update_view_model(survey_name, setup_id)
        self.setup_model.layoutChanged.emit()  # Show changes in table view
        self.tableView_observations_setups.resizeColumnsToContents()

    def compute_setup_data_for_campaign(self):
        """Compute setup data for the campaign."""
        try:
            self.campaign.calculate_setup_data(obs_type='reduced', set_epoch_of_first_obs_as_reference=True,
                                               active_obs_only_for_ref_epoch=settings.ACTIVE_OBS_ONLY_FOR_REF_EPOCH,
                                               verbose=IS_VERBOSE)
        except AssertionError as e:
            QMessageBox.critical(self, 'Error!', str(e))
            self.statusBar().showMessage(f"Error! No setup data computed")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
            self.statusBar().showMessage(f"Error! No setup data computed")
        else:
            # No errors when computing the setup data:
            self.statusBar().showMessage(f"Setup data computed successfully!")

    def on_pushbutton_obs_run_estimation(self):
        """Invoked when pushing the button 'on_pushbutton_obs_run_estimation'."""
        self.run_parameter_estimation()

    def run_parameter_estimation(self):
        """Run the parameter estimation process according to the defined estimation settings."""
        try:
            # Check if setup data is prepared:
            flag_setup_data_available = False
            for survey_name, survey in self.campaign.surveys.items():
                if survey.setup_df is not None:
                    if len(survey.setup_df) > 0:
                        flag_setup_data_available = True
            if not flag_setup_data_available:
                raise AssertionError('No setup data available!')

            # Get estimation settings
            lsm_method_value = self.dlg_estimation_settings.comboBox_adjustment_method.currentText()
            for key, value in settings.ADJUSTMENT_METHODS.items():
                if value == lsm_method_value:
                    lsm_method = key
            degree_drift_polynomial = self.dlg_estimation_settings.spinBox_degree_drift_polynomial.value()
            comment = self.dlg_estimation_settings.lineEdit_comment.text()
            sig0 = self.dlg_estimation_settings.doubleSpinBox_sig0.value()
            weight_factor_datum = self.dlg_estimation_settings.doubleSpinBox_weight_factor_datum.value()
            confidence_level_chi_test = self.dlg_estimation_settings.doubleSpinBox_conf_level_chi.value()
            confidence_level_tau_test = self.dlg_estimation_settings.doubleSpinBox_conf_level_tau.value()

            # Initialize LSM object and add it to the campaign object:
            self.campaign.initialize_and_add_lsm_run(lsm_method=lsm_method, comment=comment, write_log=True)

            # Run the estimation:
            if lsm_method == 'LSM_diff':
                self.campaign.lsm_runs[-1].adjust(drift_pol_degree=degree_drift_polynomial,
                                                  sig0_mugal=sig0,
                                                  scaling_factor_datum_observations=weight_factor_datum,
                                                  confidence_level_chi_test=confidence_level_chi_test,
                                                  confidence_level_tau_test=confidence_level_tau_test,
                                                  verbose=IS_VERBOSE)
                # TODO: Test drift plot!!!::
                self.campaign.lsm_runs[-1].create_drift_plot_matplotlib()
            elif lsm_method == 'MLR_BEV':
                self.campaign.lsm_runs[-1].adjust(drift_pol_degree=degree_drift_polynomial,
                                                  verbose=IS_VERBOSE)

        except AssertionError as e:
            QMessageBox.critical(self, 'Error!', str(e))
            self.statusBar().showMessage(f"Error! No parameters estimated.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
            self.statusBar().showMessage(f"Error! No parameters estimated.")
        else:
            # No errors when computing the setup data:
            self.statusBar().showMessage(f"Parameters estimated successfully!")
            # Update list of lsm runs in results tab of the GUI:
            self.update_results_tab(select_latest_item=True)
        pass

    def on_apply_autoselection(self):
        """Appply autoselection on the currently selected setup or survey according to the predefined setttings."""

        # Get autoselect parameters from settings dialog
        flag_apply_tilt = self.dlg_autoselect_settings.checkBox_tilt.isChecked()
        flag_apply_g_sd = self.dlg_autoselect_settings.checkBox_sd.isChecked()
        flag_apply_delta_g = self.dlg_autoselect_settings.checkBox_delta_g.isChecked()

        treshold_g_sd_mugal = int(self.dlg_autoselect_settings.spinBox_sd.text())
        treshold_tilt_arcsec = int(self.dlg_autoselect_settings.spinBox_tilt.text())
        treshold_delta_sd_mugal = int(self.dlg_autoselect_settings.spinBox_delta_g.text())
        delta_g_number_of_points = int(self.dlg_autoselect_settings.spinBox_n.text())

        if self.dlg_autoselect_settings.radioButton_ref_data_reduced_observations.isChecked():
            reference_data = 'reduced'
        else:
            reference_data = 'observed'

        # Get selected setup or survey
        survey_name, setup_id = self.get_obs_tree_widget_selected_item()

        # Apply autoselection on selected survey/setup (campaign data):
        surv = self.campaign.surveys[survey_name]

        # Check if reduced observations are available, if they are required as reference:
        if reference_data == 'reduced':
            if surv.obs_df['g_red_mugal'].isnull().any() or surv.obs_df['sd_g_red_mugal'].isnull().any():
                if IS_VERBOSE:
                    print('Reduced observations not available!')
                QMessageBox.critical(self, 'Error!',
                                     'Reduced observations (reference for autoselection) are not avialable yet!')
                return

        if flag_apply_tilt:
            surv.autselect_tilt(threshold_arcsec=treshold_tilt_arcsec, setup_id=setup_id)
        if flag_apply_g_sd:
            surv.autselect_g_sd(threshold_mugal=treshold_g_sd_mugal, obs_type=reference_data, setup_id=setup_id,
                                verbose=IS_VERBOSE)
        if flag_apply_delta_g:
            surv.autselect_delta_g(threshold_mugal=treshold_delta_sd_mugal, n_obs=delta_g_number_of_points,
                                   obs_type=reference_data, setup_id=setup_id, verbose=IS_VERBOSE)

        # Update data visualization in GUI
        self.update_obs_table_view(survey_name, setup_id)
        self.plot_observations(survey_name, setup_id)
        self.update_obs_tree_widgget_from_observation_model()
        pass

    def update_obs_tree_widgget_from_observation_model(self):
        """Update observation tree widget by checking data in the observation model."""
        self.treeWidget_observations.blockSignals(True)  # Block any signals when changing the checked state
        survey_name = self.observation_model.data_survey_name
        setup_ids = self.observation_model.get_data.loc[:, 'setup_id'].unique()
        for setup_id in setup_ids:
            keep_obs_flags_of_setup = self.observation_model.get_data.loc[
                self.observation_model.get_data['setup_id'] == setup_id, 'keep_obs']
            if all(keep_obs_flags_of_setup):
                tree_item_check_state = Qt.Checked
            elif any(keep_obs_flags_of_setup):
                tree_item_check_state = Qt.PartiallyChecked
            else:  # all are = False
                tree_item_check_state = Qt.Unchecked
            # Set Checked state:
            for tree_item_idx in range(self.treeWidget_observations.topLevelItemCount()):
                if self.treeWidget_observations.topLevelItem(tree_item_idx).text(0) == survey_name:
                    item = self.treeWidget_observations.topLevelItem(tree_item_idx)
                    for child_item_idx in range(item.childCount()):
                        if int(item.child(child_item_idx).text(0)) == setup_id:
                            setup_item = item.child(child_item_idx)
                            setup_item.setCheckState(0, tree_item_check_state)  # finally set checked state!
                            break
        self.treeWidget_observations.blockSignals(False)
        pass

    def get_obs_tree_widget_selected_item(self):
        """Returns survey name and setup nam of the selected items in the observation tree widget."""
        items = self.treeWidget_observations.selectedItems()
        if len(items) == 1:  # Only one item in tree view selected
            item = items[0]
            if item.parent() is None:  # Is a survey
                survey_name = item.text(0)  # Column 0 = Survey name
                setup_id = None  # No setup selected
            else:
                parent = item.parent()
                survey_name = parent.text(0)  # Column 0 = Survey name
                setup_id = int(item.text(0))
            self.update_obs_table_view(survey_name, setup_id)
            self.plot_observations(survey_name, setup_id)
        else:
            if IS_VERBOSE:
                print('No item or multiple items selected!')
            survey_name = None
            setup_id = None
        return survey_name, setup_id

    def on_menu_observations_autoselection_settings(self):
        """Launch dialog for defining the autoselection settings."""
        return_value = self.dlg_autoselect_settings.exec()
        pass

    def on_menu_estimation_settings(self):
        """Launch dialog for defining the estimation settings."""
        return_value = self.dlg_estimation_settings.exec()
        pass

    def set_up_obseration_plots_widget(self):
        """Set up `self.GraphicsLayoutWidget_observations`."""
        l = self.GraphicsLayoutWidget_observations
        l.setBackground('w')  # white background color

        date_axis = TimeAxisItem(orientation='bottom')

        # Create sub-plots:
        # Gravity g [µGal]
        self.plot_obs_g = l.addPlot(0, 0, name='plot_obs_g', axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_g.setLabel(axis='left', text='g [µGal]')
        self.plot_obs_g.addLegend()

        # Standard deviation of gravity g [µGal]
        self.plot_obs_sd_g = l.addPlot(1, 0, name='plot_obs_sd_g',
                                       axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_sd_g.setLabel(axis='left', text='sd_g [µGal]')
        self.plot_obs_sd_g.addLegend()
        self.plot_obs_sd_g.setXLink(self.plot_obs_g)

        # Instrument tilt in X and Y directions [arcsec]
        self.plot_obs_tilt = l.addPlot(2, 0, name='plot_obs_tilt',
                                       axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_tilt.setLabel(axis='left', text='tilt [arcsec]')
        self.plot_obs_tilt.addLegend()
        self.plot_obs_tilt.setXLink(self.plot_obs_g)

        # Observation corrections [µGal]
        self.plot_obs_corrections = l.addPlot(3, 0, name='plot_obs_corrections',
                                              axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_corrections.setLabel(axis='left', text='Corrections [µGal]')
        self.plot_obs_corrections.addLegend()
        self.plot_obs_corrections.setXLink(self.plot_obs_g)

    def plot_observations(self, survey_name, setup_id):
        """Plots observation data to the GraphicsLayoutWidget."""
        obs_df = self.observation_model.get_data
        setup_df = self.observation_model.get_setup_data
        obs_epoch_timestamps = (obs_df['obs_epoch'].values - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1,
                                                                                                                     's')
        # Distinguish between survey and setup data here, in necessary!
        # if setup_id is None:  # survey selected:
        #     pass
        # else:  # setup selected:
        #     pass

        # Plot reduced or unreduced observations:
        flag_show_reduced_observations = False
        if self.checkBox_obs_plot_reduced_observations.checkState() == Qt.Checked and not any(
                obs_df['g_red_mugal'].isnull()):
            # Reduced observations are available and will be shown:
            flag_show_reduced_observations = True
        elif self.checkBox_obs_plot_reduced_observations.checkState() == Qt.Checked and any(
                obs_df['g_red_mugal'].isnull()):
            flag_show_reduced_observations = False
            QMessageBox.warning(self, 'Warning!', 'Reduced observations are not available!')
            self.checkBox_obs_plot_reduced_observations.setChecked(Qt.Unchecked)

            # Get data:
        if flag_show_reduced_observations:
            g_mugal = obs_df['g_red_mugal'].values
            sd_g_mugal = obs_df['sd_g_red_mugal'].values
            corr_tide = obs_df['corr_tide_red_mugal'].values
            corr_tide_name = self.campaign.surveys[survey_name].red_tide_correction_type
            ref_height_name = self.campaign.surveys[survey_name].red_reference_height_type
        else:
            g_mugal = obs_df['g_obs_mugal'].values
            sd_g_mugal = obs_df['sd_g_obs_mugal'].values
            corr_tide = obs_df['corr_tide_mugal'].values
            corr_tide_name = self.campaign.surveys[survey_name].obs_tide_correction_type
            ref_height_name = self.campaign.surveys[survey_name].obs_reference_height_type

        # Gravity g [µGal]
        # - Plot with marker symbols according to their 'keep_obs' states and connect the 'sigPointsClicked' event.
        self.plot_obs_g.clear()
        pen = pg.mkPen(color='b')
        flags_keep_obs = obs_df['keep_obs'].values
        symbol_brushes = []
        for flag in flags_keep_obs:
            if flag:
                symbol_brushes.append(self.BRUSH_ACTIVE_OBS)
            else:
                symbol_brushes.append(self.BRUSH_INACTIVE_OBS)

        # setup data: g
        if setup_df is not None and self.checkBox_obs_plot_setup_data.isChecked():
            self.plot_xy_data(self.plot_obs_g, setup_df['epoch_unix'].values, setup_df['g_mugal'].values,
                              plot_name='setup', color='k', symbol='x', symbol_size=25)

        # Type of 'self.plot_obs_g_data_item': PlotDataItem
        self.plot_obs_g_data_item = self.plot_obs_g.plot(obs_epoch_timestamps, g_mugal, name=f'Ref.: {ref_height_name}',
                                                         pen=pen, symbol='o', symbolSize=10, symbolBrush=symbol_brushes)
        self.plot_obs_g_data_item.sigPointsClicked.connect(self.on_observation_plot_data_item_clicked)
        self.plot_obs_g.showGrid(x=True, y=True)
        self.plot_obs_g.autoRange()

        # Standard deviation of gravity g [µGal]
        self.plot_obs_sd_g.clear()
        # setup data: sd_g
        if setup_df is not None and self.checkBox_obs_plot_setup_data.isChecked():
            self.plot_xy_data(self.plot_obs_sd_g, setup_df['epoch_unix'].values, setup_df['sd_g_mugal'].values,
                              plot_name='setup', color='k', symbol='x', symbol_size=25)
        self.plot_xy_data(self.plot_obs_sd_g, obs_epoch_timestamps, sd_g_mugal, plot_name='sd_g_mugal',
                          color='b', symbol='o', symbol_size=10)

        self.plot_obs_sd_g.showGrid(x=True, y=True)
        self.plot_obs_sd_g.autoRange()

        # Instrument tilt in X and Y directions [arcsec]
        self.plot_obs_tilt.clear()
        tilt_x = obs_df['tiltx'].values
        tilt_y = obs_df['tilty'].values
        self.plot_xy_data(self.plot_obs_tilt, obs_epoch_timestamps, tilt_x, plot_name='X', color='b', symbol='o',
                          symbol_size=10)
        self.plot_xy_data(self.plot_obs_tilt, obs_epoch_timestamps, tilt_y, plot_name='Y', color='r', symbol='t',
                          symbol_size=10)
        self.plot_obs_tilt.showGrid(x=True, y=True)
        self.plot_obs_tilt.autoRange()

        # Observation corrections [µGal]
        self.plot_obs_corrections.clear()
        self.plot_xy_data(self.plot_obs_corrections, obs_epoch_timestamps, corr_tide,
                          plot_name=f'tides ({corr_tide_name})', color='b', symbol='o', symbol_size=10)
        self.plot_obs_corrections.showGrid(x=True, y=True)
        self.plot_obs_corrections.autoRange()

        self.plot_obs_g.autoRange()  # Finally adjust data range to g values!

    def plot_xy_data(self, plot_item, x, y, plot_name, color='k', symbol='o', symbol_size=10):
        """Plot XY-data."""
        pen = pg.mkPen(color=color)
        plot_item.plot(x, y, name=plot_name, pen=pen, symbol=symbol, symbolSize=symbol_size, symbolBrush=(color))

    def set_keep_obs_markers_in_obs_plot(self, index: int, keep_obs_flag: bool):
        """Set a marker (scatter symbols) in the observation plot according to the `keep_obs` flag.

        Notes
        -----
        This method only changes ONE single marker symbol! Hence, `index` is a scalar.

        Parameters
        ----------
        index : int
            This is the index of the marker symbol to be changed in the current plot (and data in
            `self.observation_model.get_data` respectively).
        keep_obs_flag : bool
            Indicates whether the observation should ba active or inactive. The symbol brush is set accordingly.
        """

        if keep_obs_flag:
            self.plot_obs_g_data_item.scatter.points()[index].setBrush(self.BRUSH_ACTIVE_OBS)
        else:
            self.plot_obs_g_data_item.scatter.points()[index].setBrush(self.BRUSH_INACTIVE_OBS)

    def on_observation_plot_data_item_clicked(self, points, ev):
        """Invoked whenever a data point in the observation time series plot is clicked.

        Notes
        -----
        Whenever an observation data point is clicked the `keep_obs` flag toggles (True/False) and the marker (brush)
        is set accordingly to visualize the situation in the observation timeseries plot. Additionally this flag is
        changed in the observation model data (dataframe accessed via `self.observation_model.get_data`) and the
        `dataChanged` signal is emitted. The signal triggers the method `self.on_observation_model_data_changed`
        where the observation tree view is changed accordingly and the `keep_ob` fag is set in the survey in the
        campaign data (method: `self.campaign.surveys[<survey_nanem>].activate_observation()`).
        """
        # print(points, ev)
        # print(f'Number of clicked points: {len(ev)}')  # Initially ALL points under the mouse cursor are selected!
        # print(f'-----------------------------')

        # Get first selected point (only select ONE point!):
        spot_item = ev[0]  # <class 'pyqtgraph.graphicsItems.ScatterPlotItem.SpotItem'>
        # Select item in observation_model and toggel the "keep_obs" state:
        row = self.observation_model.get_data.index[spot_item._index]
        if self.observation_model.get_data.at[row, 'keep_obs']:
            self.observation_model.get_data.at[row, 'keep_obs'] = False
            # Handled in `self.on_observation_model_data_changed` with `self.set_keep_obs_markers_in_obs_plot` instead:
            # spot_item.setBrush(self.BRUSH_INACTIVE_OBS)
        else:
            self.observation_model.get_data.at[row, 'keep_obs'] = True
            # Handled in `self.on_observation_model_data_changed` with `self.set_keep_obs_markers_in_obs_plot` instead:
            # spot_item.setBrush(self.BRUSH_ACTIVE_OBS)
        index = self.tableView_observations.model().index(spot_item._index, Survey.get_obs_df_column_index('keep_obs'))
        self.observation_model.dataChanged.emit(index, index)  # Triggers all following events...

    def on_observation_model_data_changed(self, topLeft, bottomRight, role):
        """Invoked whenever data in the observation table view changed."""
        if IS_VERBOSE:
            # print('TableModel: dataChanged: ', topLeft, bottomRight, role)
            pass
        if topLeft == bottomRight:  # Only one item selected?
            index = topLeft
            # Get column and row indices for dataframe:
            row = self.observation_model.get_data.index[index.row()]
            col = self.observation_model.get_data.columns[index.column()]
            if col == 'keep_obs':
                flag_keep_obs = self.observation_model.get_data.at[row, col]
                # Change data in `survey.obs_df`:
                survey_name = self.observation_model.data_survey_name
                survey = self.campaign.surveys[survey_name]
                survey.activate_observation(row, flag_keep_obs)

                # Change check state of obs tree widget according to "keep_obs" flags of setup
                setup_id = self.observation_model.get_data.at[row, 'setup_id']
                keep_obs_flags_of_setup = self.observation_model.get_data.loc[
                    self.observation_model.get_data['setup_id'] == setup_id, 'keep_obs']
                if all(keep_obs_flags_of_setup):
                    tree_item_check_state = Qt.Checked
                elif any(keep_obs_flags_of_setup):
                    tree_item_check_state = Qt.PartiallyChecked
                else:  # all are = False
                    tree_item_check_state = Qt.Unchecked
                # Set Checked state:
                self.treeWidget_observations.blockSignals(True)  # Block any signals when changing the checked state
                for tree_item_idx in range(self.treeWidget_observations.topLevelItemCount()):
                    if self.treeWidget_observations.topLevelItem(tree_item_idx).text(0) == survey_name:
                        item = self.treeWidget_observations.topLevelItem(tree_item_idx)
                        for child_item_idx in range(item.childCount()):
                            if int(item.child(child_item_idx).text(0)) == setup_id:
                                setup_item = item.child(child_item_idx)
                                setup_item.setCheckState(0, tree_item_check_state)  # finally set checked state!
                                break
                self.treeWidget_observations.blockSignals(False)

                # Plot 'keep_obs' marker symbols in the observation timeseries plot:
                self.set_keep_obs_markers_in_obs_plot(index.row(), flag_keep_obs)
        else:
            pass  # More than one items selected/changed.

    def on_tree_widget_item_changed(self, item, column):
        """Invoked whenever an item in th observation tree view is changed, e.g. if the check-state changes."""
        self.treeWidget_observations.blockSignals(True)  # To avoid recursive effects
        if IS_VERBOSE:
            print('TreeView: itemChanged: ', item, column)
            print(f' - CheckState of item "{item.text(0)}": {item.checkState(0)}')
        flag_checked_state = checked_state_to_bool(item.checkState(0))
        # Is parent (survey) or child (setup):
        if item.parent() is None:  # Is survey item
            survey_name = item.text(0)
            setup_id = None
            if flag_checked_state == True:
                self.campaign.activate_survey(survey_name, verbose=False)
            else:
                self.campaign.deactivate_survey(survey_name, verbose=False)
            self.on_obs_tree_widget_item_selected()
            if IS_VERBOSE:
                print('-----------Parent changed---------------------')
        else:  # Is a setup item
            # Update table view and data in dataframe:
            survey_name = item.parent().text(0)
            setup_id = int(item.text(0))
            self.campaign.surveys[survey_name].activate_setup(setup_id, flag_checked_state)
        self.treeWidget_observations.blockSignals(False)

    @pyqtSlot()
    def on_obs_tree_widget_item_selected(self):
        """Invoked whenever an item in the observation tree is selected."""
        items = self.treeWidget_observations.selectedItems()
        if len(items) == 1:  # Only one item in tree view selected
            item = items[0]
            if item.parent() is None:  # Is a survey
                survey_name = item.text(0)  # Column 0 = Survey name
                setup_id = None  # No setup selected
            else:
                parent = item.parent()
                survey_name = parent.text(0)  # Column 0 = Survey name
                setup_id = int(item.text(0))
            self.update_obs_table_view(survey_name, setup_id)
            self.update_setup_table_view(survey_name, setup_id)
            self.plot_observations(survey_name, setup_id)
        else:
            if IS_VERBOSE:
                print('No item or multiple items selected!')

    def update_obs_table_view(self, survey_name: str, setup_id: int):
        """Update the observation table view according to the selected survey and instrument setup."""
        if IS_VERBOSE:
            print(f'survey name: {survey_name}; setup ID: {setup_id}')
        # Update the observation table view model according to the selected
        self.observation_model.update_view_model(survey_name, setup_id)  # Show added survey in table
        self.observation_model.layoutChanged.emit()  # Show changes in table view
        self.tableView_observations.resizeColumnsToContents()

    def set_up_survey_tree_widget(self):
        """Set up the survey tree widget."""
        self.treeWidget_observations.setColumnCount(3)
        self.treeWidget_observations.setHeaderLabels(['Survey', 'Station', '#Obs'])
        header = self.treeWidget_observations.header()
        header.setVisible(True)
        header.setSectionResizeMode(QHeaderView.ResizeToContents)
        header.setStretchLastSection(True)
        # The following line raises the following error:
        # "Process finished with exit code 139 (interrupted by signal 11: SIGSEGV)"
        # header.setSectionResizeMode(5, QHeaderView.Stretch)

    def set_up_setup_view_model(self):
        """Set up the view model for the setup data table view."""
        try:
            self.setup_model = SetupTableModel(self.campaign.surveys)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No surveys available!')
            self.statusBar().showMessage(f"No surveys available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.tableView_observations_setups.setModel(self.setup_model)
            self.tableView_observations_setups.resizeColumnsToContents()

    @pyqtSlot()
    def set_up_observation_view_model(self):
        """Set up observation data view model and show observation data table view."""
        # Set model:
        try:
            self.observation_model = ObservationTableModel(self.campaign.surveys)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No surveys available!')
            self.statusBar().showMessage(f"No surveys available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.tableView_observations.setModel(self.observation_model)
            # self.connect_station_model_to_table_view()  # Set view
            self.tableView_observations.resizeColumnsToContents()
            self.statusBar().showMessage(f"{self.campaign.number_of_surveys} surveys in current campaign.")

            self.observation_model.dataChanged.connect(self.on_observation_model_data_changed)

    def populate_survey_tree_widget(self):
        """Populate the survey tree widget."""
        # Delete existing items:
        self.treeWidget_observations.blockSignals(True)
        self.delete_all_items_from_survey_tree_widget()

        # Add new items:
        # Parent items (surveys):
        for survey_name, survey in self.campaign.surveys.items():
            obs_df = survey.obs_df.sort_values('obs_epoch').copy(deep=True)
            num_of_obs_in_survey = survey.get_number_of_observations()

            parent = QTreeWidgetItem(self.treeWidget_observations)
            parent.setText(0, survey_name)
            parent.setText(2, str(num_of_obs_in_survey))
            parent.setFlags(parent.flags() | Qt.ItemIsTristate | Qt.ItemIsUserCheckable)
            if survey.keep_survey:
                parent.setCheckState(0, Qt.Checked)
            else:
                parent.setCheckState(0, Qt.Unchecked)

            # Loop over instrument setups in survey:
            setup_ids = obs_df['setup_id'].unique()
            for setup_id in setup_ids:
                # get all observations
                setup_keep_obs_flags = obs_df.loc[obs_df['setup_id'] == setup_id, 'keep_obs']
                setup_station_names = obs_df.loc[obs_df['setup_id'] == setup_id, 'station_name']
                num_of_obs_in_setup = len(setup_keep_obs_flags)

                # get station name (has to be unique within one setup!):
                if len(setup_station_names.unique()) > 1:
                    QMessageBox.critical(self, 'Error!', 'Setup {} includes observations to more than one '
                                                         'station ({})!'.format(setup_id,
                                                                                ', '.join(
                                                                                    setup_station_names.unique().tolist())))
                    self.delete_all_items_from_survey_tree_widget()
                    self.treeWidget_observations.blockSignals(False)
                    return
                else:
                    setup_station_name = setup_station_names.unique()[0]
                    child = QTreeWidgetItem(parent)
                    child.setFlags(child.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsTristate)
                    child.setText(0, str(setup_id))
                    child.setText(2, str(num_of_obs_in_setup))
                    child.setText(1, setup_station_name)

                    if setup_keep_obs_flags.all():
                        child.setCheckState(0, Qt.Checked)
                    elif setup_keep_obs_flags.any():
                        child.setCheckState(0, Qt.PartiallyChecked)
                    else:
                        child.setCheckState(0, Qt.Unhecked)
            parent.setExpanded(True)  # Expand the current parent
        self.treeWidget_observations.show()
        self.treeWidget_observations.blockSignals(False)
        # self.treeWidget_observations.expandToDepth(0)  # Expand items to a certain depth

    def delete_all_items_from_survey_tree_widget(self):
        """Delete all items from the survey tree widget."""
        self.treeWidget_observations.clear()

    @pyqtSlot(str)
    def on_lineEdit_filter_stat_name_textChanged(self, text):
        """Only display the observed stations."""
        search = QRegExp(text, Qt.CaseInsensitive, QRegExp.RegExp)
        try:
            self.proxy_station_model.setFilterKeyColumn(self.campaign.stations._STAT_DF_COLUMNS.index('station_name'))
            self.proxy_station_model.setFilterRegExp(search)
        except:
            if IS_VERBOSE:
                print('No filter proxy model connected.')

    @pyqtSlot(int)  # Required, because 2 signals are emitted and one (int) has to be selected!
    def on_checkBox_filter_observed_stat_only_toggled(self, state):
        """Event handler for the filter observed stations only checkbox."""
        self.proxy_station_model.setFilterKeyColumn(self.campaign.stations._STAT_DF_COLUMNS.index('is_observed'))
        if state == Qt.Checked:
            self.proxy_station_model.setFilterFixedString('True')
        else:
            self.proxy_station_model.setFilterFixedString('')

    @pyqtSlot()
    def exit_application(self):
        """Exit the application."""
        sys.exit()

    def on_menu_observations_corrections(self):
        """Launch diaglog to select and apply observation corrections."""
        # dlg = DialogCorrections()
        return_value = self.dlg_corrections.exec()
        if return_value == QDialog.Accepted:
            flag_corrections_ok, error_msg = self.apply_observation_corrections()
            if flag_corrections_ok:
                # Load survey from campaing data to observations vie model:
                self.observation_model.load_surveys(self.campaign.surveys)
                self.on_obs_tree_widget_item_selected()
                self.statusBar().showMessage(f"Observation corrections applied.")
            else:
                QMessageBox.critical('Error!', error_msg)
                self.statusBar().showMessage(f"Error: No observation corrections applied.")
        else:
            self.statusBar().showMessage(f"No observation corrections applied.")

    def apply_observation_corrections(self):
        """Apply observation corrections according to the selected settings."""
        flag_selection_ok = True
        error_msg = ''
        if self.dlg_corrections.radioButton_corr_ref_heights_ground.isChecked():
            target_ref_height = 'ground'
        elif self.dlg_corrections.radioButton_corr_ref_heights_control_point.isChecked():
            target_ref_height = 'control_point'
        elif self.dlg_corrections.radioButton_corr_ref_heights_sensor.isChecked():
            target_ref_height = 'sensor_height'
        elif self.dlg_corrections.radioButton_corr_ref_heights_instrument_top.isChecked():
            target_ref_height = 'instrument_top'
        else:
            flag_selection_ok = False
            error_msg = f'Invalid selection of reference height in GUI (observation corrections dialog).'
            # QMessageBox.critical('Error!', error_msg)

        if self.dlg_corrections.radioButton_corr_tides_no_correction.isChecked():
            target_tide_corr = 'no_tide_corr'
        elif self.dlg_corrections.radioButton_corr_tides_cg5_model.isChecked():
            target_tide_corr = 'cg5_longman1959'
        else:
            flag_selection_ok = False
            error_msg = f'Invalid selection of tidal correction in GUI (observation corrections dialog).'
            # QMessageBox.critical('Error!', error_msg)

        if flag_selection_ok:
            flag_corrections_ok, error_msg = self.campaign.reduce_observations_in_all_surveys(
                target_ref_height=target_ref_height,
                target_tide_corr=target_tide_corr,
                verbose=IS_VERBOSE)
        else:
            flag_corrections_ok = False

        return flag_corrections_ok, error_msg

    @pyqtSlot()
    def on_menu_file_new_campaign(self):
        """Launching dialog to create a new campaign."""
        dlg = DialogNewCampaign(old_campaign=self.campaign)
        return_value = dlg.exec()
        if return_value == QDialog.Accepted:
            # Create Campaign object:
            self.campaign = Campaign(
                campaign_name=dlg.lineEdit_campaign_name.text(),
                output_directory=dlg.lineEdit_output_directory.text(),
                surveys=None,  # Always use non-mutable default arguments!
                stations=None,  # Always use non-mutable default arguments!
            )
            # Enable/disable main menu items:
            self.menuAdd_Survey.setEnabled(True)
            self.action_Add_Stations.setEnabled(True)

            self.statusBar().showMessage(f"New Campaign created (name: {self.campaign.campaign_name}, "
                                         f"output directory: {self.campaign.output_directory})")
            self.setWindowTitle('GravTools - Campaign: ' + self.campaign.campaign_name)
            self.enable_menu_observations_based_on_campaign_data()  # Disable when starting a new campaign (no surveys).

            # Set up view models and views for this campaign:
            self.set_up_station_view_model()
            self.set_up_observation_view_model()
            self.set_up_setup_view_model()
            self.set_up_results_stations_view_model()
            self.set_up_results_observations_view_model()
            self.set_up_results_drift_view_model()

        elif return_value == QDialog.Rejected:
            self.statusBar().showMessage(f"Canceled creating new campaign.")

    def enable_menu_observations_based_on_campaign_data(self):
        """Enable the main menu item `observations` is the campaign contains at least one survey."""
        if self.campaign.number_of_surveys > 0:
            self.menu_Observations.setEnabled(True)
            self.groupBox_obs_view_options.setEnabled(True)
            self.groupBox_obs_data_manipulation.setEnabled(True)
        else:
            self.menu_Observations.setEnabled(False)
            self.groupBox_obs_view_options.setEnabled(False)
            self.groupBox_obs_data_manipulation.setEnabled(False)

    @pyqtSlot()
    def on_menu_file_load_stations(self):
        """Launch dialog to load stations to project."""
        dlg = DialogLoadStations()
        return_value = dlg.exec()
        if return_value == QDialog.Accepted:
            # Load Stations to campaign
            number_of_stations_old = self.campaign.stations.get_number_of_stations
            try:
                # WARNING: The following line is required in order to prevent a memory error with
                # "QSortFilterProxyModelPrivate::proxy_to_source()"
                # => When adding stations to the model do the follwong steps:
                # 1.) connect the station model ("self.station_model")
                # 2.) Add station data
                # 3.) Set up an connect the proxy model for sorting and filtering ("self.set_up_proxy_station_model")
                self.connect_station_model_to_table_view()

                # WARNING: Whenever new station data is loaded with "self.campaign.add_stations_from_oesgn_table_file",
                # the reference between the "view model" () and the "campaign stations dataframe"
                # (self.campaign.stations.stat_df) is broken!
                # => Therefore, these two have to be reassigned in order to mirror the same data! This is done by
                #    calling the method "self.station_model.load_stat_df(self.campaign.stations.stat_df)".
                self.campaign.add_stations_from_oesgn_table_file(dlg.lineEdit_oesgn_table_file_path.text(),
                                                                 verbose=IS_VERBOSE)
                self.campaign.synchronize_stations_and_surveys(verbose=IS_VERBOSE)
                self.refresh_stations_table_model_and_view()
                self.set_up_proxy_station_model()

            except FileNotFoundError:
                QMessageBox.critical(self, 'File not found error', f'"{dlg.lineEdit_oesgn_table_file_path.text()}" '
                                                                   f'not found.')
                self.statusBar().showMessage(f"No stations added.")
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f"No stations added.")
            else:
                self.enable_station_view_options_based_on_model()
                number_of_stations_added = self.campaign.stations.get_number_of_stations - number_of_stations_old
                self.statusBar().showMessage(f"{number_of_stations_added} stations added.")
        else:
            self.statusBar().showMessage(f"No stations added.")

    def enable_station_view_options_based_on_model(self):
        """Enable or disable the station view options based on the number of stations in the model."""
        if len(self.station_model._data) > 0:
            self.groupBox_filter_options.setEnabled(True)
            self.groupBox_edit_options.setEnabled(True)
        else:
            self.groupBox_filter_options.setEnabled(False)
            self.groupBox_edit_options.setEnabled(False)

    def refresh_stations_table_model_and_view(self):
        """Refrech the station table model and the respective table view."""
        self.station_model.load_stat_df(self.campaign.stations.stat_df)
        self.station_model.layoutChanged.emit()  # Refresh the table view, after changing the model's size.
        self.tableView_Stations.resizeColumnsToContents()

    @pyqtSlot()
    def set_up_station_view_model(self):
        """Set up station data view model and show station data table view."""
        # Set model:
        try:
            self.station_model = StationTableModel(self.campaign.stations.stat_df)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No stations available!')
            self.statusBar().showMessage(f"No stations available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.connect_station_model_to_table_view()  # Set view
            self.tableView_Stations.resizeColumnsToContents()
            self.statusBar().showMessage(f"{self.campaign.number_of_stations} stations in current campaign.")

    def connect_station_model_to_table_view(self):
        """Connect model from table view for stations."""
        self.tableView_Stations.setModel(self.station_model)

    def set_up_proxy_station_model(self):
        """Set up and activate sort and filter proxy model for the station view model."""
        self.proxy_station_model = QSortFilterProxyModel()
        # self.proxy_station_model.setDynamicSortFilter(True)
        # self.proxy_station_model.setDynamicSortFilter(False)
        self.proxy_station_model.setSourceModel(self.station_model)

        # Enable sorting the table:
        self.tableView_Stations.setSortingEnabled(True)

        # Set filter:
        self.proxy_station_model.setFilterKeyColumn(0)  # Select column for the filtering  # TODO: set dynamically!
        # self.proxy_station_model.setFilterFixedString('True')
        self.connect_proxy_station_model_to_table_view()  # Set view

    def connect_proxy_station_model_to_table_view(self):
        """Connect proxy model from table view for stations."""
        self.tableView_Stations.setModel(self.proxy_station_model)

    @pyqtSlot()
    def on_menu_file_load_survey_from_cg5_observation_file(self):
        """Launch file selection dialog to select a CG5 observation file and load the data to the campaign."""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        cg5_obs_file_filename, _ = QFileDialog.getOpenFileName(self,
                                                               'Select CG5 observation file',
                                                               DEFAULT_CG5_OBS_FILE_PATH,
                                                               "CG5 observation file (*.TXT)",
                                                               options=options)
        if cg5_obs_file_filename:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            cg5_obs_file_filename = QDir.toNativeSeparators(cg5_obs_file_filename)
            # Add survey data to Campaign:
            try:
                new_cg5_survey = Survey.from_cg5_obs_file(cg5_obs_file_filename)
                flag_no_duplicate = self.campaign.add_survey(survey_add=new_cg5_survey, verbose=IS_VERBOSE)
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f"No survey data added.")
            else:
                if flag_no_duplicate:
                    # No problems occurred:
                    # Check if all expected information is available in the survey and raise warning:
                    if new_cg5_survey.obs_tide_correction_type == 'unknown':
                        QMessageBox.warning(self, 'Warning!',
                                            'Type of tidal correction is unknown! '
                                            'Check, if the "CG-5 OPTIONS" block in the input file is missing. '
                                            ' => Reduced observations not calculated yet!')
                    else:
                        # Calculate reduced observation data:
                        if settings.CALCULATE_REDUCED_OBS_WHEN_LOADING_DATA:
                            flag_corrections_ok, error_msg = self.apply_observation_corrections()
                            if flag_corrections_ok:
                                # Load survey from campaing data to observations vie model:
                                # self.observation_model.load_surveys(self.campaign.surveys)
                                # self.on_obs_tree_widget_item_selected()
                                self.statusBar().showMessage(f"Observation corrections applied.")
                            else:
                                QMessageBox.critical('Error!', error_msg)
                                self.statusBar().showMessage(f"Error: No observation corrections applied.")
                    self.connect_station_model_to_table_view()  # Disconnect sort & filter proxy model from station view.
                    self.campaign.synchronize_stations_and_surveys(verbose=IS_VERBOSE)
                    self.refresh_stations_table_model_and_view()
                    self.populate_survey_tree_widget()
                    # Select the added survey in the tree view:
                    for tree_item_idx in range(self.treeWidget_observations.topLevelItemCount()):
                        if self.treeWidget_observations.topLevelItem(tree_item_idx).text(0) == new_cg5_survey.name:
                            self.treeWidget_observations.topLevelItem(tree_item_idx).setSelected(True)
                    self.enable_menu_observations_based_on_campaign_data()
                    self.statusBar().showMessage(f"Survey {new_cg5_survey.name} "
                                                 f"({new_cg5_survey.get_number_of_observations()} observations) added.")
                else:
                    QMessageBox.warning(self,
                                        'Warning!',
                                        f'the current campaign already contains a survey named {new_cg5_survey.name}. '
                                        f'Survey names have to be unique within a campaign.')
                    self.statusBar().showMessage(f"No survey data added.")
            finally:
                self.set_up_proxy_station_model()  # Re-connect the sort & filter proxy model to the station view.
                self.enable_station_view_options_based_on_model()
        else:
            self.statusBar().showMessage(f"No survey data added.")

    def closeEvent(self, event):
        """Ask the user whether to close the window or not!"""
        reply = QMessageBox.question(self, 'Message', "Are you sure to quit?",
                                     QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


class DialogNewCampaign(QDialog, Ui_Dialog_new_Campaign):
    """New campaign dialog."""

    def __init__(self, parent=None, old_campaign=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)
        # connect signals and slots:
        self.pushButton_change_output_directory.clicked.connect(self.get_output_directory_dialog)
        # self.accepted.connect(self.test_fun)

        # Initialize the name and the output directory from an existing campaign (is it exists):
        # self.old_campaign = old_campaign
        if isinstance(old_campaign, Campaign):  # not None!
            self.lineEdit_output_directory.setText(old_campaign.output_directory)
            self.lineEdit_campaign_name.setText(old_campaign.campaign_name)
        else:
            self.lineEdit_output_directory.setText(DEFAULT_OUTPUT_DIR)
            self.lineEdit_campaign_name.setText('')

    def get_output_directory_dialog(self):
        """Open dialog to get the output directory."""

        initial_folder_path = self.lineEdit_output_directory.text()
        output_dir_name = QFileDialog.getExistingDirectory(self, 'Select a directory', initial_folder_path)

        if output_dir_name:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            output_dir_name = QDir.toNativeSeparators(output_dir_name)
            if IS_VERBOSE:
                print(output_dir_name)

        # Check, if path exists:
        if os.path.isdir(output_dir_name):
            self.lineEdit_output_directory.setText(output_dir_name)

    def done(self, a0: int) -> None:
        """Entered, whenever the dialog is about to be closed.
        Examples: https://programtalk.com/python-examples/PyQt5.Qt.QDialog.done/"""
        # Check form input:
        if a0 == 1:  # If the "OK" button was pressed
            if len(self.lineEdit_campaign_name.text()) == 0:
                self.lineEdit_campaign_name.setFocus()
                msg_critical = QMessageBox.critical(self, "ERROR", "Please enter a Campaign name!")
                return
            if not os.path.isdir(self.lineEdit_output_directory.text()):
                self.lineEdit_output_directory.setFocus()
                msg_critical = QMessageBox.critical(self, "ERROR", "The output directory does not exist!")
                return
        return QDialog.done(self, a0)  # Everything is OK!


class DialogLoadStations(QDialog, Ui_Dialog_load_stations):
    """Dialog to load stations to project."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)
        # connect signals and slots:
        self.pushButton_select_oesgn_table.clicked.connect(self.select_oesgn_table_file_dialog)

    def select_oesgn_table_file_dialog(self):
        """Open file selection dialog."""
        # initial_oesgn_file_path = self.lineEdit_oesgn_table_file_path.text()
        if self.lineEdit_oesgn_table_file_path.text():
            initial_oesgn_file_path = os.path.dirname(os.path.abspath(self.lineEdit_oesgn_table_file_path.text()))
        else:
            initial_oesgn_file_path = self.lineEdit_oesgn_table_file_path.text()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        oesgn_filename, _ = QFileDialog.getOpenFileName(self, 'Select OESGN table file', initial_oesgn_file_path,
                                                        "OESGN table (*.TAB)", options=options)
        if oesgn_filename:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            oesgn_filename = QDir.toNativeSeparators(oesgn_filename)
            self.lineEdit_oesgn_table_file_path.setText(oesgn_filename)


class DialogCorrections(QDialog, Ui_Dialog_corrections):
    """Dialog to select and apply observation corrections."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)
        # connect signals and slots:
        pass


class DialogAutoselectSettings(QDialog, Ui_Dialog_autoselection_settings):
    """Dialog to define the autoselect settings."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)
        # connect signals and slots:
        pass


class DialogEstimationSettings(QDialog, Ui_Dialog_estimation_settings):
    """Dialog to define the estimation settings."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)




if __name__ == "__main__":
    """Main Program."""
    # Create the application
    app = QApplication(sys.argv)

    # Create and show the application's main window
    main_window = MainWindow()
    main_window.show()
    # Run the application's main loop:
    sys.exit(
        app.exec())  # exit or error code of Qt (app.exec_) is passed to sys.exit. Terminates pgm with standard python method
