"""Graphical user interface of GravTools written with PyQt.

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

Notes
-----
The graphical layout of the GUI was created by using the Qt Designer (<https://www.qt.io/>).

"""

import sys
import os
import warnings
from PyQt5.QtWidgets import QApplication, QDialog, QMainWindow, QFileDialog, QMessageBox, QTreeWidgetItem, \
    QHeaderView, QInputDialog, QMenu
from PyQt5.QtCore import QDir, Qt, QSortFilterProxyModel, pyqtSlot, QRegExp, QPoint
from PyQt5 import QtGui
import datetime as dt
import pyqtgraph as pg
import pyqtgraph.exporters
import numpy as np
import pandas as pd
import pytz

# optional imports:
try:
    import geopandas
except ImportError:
    _has_geopandas = False
    warnings.warn('The optional dependency "geopandas" is not installed. Some features, e.g. GIS file export,'
                  ' will not be available.', UserWarning)
else:
    _has_geopandas = True

from gravtools.gui.MainWindow import Ui_MainWindow
from gravtools.gui.dialog_new_campaign import Ui_Dialog_new_Campaign
from gravtools.gui.dialog_load_stations import Ui_Dialog_load_stations
from gravtools.gui.dialog_autoselection_settings import Ui_Dialog_autoselection_settings
from gravtools.gui.dialog_estimation_settings import Ui_Dialog_estimation_settings
from gravtools.gui.dialog_setup_data import Ui_Dialog_setup_data
from gravtools.gui.dialog_export_results import Ui_Dialog_export_results
from gravtools.gui.dialog_options import Ui_Dialog_options
from gravtools.gui.dialog_about import Ui_Dialog_about
from gravtools.gui.dialog_gis_export_settings import Ui_Dialog_gis_settings
from gravtools.gui.gui_model_station_table import StationTableModel
from gravtools.gui.gui_model_setup_table import SetupTableModel
from gravtools.gui.gui_model_observation_table import ObservationTableModel
from gravtools.gui.gui_model_results_stat_table import ResultsStationModel
from gravtools.gui.gui_model_results_obs_table import ResultsObservationModel
from gravtools.gui.gui_model_results_drift_table import ResultsDriftModel
from gravtools.gui.gui_model_results_correlation_matrix_table import ResultsCorrelationMatrixModel
from gravtools.gui.gui_model_results_vg_table import ResultsVGModel
from gravtools.gui.gui_model_survey_table import SurveyTableModel
from gravtools.gui.gui_misc import get_station_color_dict, checked_state_to_bool, resize_table_view_columns
from gravtools.gui.dlg_correction_time_series import DialogCorrectionTimeSeries
from gravtools.gui.dlg_corrections import DialogCorrections
from gravtools.gui.survey_table_handler import SurveyTableHandler
from gravtools.gui.cumstom_widgets import ScrollMessageBox
from gravtools import __version__, __author__, __git_repo__, __email__, __copyright__, __pypi_repo__

from gravtools.models.survey import Survey
from gravtools.models.campaign import Campaign
from gravtools.models.misc import time_it, conditional_decorator, unique_ordered_list
from gravtools import settings

DEFAULT_OUTPUT_DIR = os.path.abspath(os.getcwd())  # Current working directory
DEFAULT_WORK_DIR_PATH = os.path.abspath(os.getcwd())  # Current working directory
IS_VERBOSE = True  # Define, whether screen output is enabled.

MARKER_SYMBOL_ORDER = ('o', 't', 'x', 's', 'star', '+', 'd', 't1', 'p', 't2', 'h', 't3')
MARKER_COLOR_ODER = ('b', 'r', 'g', 'c', 'm', 'y')


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
        self.action_Corrections.triggered.connect(self.on_menu_observations_corrections)
        self.action_Flag_observations.triggered.connect(self.on_manu_observations_flag_observations)
        self.action_Autoselection_settings.triggered.connect(self.on_menu_observations_autoselection_settings)
        self.action_Estimation_settings.triggered.connect(self.on_menu_estimation_settings)
        self.action_Gis_Export_settings.triggered.connect(self.on_menu_gis_export_settings)
        self.action_Setup_data_options.triggered.connect(self.on_menu_observations_setup_data)
        self.action_Export_Results.triggered.connect(self.on_menu_file_export_results)
        self.action_Options.triggered.connect(self.on_menu_file_options)
        self.action_About.triggered.connect(self.on_menu_help_about)
        self.pushButton_obs_apply_autoselect_current_data.pressed.connect(self.on_apply_autoselection)
        self.pushButton_obs_comp_setup_data.pressed.connect(self.on_pushbutton_obs_comp_setup_data)
        self.pushButton_obs_run_estimation.pressed.connect(self.on_pushbutton_obs_run_estimation)
        self.pushButton_results_delete_lsm_run.pressed.connect(self.on_pushbutton_results_delete_lsm_run)
        self.pushButton_results_export_shapefile.pressed.connect(self.on_pushButton_results_export_shapefile)
        self.action_from_CG5_observation_file.triggered.connect(self.on_menu_file_load_survey_from_cg5_observation_file)
        self.action_from_oesgn_table.triggered.connect(self.on_menu_file_load_stations_from_oesgn_table)
        self.action_from_csv_file.triggered.connect(self.on_menu_file_load_stations_from_csv_file)
        self.action_Correction_time_series.triggered.connect(self.action_correction_time_series_triggered)
        self.lineEdit_filter_stat_name.textChanged.connect(self.on_lineEdit_filter_stat_name_textChanged)
        self.checkBox_filter_observed_stat_only.stateChanged.connect(self.on_checkBox_filter_observed_stat_only_toggled)
        self.checkBox_obs_plot_setup_data.stateChanged.connect(self.on_checkBox_obs_plot_setup_data_state_changed)
        self.treeWidget_observations.itemSelectionChanged.connect(self.on_obs_tree_widget_item_selected)
        self.treeWidget_observations.itemChanged.connect(self.on_tree_widget_item_changed)
        self.checkBox_obs_plot_reduced_observations.clicked.connect(self.on_obs_tree_widget_item_selected)
        self.pushButton_obs_collaps_all.pressed.connect(self.survey_tree_widget_collapse_all)
        self.pushButton_obs_expand_all.pressed.connect(self.survey_tree_widget_expand_all)
        self.comboBox_results_lsm_run_selection.currentIndexChanged.connect(
            self.on_comboBox_results_lsm_run_selection_current_index_changed)
        self.comboBox_results_selection_station.currentIndexChanged.connect(
            self.on_comboBox_results_selection_station_current_index_changed)
        self.comboBox_results_selection_survey.currentIndexChanged.connect(
            self.on_comboBox_results_selection_survey_current_index_changed)
        self.comboBox_results_obs_plot_select_data_column.currentIndexChanged.connect(
            self.on_comboBox_results_obs_plot_select_data_column_current_index_changed)
        self.spinBox_results_drift_plot_v_offset.valueChanged.connect(
            self.on_spinBox_results_drift_plot_v_offset_value_changed)
        self.checkBox_stations_map_show_stat_name_labels.stateChanged.connect(
            self.on_checkBox_stations_map_show_stat_name_labels_state_changed)
        self.radioButton_results_vg_plot_details.toggled.connect(self.on_results_vg_plot_type_radiobuttons_changed)
        self.radioButton_results_vg_plot_full_polynomial.toggled.connect(self.on_results_vg_plot_type_radiobuttons_changed)
        self.checkBox_stations_map_show_stat_name_labels.stateChanged.connect(
            self.on_checkBox_stations_map_show_stat_name_labels_state_changed)
        self.checkBox_results_vg_plot_show_residuals.stateChanged.connect(
            self.checkBox_results_vg_plot_show_residuals_state_changed)
        self.radioButton_results_obs_plot_timeseries.toggled.connect(self.update_results_obs_plots)
        self.radioButton_results_obs_plot_histogram.toggled.connect(self.update_results_obs_plots)
        self.comboBox_results_obs_plot_hist_method.currentIndexChanged.connect(self.update_results_obs_plots)
        self.comboBox_results_obs_plot_hist_method.currentIndexChanged.connect(
            self.on_histogram_bin_method_currentIndexChanged)
        self.spinBox_results_obs_plot_number_bins.valueChanged.connect(self.update_results_obs_plots)
        self.comboBox_results_stations_statistics_select_col.currentIndexChanged.connect(
            self.display_station_results_statistics)
        self.tableView_surveys.customContextMenuRequested.connect(self.tableView_surveys_right_click)
        # self.action_Load_Campaign.triggered.connect(self.on_action_Load_Campaign_triggered)  # Not needed!?!
        # self.action_Change_output_directory.triggered.connect(self.on_action_Change_output_directory_triggered)  # Not needed!?!

        # Set up GUI items and widgets:
        self.set_up_survey_tree_widget()
        self.set_up_obseration_plots_widget()
        self.set_up_obseration_results_plots_widget()
        self.set_up_obseration_results_plots_hist_method_comboBox()
        self.set_up_drift_plot_widget()
        self.set_up_stations_map()
        self.set_up_vg_plot_widget()
        # self.observations_splitter.setSizes([1000, 10])

        # Initialize dialogs if necessary at the start of the application:
        self.dlg_corrections = DialogCorrections(self)
        self.dlg_autoselect_settings = DialogAutoselectSettings(self)
        self.dlg_estimation_settings = DialogEstimationSettings(self)
        self.dlg_gis_export_settings = DialogGisExportSettings(self)
        self.dlg_options = DialogOptions(self)
        self.dlg_setup_data = DialogSetupData(self)
        self.dlg_about = DialogAbout(self)
        self.dlg_correction_time_series = DialogCorrectionTimeSeries(self)
        # self.dlg_about.label_author.setText(__author__)
        self.dlg_about.label_version.setText(__version__)
        self.dlg_about.label_git_repo.setText(__git_repo__)
        self.dlg_about.label_pypi_repo.setText(__pypi_repo__)
        self.dlg_about.label_email.setText(__email__)
        self.dlg_about.label_copyright.setText(__copyright__)

        # Estimation settings GUI:
        self.dlg_estimation_settings.comboBox_adjustment_method.currentIndexChanged.connect(
            self.on_dlg_estimation_settings_comboBox_adjustment_method_current_index_changed)
        self.dlg_estimation_settings.comboBox_iteration_approach.currentIndexChanged.connect(
            self.on_dlg_estimation_settings_comboBox_iteration_approach_current_index_changed)

        # Overwrite/change setting from ui file, if necessary:
        self.dlg_estimation_settings.comboBox_adjustment_method.addItems(settings.ADJUSTMENT_METHODS.values())
        self.dlg_estimation_settings.comboBox_iteration_approach.addItems(settings.ITERATION_APPROACHES.keys())

        # Configure GUI according to optional dependencies:
        if _has_geopandas:
            self.action_Gis_Export_settings.setEnabled(True)
            self.groupBox_gis_data.setEnabled(True)
            self.dlg_gis_export_settings.lineEdit_stat_coord_epsg.setText(f'{settings.DEFAULT_EPSG_CODE:d}')
            self.dlg_gis_export_settings.lineEdit_filename_obs_results_shp.setText(f'{settings.DEFUALT_FILENAME_OBERVATION_RESULTS_SHP}')
            self.dlg_gis_export_settings.lineEdit_filename_stat_results_shp.setText(f'{settings.DEFUALT_FILENAME_STATION_RESULTS_SHP}')
        else:
            self.groupBox_gis_data.setEnabled(False)
            self.action_Gis_Export_settings.setEnabled(False)

        # Init models:
        self.station_model = None
        self.observation_model = None
        self.setup_model = None
        self.results_station_model = None
        self.results_observation_model = None
        self.results_drift_model = None
        self.survey_model = None

        # Get system fonts:
        self.system_default_fixed_width_font = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont)

        # Set fonts:
        self.plainTextEdit_results_log.setFont(self.system_default_fixed_width_font)  # Monospace font

        # Inits misc:
        self.station_colors_dict_results = {}  # set in self.update_results_tab()

    @pyqtSlot(QPoint)
    def tableView_surveys_right_click(self, position):
        """Invoked on right mouse click on a field in the stations results table view."""
        idx = self.tableView_surveys.indexAt(position)  # Returns QModelIndex
        if not idx.isValid():
            return
        survey_name = self.survey_model.get_survey_name_by_row_index(idx.row())
        ctx_menu = QMenu()
        ctx_menu.setToolTip('Remove the selected survey from the campaign. Existing adjustment results containing '
                            'observations of this survey will remain present.')
        delete_survey = ctx_menu.addAction(f'Delete Survey {survey_name}')
        action = ctx_menu.exec_(QtGui.QCursor.pos())
        if action == delete_survey:
            reply = QMessageBox.question(self, 'Message', f'Remove survey "{survey_name}" from the campaign?',
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.No:
                return
            self.campaign.remove_survey(survey_name=survey_name, verbose=IS_VERBOSE)

            # Stations tab:
            self.update_comboBox_stations_selection_surrvey(survey_names=self.campaign.survey_names)
            self.campaign.synchronize_stations_and_surveys(verbose=IS_VERBOSE)
            self.refresh_stations_table_model_and_view()
            # self.set_up_proxy_station_model()  # Re-connect the sort & filter proxy model to the station view.  # Required?
            self.on_checkBox_filter_observed_stat_only_toggled(
                state=self.checkBox_filter_observed_stat_only.checkState())
            self.enable_station_view_options_based_on_model()
            # Show observed stations only based on Checkbox state:
            self.on_checkBox_filter_observed_stat_only_toggled(self.checkBox_filter_observed_stat_only.checkState(),
                                                               auto_range_stations_plot=True)

            # Observations tab:
            selected_survey_name, selected_setup_id = self.get_obs_tree_widget_selected_item()

            self.populate_survey_tree_widget()
            # Select the previously selected survey again, if still available.
            # - Select nothing, if no survey left => Clear Plots and tables in observations tab
            # - Select the first survey if the previously selected one was deleted
            if self.campaign.number_of_surveys == 0:
                self.update_obs_table_view(survey_name=None, setup_id=None)
                self.update_setup_table_view(survey_name=None, setup_id=None)
                self.plot_observations(survey_name=None)
                # self.on_obs_tree_widget_item_selected()
            elif selected_survey_name == survey_name:  # Select first survey in the tree widget
                self.treeWidget_observations.topLevelItem(0).setSelected(True)
            else:  # Restore the previous selection
                for tree_item_idx in range(self.treeWidget_observations.topLevelItemCount()):
                    if self.treeWidget_observations.topLevelItem(tree_item_idx).text(0) == selected_survey_name:
                        if selected_setup_id is None:
                            self.treeWidget_observations.topLevelItem(tree_item_idx).setSelected(True)
                        else:
                            for child_item_idx in range(self.treeWidget_observations.topLevelItem(tree_item_idx).childCount()):
                                if self.treeWidget_observations.topLevelItem(tree_item_idx).child(child_item_idx).text(0) == selected_setup_id:
                                    self.treeWidget_observations.topLevelItem(tree_item_idx).child(child_item_idx).setSelected()
            self.enable_menu_observations_based_on_campaign_data()
            self.set_up_survey_view_model()

            # TODO: Handle results!
            # Problem: Wenn man einen Survey entfernt, so kommt es zu Fehlermeldungen, wenn man bestehende Ergebnisse
            # dieses Surveys im Results tab anzeigen möchte! Offenbar wird hierbei auf campaing.surveys zugegriffen.
            # Das muss geändert werden, so dass sämtliche Daten, die zur Visualisierung der Ergebnisse nötig sind,
            # im LSM object gespeichert werden und dort zur Verfügung stehen.

            self.statusBar().showMessage(
                f'Survey "{survey_name}" removed from campaign.')

    def action_correction_time_series_triggered(self):
        """Launch dialog for managing correction time series data."""
        _ = self.dlg_correction_time_series.exec()

    @pyqtSlot()
    def on_pushButton_results_export_shapefile(self):
        """Invoked whenever pressing the button."""
        # Get the currently selected lsm run object:
        idx, time_str = self.get_selected_lsm_run()
        if idx == -1:  # invalid index
            self.statusBar().showMessage(f'No data selected for export to shapefiles...')
            return
        lsm_run = self.campaign.lsm_runs[idx]
        try:
            epsg_code = int(self.dlg_gis_export_settings.lineEdit_stat_coord_epsg.text())
        except ValueError:
            QMessageBox.critical(self, 'Error!', 'Invalid EPSG code. Need to be an integer value.')
            return

        # Get output directory:
        if self.dlg_gis_export_settings.radioButton_campaign_output_dir.isChecked():
            if self.dlg_gis_export_settings.lineEdit_output_subdir.text():
                gis_output_dir = os.path.join(self.campaign.output_directory,
                                              self.dlg_gis_export_settings.lineEdit_output_subdir.text())
                if not os.path.isdir(gis_output_dir):
                    try:
                        os.mkdir(gis_output_dir)
                    except PermissionError:
                        QMessageBox.critical(self, 'Error!',
                                             f'Cannot create the output directory for GIS results: {gis_output_dir}')
                        return
            else:
                gis_output_dir = self.campaign.output_directory
        else:
            gis_output_dir = self.dlg_gis_export_settings.lineEdit_gis_output_dir.text()
        if not os.path.isdir(gis_output_dir):
            QMessageBox.critical(self, 'Error!',
                                 f'Invalid output directory for GIS files: {gis_output_dir}')
        else:
            # Export station results:
            if self.dlg_gis_export_settings.checkBox_export_stat_results_shp.checkState() == Qt.Checked:
                if self.dlg_gis_export_settings.checkBox_add_lsm_method_filename.checkState() == Qt.Checked:
                    filename = self.dlg_gis_export_settings.lineEdit_filename_stat_results_shp.text() + lsm_run.lsm_method + '.shp'
                else:
                    filename = self.dlg_gis_export_settings.lineEdit_filename_stat_results_shp.text() + '.shp'
                filename = os.path.join(gis_output_dir, filename)
                try:
                    lsm_run.export_stat_results_shapefile(filename=filename, epsg_code=epsg_code)
                except AttributeError:
                    QMessageBox.warning(self, f'Export not available!', f'Export of station results to a shapefile is not '
                                                                       f'supported by the lsm method {lsm_run.lsm_method}.')
                except Exception as e:
                    QMessageBox.critical(self, 'Error!', str(e))
                else:
                    self.statusBar().showMessage(f'Save station results of lsm run "{lsm_run.comment}" to: {filename}')

            # Export observations results:
            if self.dlg_gis_export_settings.checkBox_export_obs_results_shp.checkState() == Qt.Checked:
                if self.dlg_gis_export_settings.checkBox_add_lsm_method_filename.checkState() == Qt.Checked:
                    filename = self.dlg_gis_export_settings.lineEdit_filename_obs_results_shp.text() + lsm_run.lsm_method + '.shp'
                else:
                    filename = self.dlg_gis_export_settings.lineEdit_filename_obs_results_shp.text() + '.shp'
                filename = os.path.join(gis_output_dir, filename)
                try:
                    lsm_run.export_obs_results_shapefile(filename=filename, epsg_code=epsg_code)
                except AttributeError:
                    QMessageBox.warning(self, f'Export not available!', f'Export of observation results to a shapefile is not '
                                                                        f'supported by the lsm method {lsm_run.lsm_method}.')
                except Exception as e:
                    QMessageBox.critical(self, 'Error!', str(e))
                else:
                    self.statusBar().showMessage(f'Save observations results of lsm run "{lsm_run.comment}" to: {filename}')


    @pyqtSlot()
    def on_action_Change_Campaign_name_triggered(self):
        """Invoked whenever the menu item change campaign name is pressed."""
        # Get name:
        new_campaign_name, flag_done = QInputDialog.getText(self, 'Enter the new campaign name',
                                                            'New campaign name (No blanks or special characters allowed):')
        if flag_done:
            try:
                self.campaign.change_campaign_name(new_campaign_name)
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f'No new campaign name set.')
            else:  # Valid input
                # Set new window title:
                self.setWindowTitle('GravTools - Campaign: ' + self.campaign.campaign_name)
                self.statusBar().showMessage(f'New campaign name set to: {self.campaign.campaign_name}.')
        else:
            self.statusBar().showMessage(f'No new campaign name set.')

    @pyqtSlot()
    def on_action_Change_output_directory_triggered(self):
        """Invoked whenever the menu item change output directory is pressed."""
        self.change_campaign_output_directory()

    def change_campaign_output_directory(self):
        """Change the output directory of the current campaign."""
        initial_folder_path = self.campaign.output_directory
        output_dir_name = QFileDialog.getExistingDirectory(self, 'Select a directory', initial_folder_path,
                                                           QFileDialog.ShowDirsOnly)
        if output_dir_name:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            output_dir_name = QDir.toNativeSeparators(output_dir_name)
            # Check, if path exists:
            if os.path.isdir(output_dir_name):
                try:
                    self.campaign.set_output_directory(output_dir_name)
                except Exception as e:
                    QMessageBox.critical(self, 'Error!', str(e))
                    self.statusBar().showMessage(f"No valid output directory selected.")
                else:
                    if IS_VERBOSE:
                        print(f'New output directory: {output_dir_name}')
                    self.statusBar().showMessage(f'New output directory: {output_dir_name}')
            else:
                self.statusBar().showMessage(f'Output directory "{output_dir_name}" does not exist!')
                QMessageBox.critical(self, 'Error!', f'Directory "{output_dir_name}" does not exist!')
        else:
            self.statusBar().showMessage(f'No output directory selected.')

    @pyqtSlot(int)
    def on_dlg_estimation_settings_comboBox_iteration_approach_current_index_changed(self, index: int):
        """Invoked whenever the iteration approach in the estimation settings dialog changed."""
        selected_iteration_approach = self.dlg_estimation_settings.comboBox_iteration_approach.currentText()
        if selected_iteration_approach == 'Multiplicative':
            self.dlg_estimation_settings.groupBox_iterative_scaling_additive.setEnabled(False)
            self.dlg_estimation_settings.groupBox_iterative_scaling_multiplicative.setEnabled(True)
        elif selected_iteration_approach == 'Additive':
            self.dlg_estimation_settings.groupBox_iterative_scaling_additive.setEnabled(True)
            self.dlg_estimation_settings.groupBox_iterative_scaling_multiplicative.setEnabled(False)
        else:
            self.dlg_estimation_settings.groupBox_iterative_scaling_additive.setEnabled(False)
            self.dlg_estimation_settings.groupBox_iterative_scaling_multiplicative.setEnabled(False)
            QMessageBox.warning(self, 'Warning!', 'Unknown iteration approach selected!')
            self.statusBar().showMessage(f"Unknown iteration approach selected!")

    @pyqtSlot(int)
    def on_dlg_estimation_settings_comboBox_adjustment_method_current_index_changed(self, index: int):
        """Invoked whenever the adjustment method changed in the estimation settings dialog."""
        # enable/disable GUI elements in the estimation settings dialog according to the selected method:
        selected_method = self.dlg_estimation_settings.comboBox_adjustment_method.currentText()
        if selected_method == 'MLR (BEV legacy processing)':
            self.dlg_estimation_settings.groupBox_constraints.setEnabled(False)
            self.dlg_estimation_settings.groupBox_statistical_tests.setEnabled(False)
            self.dlg_estimation_settings.groupBox_observations.setEnabled(False)
            self.dlg_estimation_settings.doubleSpinBox_sig0.setEnabled(False)
            self.dlg_estimation_settings.label_sig0.setEnabled(False)
            self.dlg_estimation_settings.groupBox_iterative_scaling.setEnabled(False)
            self.dlg_estimation_settings.groupBox_drift_polynomial_advanced.setEnabled(False)
            self.dlg_estimation_settings.groupBox_vg_polynomial.setEnabled(False)
            self.dlg_estimation_settings.checkBox_iterative_s0_scaling.setEnabled(False)
            self.dlg_estimation_settings.groupBox_se_determination.setEnabled(False)
        elif selected_method == 'LSM (differential observations)':
            self.dlg_estimation_settings.groupBox_constraints.setEnabled(True)
            self.dlg_estimation_settings.groupBox_statistical_tests.setEnabled(True)
            self.dlg_estimation_settings.groupBox_observations.setEnabled(True)
            self.dlg_estimation_settings.doubleSpinBox_sig0.setEnabled(True)
            self.dlg_estimation_settings.label_sig0.setEnabled(True)
            self.dlg_estimation_settings.groupBox_iterative_scaling.setEnabled(True)
            self.dlg_estimation_settings.groupBox_drift_polynomial_advanced.setEnabled(True)
            self.dlg_estimation_settings.groupBox_vg_polynomial.setEnabled(False)
            self.dlg_estimation_settings.checkBox_iterative_s0_scaling.setEnabled(True)
            self.dlg_estimation_settings.groupBox_se_determination.setEnabled(True)
        elif selected_method == 'LSM (non-differential observations)':
            self.dlg_estimation_settings.groupBox_constraints.setEnabled(True)
            self.dlg_estimation_settings.groupBox_statistical_tests.setEnabled(True)
            self.dlg_estimation_settings.groupBox_observations.setEnabled(True)
            self.dlg_estimation_settings.doubleSpinBox_sig0.setEnabled(True)
            self.dlg_estimation_settings.label_sig0.setEnabled(True)
            self.dlg_estimation_settings.groupBox_iterative_scaling.setEnabled(True)
            self.dlg_estimation_settings.groupBox_drift_polynomial_advanced.setEnabled(True)
            self.dlg_estimation_settings.groupBox_vg_polynomial.setEnabled(False)
            self.dlg_estimation_settings.checkBox_iterative_s0_scaling.setEnabled(True)
            self.dlg_estimation_settings.groupBox_se_determination.setEnabled(True)
        elif selected_method == 'VG LSM (non-differential observations)':
            self.dlg_estimation_settings.groupBox_constraints.setEnabled(False)
            self.dlg_estimation_settings.groupBox_statistical_tests.setEnabled(True)
            self.dlg_estimation_settings.groupBox_observations.setEnabled(True)
            self.dlg_estimation_settings.doubleSpinBox_sig0.setEnabled(True)
            self.dlg_estimation_settings.label_sig0.setEnabled(True)
            self.dlg_estimation_settings.groupBox_iterative_scaling.setEnabled(False)
            self.dlg_estimation_settings.groupBox_drift_polynomial_advanced.setEnabled(True)
            self.dlg_estimation_settings.groupBox_vg_polynomial.setEnabled(True)
            self.dlg_estimation_settings.checkBox_iterative_s0_scaling.setEnabled(False)
            self.dlg_estimation_settings.groupBox_se_determination.setEnabled(False)
        else:
            # Enable all and show warning:
            self.dlg_estimation_settings.groupBox_constraints.setEnabled(True)
            self.dlg_estimation_settings.groupBox_statistical_tests.setEnabled(True)
            self.dlg_estimation_settings.groupBox_observations.setEnabled(True)
            self.dlg_estimation_settings.doubleSpinBox_sig0.setEnabled(True)
            self.dlg_estimation_settings.label_sig0.setEnabled(True)
            self.dlg_estimation_settings.groupBox_iterative_scaling.setEnabled(True)
            self.dlg_estimation_settings.groupBox_drift_polynomial_advanced.setEnabled(True)
            self.dlg_estimation_settings.groupBox_vg_polynomial.setEnabled(True)
            self.dlg_estimation_settings.checkBox_iterative_s0_scaling.setEnabled(True)
            self.dlg_estimation_settings.groupBox_se_determination.setEnabled(True)
            QMessageBox.warning(self, 'Warning!', 'Unknown estimation method selected!')
            self.statusBar().showMessage(f"Unknown estimation method selected!")

    @pyqtSlot()
    def on_action_Load_Campaign_triggered(self):
        """Invoked whenever the menu item load campaign is pressed."""
        self.select_campaign_data_file_pickle()

    def select_campaign_data_file_pickle(self):
        """Launch file selection dialog to select a pickle file with saved campaign data."""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        pkl_file_filename, _ = QFileDialog.getOpenFileName(self,
                                                           'Select pkl file with campaign data',
                                                           DEFAULT_WORK_DIR_PATH,
                                                           "Pickle file (*.pkl)",
                                                           options=options)
        if pkl_file_filename:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            pkl_file_filename = QDir.toNativeSeparators(pkl_file_filename)
            # Add survey data to Campaign:
            try:
                self.load_campaign_from_pickle(pkl_file_filename, verbose=IS_VERBOSE)
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f"No campaign data loaded.")
            else:
                # Update GUI:

                # Enable/disable main menu items:
                self.menuAdd_Survey.setEnabled(True)
                self.menu_Add_Stations.setEnabled(True)
                self.action_Export_Results.setEnabled(True)
                self.action_Save_Campaign.setEnabled(True)
                self.action_Change_output_directory.setEnabled(True)
                self.action_Change_Campaign_name.setEnabled(True)

                # Set up view models and views for this campaign:
                # - Stations tab:
                self.set_up_station_view_model()
                self.enable_station_view_options_based_on_model()
                self.set_up_proxy_station_model()  # Time consuming!!!!!
                self.on_checkBox_filter_observed_stat_only_toggled(
                    state=self.checkBox_filter_observed_stat_only.checkState())
                self.update_comboBox_stations_selection_surrvey(survey_names=self.campaign.survey_names)
                self.update_stations_map(auto_range=True)
                # - Observations tab
                self.set_up_observation_view_model()
                self.enable_menu_observations_based_on_campaign_data()
                self.populate_survey_tree_widget()
                self.set_up_setup_view_model()
                if self.treeWidget_observations.topLevelItemCount() > 0:
                    self.treeWidget_observations.topLevelItem(0).setSelected(True)
                self.set_up_survey_view_model()
                # - Results tab:
                self.set_up_results_stations_view_model()
                self.set_up_results_correlation_matrix_view_model()
                self.set_up_results_observations_view_model()
                self.set_up_results_drift_view_model()
                self.set_up_results_vg_view_model()
                self.update_results_tab(select_latest_item=True)
                # - Corrections time series dialog:
                self.dlg_correction_time_series.check_correction_time_series_object(parent=self)
                self.dlg_correction_time_series.reset_update_gui()

                self.statusBar().showMessage(
                    f"Previously saved campaign loaded rom pickle file (name: {self.campaign.campaign_name}, "
                    f"output directory: {self.campaign.output_directory})")
                self.setWindowTitle('GravTools - Campaign: ' + self.campaign.campaign_name)
        else:
            self.statusBar().showMessage(f"No campaign data loaded.")

    def load_campaign_from_pickle(self, filename, verbose):
        """Load campaign data from pickle file."""
        # Load campaign object and replace previous one:
        if self.campaign is not None and verbose:
            print(f'Data of current campaign "{self.campaign.campaign_name}" will be overwritten!')
        self.campaign = Campaign.from_pkl(filename, verbose=verbose)

    @pyqtSlot()
    def on_action_Save_Campaign_triggered(self):
        """Invoked whenever the menu item save campaign is pressed."""
        self.save_campaign_to_pickle()

    def save_campaign_to_pickle(self):
        """Save campaign data (object) to a pickle file using the default path and filename."""
        try:
            filename = self.campaign.save_to_pickle(filename=None, verbose=IS_VERBOSE)
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            QMessageBox.information(self, 'Save Campaign Data', f'Campaign data saved as "{filename}".')

    @pyqtSlot()
    def on_spinBox_results_drift_plot_v_offset_value_changed(self):
        """Invoked whenever the value of the spin box changed."""
        self.update_drift_plot()

    def on_station_model_data_changed(self, topLeft, bottomRight, role):
        """Invoked, whenever data in the station view model changed."""
        # Did the "is_datum" flag change? => Update stations map!
        if not topLeft == bottomRight:
            QMessageBox.critical(self, 'Error!', 'Selection of multiple stations is not allowed!')
        else:
            is_datum_idx = self.station_model.get_data.columns.to_list().index('is_datum')
            if bottomRight.column() >= is_datum_idx and topLeft.column() <= is_datum_idx:
                self.update_stations_map(auto_range=False)
                # Set datum status ind the campaign data:
                station_record = self.station_model._data.iloc[topLeft.row()]
                if not (isinstance(station_record.is_datum, bool) or isinstance(station_record.is_datum, np.bool_)):
                    QMessageBox.critical(self, 'Error!', 'The is_datum flag is not a bool type!')
                else:
                    self.campaign.stations.set_datum_stations([station_record.station_name],
                                                              is_datum=station_record.is_datum, verbose=IS_VERBOSE)

    def set_up_stations_map(self):
        """Set up `self.GraphicsLayoutWidget_stations_map` widget."""
        self.glw_stations_map = self.GraphicsLayoutWidget_stations_map
        self.glw_stations_map.setBackground('w')  # white background color
        # Create sub-plots:
        self.stations_map = self.glw_stations_map.addPlot(0, 0, name='stations_map')
        self.stations_map.setLabel(axis='left', text='Latitude [°]')
        self.stations_map.setLabel(axis='bottom', text='Longitude [°]')
        self.stations_map.setTitle('')
        # self.stations_map.addLegend()

    def update_stations_map(self, auto_range=True):
        """Update the stations map in the stations tab.

        This method is used as slot. Hence, it will be invoked by signals from various GUI widgets that change the
        stations.

        Parameters
        ----------
        auto_range : bool, optional (deault = `True`)
            `True` indicates that the stations map is auto-ranged in order to view all content items.
        """
        SCATTER_PLOT_SYMBOL_SIZE = 10
        SCATTER_PLOT_PEN_WIDTH = 3
        SCATTER_PLOT_PEN_COLOR = 'k'
        SCATTER_PLOT_PEN_COLOR_DATUM = 'r'
        SCATTER_PLOT_BRUSH_COLOR_OBSERVED = QtGui.QColor('cyan')
        SCATTER_PLOT_BRUSH_COLOR_NOT_OBSERVED = 'g'
        STATION_LABEL_TEXT_SIZE = 10
        STATION_LABEL_TEXT_COLOR = 'k'

        self.stations_map.clear()
        self.stations_map.setTitle('')

        # Get list of stations from filter proxy model:
        station_name_list = []
        col_idx_stat_name = self.station_model.get_data.columns.to_list().index('station_name')
        for row_idx in range(self.proxy_station_model.rowCount()):
            station_name_list.append(self.proxy_station_model.index(row_idx, col_idx_stat_name).data())

        number_of_observed_stations = len(self.station_model.get_data.loc[self.station_model.get_data['is_observed']])
        number_of_stations_total = len(self.station_model.get_data)

        # Get stat_df from station view model and filter for stations in filter proxy model:
        filter_tmp = self.station_model.get_data['station_name'].isin(station_name_list)
        stat_df_filtered = self.station_model.get_data.loc[filter_tmp]

        # Plot all stations in filtered dataframe:
        # - Example: https://www.geeksforgeeks.org/pyqtgraph-different-colored-spots-on-scatter-plot-graph/
        scatter = pg.ScatterPlotItem()
        spots = []
        # - prep. data for scatterplot:
        for index, row in stat_df_filtered.iterrows():
            if row['is_observed']:
                brush_color = SCATTER_PLOT_BRUSH_COLOR_OBSERVED
            else:
                brush_color = SCATTER_PLOT_BRUSH_COLOR_NOT_OBSERVED
            if row['is_datum']:
                pen_color = SCATTER_PLOT_PEN_COLOR_DATUM
            else:
                pen_color = SCATTER_PLOT_PEN_COLOR
            spot_dic = {'pos': (row['long_deg'], row['lat_deg']),
                        'size': SCATTER_PLOT_SYMBOL_SIZE,
                        'pen': {'color': pen_color, 'width': SCATTER_PLOT_PEN_WIDTH},
                        'brush': brush_color,
                        'symbol': 'o'}
            spots.append(spot_dic)
        scatter.addPoints(spots)
        self.stations_map.addItem(scatter)

        # Add statione name labels:
        from PyQt5.QtGui import QPainterPath, QFont, QTransform

        # Show station name labels next to the scatter plot symbols, if enabled via the GUI:
        # See: https://pyqtgraph.readthedocs.io/en/latest/graphicsItems/scatterplotitem.html (symbol)
        # Example: https://www.geeksforgeeks.org/pyqtgraph-show-text-as-spots-on-scatter-plot-graph/
        if self.checkBox_stations_map_show_stat_name_labels.checkState() == Qt.Checked:
            spots = []
            # self.stations_map.addItem(scatter)
            for index, row in stat_df_filtered.iterrows():
                symbol = QtGui.QPainterPath()
                # creating QFont object
                f = QtGui.QFont()
                # setting font size
                f.setPointSize(STATION_LABEL_TEXT_SIZE)
                # adding text
                symbol.addText(STATION_LABEL_TEXT_SIZE, 0, f, row['station_name'])
                # getting bounding rectangle
                br = symbol.boundingRect()
                # getting scale
                # scale = min(1. / (br.width()), 1. / br.height())
                scale = 1. / (br.width() * 3)
                # getting transform object
                tr = QtGui.QTransform()
                # setting scale to transform object
                tr.scale(scale, scale)
                symbol_plot = tr.map(symbol)
                spot_dic = {'pos': (row['long_deg'], row['lat_deg']),
                            'size': STATION_LABEL_TEXT_SIZE / symbol_plot.boundingRect().height(),
                            # 'pen': {'color': pen_color, 'width': SCATTER_PLOT_PEN_WIDTH},
                            'brush': STATION_LABEL_TEXT_COLOR,
                            'symbol': symbol_plot}
                spots.append(spot_dic)
            scatter.addPoints(spots)

        # Title, grid, etc.:
        self.stations_map.setTitle(
            f'{len(station_name_list)} stations displayed ({number_of_observed_stations} observed; {number_of_stations_total} in total)')
        self.stations_map.showGrid(x=True, y=True)
        if auto_range:
            self.stations_map.autoRange()

    def on_checkBox_stations_map_show_stat_name_labels_state_changed(self):
        """Invoke, whenever the state of the checkbox changes."""
        self.update_stations_map(auto_range=False)

    def set_up_vg_plot_widget(self):
        """Set up `self.graphicsLayoutWidget_results_vg_plot` widget."""
        self.glw_vg_plot = self.graphicsLayoutWidget_results_vg_plot
        self.glw_vg_plot.setBackground('w')  # white background color
        # Create sub-plots:
        self.vg_plot = self.glw_vg_plot.addPlot(0, 0, name='vg_plot')
        self.vg_plot.setLabel(axis='left', text='Vertical gravity gradient [µGal/m]')
        self.vg_plot.setLabel(axis='bottom', text='Height over reference [m]')
        self.vg_plot.addLegend()
        self.vg_plot.showGrid(x=True, y=True)
        self.vg_plot.setTitle('')

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_vg_plot(self):
        """Update the vg plot in the results tab."""
        # Clear plot:
        self.vg_plot.clear()
        self.vg_plot.legend.clear()
        self.vg_plot.setTitle('')
        # Get GUI parameters:
        # - Selected LSM run:
        lsm_run_idx, lsm_run_time_str = self.get_selected_lsm_run()
        if lsm_run_idx != -1:
            # - Selected station:
            idx_selected_station, selected_station_name = self.get_selected_station()
            if selected_station_name == 'All stations':
                selected_station_names = None
            else:
                selected_station_names = [selected_station_name]
            # - Selected Survey:
            # idx_selected_survey, selected_survey_name = self.get_selected_survey()
            # if selected_survey_name == 'All surveys':
            #     selected_survey_names = None
            # else:
            #     selected_survey_names = [selected_survey_name]

            # Select lsm run:
            lsm_run = self.campaign.lsm_runs[lsm_run_idx]

            # Invoke plotting method according to the LSM method:
            if lsm_run.lsm_method == 'VG_LSM_nondiff':
                # Get plot settings from the GUI:
                plot_residuals = self.checkBox_results_vg_plot_show_residuals.checkState() == Qt.Checked
                if self.radioButton_results_vg_plot_details.isChecked() and not self.radioButton_results_vg_plot_full_polynomial.isChecked():
                    plot_type = 'detail'
                elif not self.radioButton_results_vg_plot_details.isChecked() and self.radioButton_results_vg_plot_full_polynomial.isChecked():
                    plot_type = 'full'
                else:
                    raise AssertionError(f'Invalid VG plot settings!')
                self.plot_vg_lsm_nondiff(lsm_run,
                                         plot_type=plot_type,
                                         plot_residuals=plot_residuals,
                                         stations=selected_station_names)
            else:
                self.vg_plot.clear()  # Clear vg plot

    def plot_vg_lsm_nondiff(self, lsm_run, plot_type='detail', plot_residuals=True, stations=None):
        """Create a VG plot for LSM estimation based on non-differential observations.

        Notes
        -----
        This method is applicable for the LSM methods VG_LSM_nondiff.

        Parameters
        ----------
        lsm_run : LSMNonDiff object.
            LSM object for VG estimation based non-differential observations.
        plot_type : str, optional (default=`detail`)
            Plot type selection. It `full` just the full vg polynomial (all degrees) is plotted. If `detail` the linear
            and the non-linear constituents are visualized separately in addition to post-fit residuals. Additionally,
            the mean setup heights are plotted.
        plot_residuals : booblean, optional (default=True)
            if `True`, the post-fit residuals are grouped per setup height (station) and plotted.
        stations : `None` (default) or list of station names (str)
            To filter for stations that will be displayed.
        """
        self.vg_plot.clear()
        self.vg_plot.legend.clear()
        self.vg_plot.setTitle('')

        res_plot_legend_str = ''

        if len(lsm_run.setups) != 1:  # Only one survey allowed!
            raise AssertionError(f'The current LSM run contains {len(lsm_run.setups)}. Only one survey allowed at '
                                 f'estimation of vertical gravity gradients!')

        # Get and prep. required data:
        vg_pol_df = lsm_run.vg_pol_df
        for survey_name, setup_data in lsm_run.setups.items():
            setup_df = setup_data['setup_df'].copy(deep=True)
            ref_epoch_delta_t_h = setup_data['ref_epoch_delta_t_h']
            ref_epoch_delta_t_campaign_h = setup_data['ref_epoch_delta_t_campaign_h']
        setup_obs_df = lsm_run.setup_obs_df.copy(deep=True)
        stat_obs_df = lsm_run.stat_obs_df.copy(deep=True)
        vg_polynomial_degree = lsm_run.vg_polynomial_degree
        vg_polynomial_ref_height_offset_m = lsm_run.vg_polynomial_ref_height_offset_m
        vg_pol_df = lsm_run.vg_pol_df

        # Get height range:
        h_max_m = settings.VG_PLOT_MAX_HEIGHT_M
        h_min_m = settings.VG_PLOT_MIN_HEIGHT_M
        if (h_max_m is None) or (stat_obs_df['dhf_sensor_max_m'].max() > h_max_m):
            h_max_m = stat_obs_df['dhf_sensor_max_m'].max() + settings.VG_PLOT_HEIGHT_DELTA_M
        if (h_min_m is None) or (stat_obs_df['dhf_sensor_min_m'].max() < h_min_m):
            h_min_m = stat_obs_df['dhf_sensor_min_m'].min() + settings.VG_PLOT_HEIGHT_DELTA_M
        h_m = np.linspace(h_min_m, h_max_m, settings.VG_PLOT_NUM_ITEMS_VG_POLYNOMIAL)

        # Evaluate the linear part of the vg polynomial (degree 2 to max.):
        coeff_list = [0.0] + vg_pol_df['coefficient'].to_list()
        del coeff_list[2:]  # Delete non-linear constituents
        coeff_list.reverse()
        vg_polynomial_linear_mugal = np.polyval(coeff_list, h_m)

        # Evaluate non-linear part of the vg polynomial (degree 2 to max.):
        if vg_polynomial_degree > 1:
            coeff_list = [0.0] + vg_pol_df['coefficient'].to_list()
            coeff_list[1] = 0.0  # Set linear constituent to 0
            coeff_list.reverse()
            vg_polynomial_nonlinear_mugal = np.polyval(coeff_list, h_m)
        else:  # Array of zeros
            vg_polynomial_nonlinear_mugal = np.zeros(settings.VG_PLOT_NUM_ITEMS_VG_POLYNOMIAL)

        # Full VG polynomial
        vg_polynomial_full_mugal = (vg_polynomial_linear_mugal + vg_polynomial_nonlinear_mugal)

        # Create plot and set title:
        if plot_type == 'full':
            self.vg_plot.setTitle(f'Full VG polynomial (degree = {vg_polynomial_degree})')
            self.vg_plot.setLabel(axis='left', text='Gravity [µGal]')

            # Plot VG polynomial
            pen = pg.mkPen(color='k', width=2)
            self.vg_plot.plot(h_m, vg_polynomial_full_mugal,
                                 name=f'VG polynomial',
                                 pen=pen)

        elif plot_type == 'detail':
            self.vg_plot.setTitle(f'Linear and non-linear components (pol. degree = {vg_polynomial_degree})')
            self.vg_plot.setLabel(axis='left', text='Gravity [µGal]')

            # Linear component (horizontal line):
            vg_linear = vg_pol_df.loc[vg_pol_df["degree"] == 1, "coefficient"].values[0]
            pen = pg.mkPen(color='k', width=2)
            self.vg_plot.plot([h_m.min(), h_m.max()], [0.0, 0.0],
                              name=f'LinearVG: {vg_linear:7.3f} µGal/m',
                              pen=pen)
            vg_linear_TextItem = pg.TextItem(text=f'{vg_linear:7.3f} µGal/m', color='k')
            text_width_px = vg_linear_TextItem.boundingRect().width()
            # TODO: For unknown reasons the command "self.vg_plot.getViewBox().viewPixelSize()" fails sometimes.
            #  Therefor the exception handling was added below.
            try:
                sx, _ = self.vg_plot.getViewBox().viewPixelSize()
                text_width_m = text_width_px * sx
                vg_linear_TextItem.setPos(h_m.max() - text_width_m, 0)
            except:
                vg_linear_TextItem.setPos(0, 0)
            self.vg_plot.addItem(vg_linear_TextItem)

            if vg_polynomial_degree == 1:
                # Plot post-fit residuals:
                if plot_residuals:
                    # Prep. residuals data for plotting:
                    setup_obs_df['res_plot_mugal'] = setup_obs_df['v_obs_est_mugal'].astype('float')
                    res_plot_legend_str = 'Post-fit residuals'

            elif vg_polynomial_degree > 1:
                pen = pg.mkPen(color='b', width=2)
                self.vg_plot.plot(h_m, vg_polynomial_nonlinear_mugal,
                                  name=f'Non-linear component',
                                  pen=pen)

                # Residuals w.r.t. estimated linear component of the VG (grouped by setup height):
                if plot_residuals:
                    self.vg_plot.setLabel(axis='left', text='Gravity [µGal], Residuals [µGal]')
                    # Calc. residuals w.r.t. linear VG:
                    mat_A = lsm_run.mat_A
                    mat_x = lsm_run.mat_x
                    mat_x_lin = mat_x[:-(vg_polynomial_degree - 1)]
                    mat_A_lin = mat_A[:, :-(vg_polynomial_degree - 1)]
                    setup_obs_df['g_obs_est_lin_vg_mugal'] = mat_A_lin @ mat_x_lin
                    setup_obs_df['res_lin_vg_mugal'] = setup_obs_df['g_obs_mugal'] - setup_obs_df['g_obs_est_lin_vg_mugal']
                    setup_obs_df['res_plot_mugal'] = setup_obs_df['res_lin_vg_mugal']
                    res_plot_legend_str = 'Residuals w.r.t. linear VG'
            else:
                raise AssertionError(f'Invalid degree of the VG polynomial: {vg_polynomial_degree}')

            # Plot residuals:
            if plot_residuals:
                scatter = pg.ScatterPlotItem()
                spots = []
                # - Prep. data for scatterplot:
                for index, row in setup_obs_df.iterrows():
                    if stations is None:
                        brush_color = self.station_colors_dict_results[row['station_name']]
                    else:
                        if row['station_name'] not in stations:
                            brush_color = 'w'
                        else:
                            brush_color = self.station_colors_dict_results[row['station_name']]
                    spot_dic = {'pos': (row['dhf_sensor_m'], row['res_plot_mugal']),
                                'size': settings.VG_PLOT_SCATTER_PLOT_SYMBOL_SIZE,
                                'pen': {'color': settings.VG_PLOT_SCATTER_PLOT_PEN_COLOR,
                                        'width': settings.VG_PLOT_SCATTER_PLOT_PEN_WIDTH},
                                'brush': brush_color}
                    spots.append(spot_dic)
                scatter.addPoints(spots)
                self.vg_plot.addItem(scatter)

                # Add residuals to legend:
                # - https://pyqtgraph.readthedocs.io/en/latest/graphicsItems/legenditem.html
                for station, color in self.station_colors_dict_results.items():
                    if (stations is None) or (station in stations):
                        s_item_tmp = pg.ScatterPlotItem()
                        s_item_tmp.setBrush(color)
                        s_item_tmp.setPen({'color': settings.VG_PLOT_SCATTER_PLOT_PEN_COLOR,
                                           'width': settings.VG_PLOT_SCATTER_PLOT_PEN_WIDTH})
                        s_item_tmp.setSize(settings.VG_PLOT_SCATTER_PLOT_SYMBOL_SIZE)
                        self.vg_plot.legend.addItem(s_item_tmp, res_plot_legend_str + f' ({station})')

                # Plot mean residuals w.r.t. linear VG and error bars per station:
                # - Calc. mean residuals per setup height (station) and their standard deviation:
                tmp_df = setup_obs_df.loc[:, ['station_name', 'res_plot_mugal']].groupby(
                    'station_name').mean().rename(
                    columns={"res_plot_mugal": "mean_res_plot_mugal"})
                stat_obs_df = stat_obs_df.merge(tmp_df, left_on='station_name', right_on='station_name', how='left')
                tmp_df = setup_obs_df.loc[:, ['station_name', 'res_plot_mugal']].groupby(
                    'station_name').std().rename(
                    columns={"res_plot_mugal": "std_res_plot_mugal"})
                stat_obs_df = stat_obs_df.merge(tmp_df, left_on='station_name', right_on='station_name', how='left')
                scatter_mean_res = pg.ScatterPlotItem()
                spots = []
                # - Plot data:
                for index, row in stat_obs_df.iterrows():
                    spot_dic = {'pos': (
                        row['dhf_sensor_mean_m'], row['mean_res_plot_mugal']),
                        'size': 20,
                        'pen': {'color': 'k',
                                'width': 2},
                        'symbol': 'x',
                        'brush': 'k'}
                    spots.append(spot_dic)
                scatter_mean_res.addPoints(spots)
                self.vg_plot.addItem(scatter_mean_res)

                pen = pg.mkPen(color='k', width=1)
                error_bar = pg.ErrorBarItem(x=stat_obs_df['dhf_sensor_mean_m'].to_numpy(),
                                            y=stat_obs_df['mean_res_plot_mugal'].to_numpy(),
                                            height=stat_obs_df['std_res_plot_mugal'].to_numpy() * 2,
                                            beam=0.1,
                                            pen=pen)
                self.vg_plot.addItem(error_bar)

                # Add item to legend
                s_item_tmp = pg.ScatterPlotItem()
                s_item_tmp.setBrush('k')
                s_item_tmp.setSymbol('x')
                s_item_tmp.setSize(20)
                s_item_tmp.setPen({'color': 'k',
                                   'width': 2})
                self.vg_plot.legend.addItem(s_item_tmp, f'Mean Residuals + Errorbars')

        else:
            raise AssertionError(f'Unknown VG plot type: {plot_type}!')

        # Check/set min/mx Y-range:
        self.vg_plot.autoRange()
        [y_min, y_max] = self.vg_plot.axes['left']['item'].range
        flag_set_y_range = False
        if y_min > settings.VG_PLOT_MIN_LOWER_L_RANGE:
            y_min = settings.VG_PLOT_MIN_LOWER_L_RANGE
            flag_set_y_range = True
        if y_max < settings.VG_PLOT_MIN_UPPER_Y_RANGE:
            y_max = settings.VG_PLOT_MIN_UPPER_Y_RANGE
            flag_set_y_range = True
        if flag_set_y_range:
            self.vg_plot.setYRange(y_min, y_max)

        # Plot mean setup heights:
        for index, row in stat_obs_df.iterrows():
            if stations is None:
                brush_color = self.station_colors_dict_results[row['station_name']]
            else:
                if row['station_name'] not in stations:
                    brush_color = '#bebebe'  # Grey, HEX-codd
                else:
                    brush_color = self.station_colors_dict_results[row['station_name']]
            pen = pg.mkPen(color=brush_color, width=1, style=Qt.DashLine)
            line = self.vg_plot.plot([row['dhf_sensor_mean_m'], row['dhf_sensor_mean_m']], [y_min, y_max],
                                     pen=pen,
                                     name=f'Mean height: {row["station_name"]}')

        # Plot reference height
        pen = pg.mkPen(color='k', width=1.5, style=Qt.DashDotDotLine)
        line = self.vg_plot.plot([vg_polynomial_ref_height_offset_m, vg_polynomial_ref_height_offset_m], [y_min, y_max],
                                 pen=pen,
                                 name=f'Ref. height: {vg_polynomial_ref_height_offset_m:5.3f} m')

        # General plot settings:
        # - Get and set X-range:
        self.vg_plot.setXRange(h_m.min(), h_m.max())
        # self.vg_plot.autoRange()

    def set_up_drift_plot_widget(self):
        """Set up `self.graphicsLayoutWidget_results_drift_plot` widget."""
        self.glw_drift_plot = self.graphicsLayoutWidget_results_drift_plot
        self.glw_drift_plot.setBackground('w')  # white background color
        # Create sub-plots:
        self.drift_plot = self.glw_drift_plot.addPlot(0, 0, name='drift_plot',
                                                      axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.drift_plot.setLabel(axis='left', text='')
        self.drift_plot.addLegend()

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_drift_plot(self):
        """Update the drift plot in the results tab.

        This method is used as slot. Hence, it will be invoked by signals from various GUI widgets that change the
        drift plot in the according plotting widget.
        """
        # Clear plot:
        self.drift_plot.clear()
        self.drift_plot.legend.clear()
        self.drift_plot.setTitle('')
        # Get GUI parameters:
        # - Selected LSM run:
        lsm_run_idx, lsm_run_time_str = self.get_selected_lsm_run()
        if lsm_run_idx != -1:
            # - Selected station:
            idx_selected_station, selected_station_name = self.get_selected_station()
            if selected_station_name == 'All stations':
                selected_station_names = None
            else:
                selected_station_names = [selected_station_name]
            # - Selected Survey:
            idx_selected_survey, selected_survey_name = self.get_selected_survey_results()
            if selected_survey_name == 'All surveys':
                selected_survey_names = None
            else:
                selected_survey_names = [selected_survey_name]

            # Select lsm run:
            lsm_run = self.campaign.lsm_runs[lsm_run_idx]

            # - Enable/disable plot settings groupBox:
            if lsm_run.lsm_method == 'LSM_diff':
                # self.groupBox_results_drift_plot.setEnabled(True)
                self.spinBox_results_drift_plot_v_offset.setEnabled(True)
                self.label_results_drift_plot_v_offset.setEnabled(True)
            else:
                # self.groupBox_results_drift_plot.setEnabled(False)
                self.spinBox_results_drift_plot_v_offset.setEnabled(False)
                self.label_results_drift_plot_v_offset.setEnabled(False)

            if lsm_run.lsm_method == 'LSM_diff':
                offset_mugal = self.spinBox_results_drift_plot_v_offset.value()
                self.plot_drift_lsm_diff(lsm_run, surveys=selected_survey_names,
                                         stations=selected_station_names,
                                         offset_user_defined_mugal=offset_mugal)
            elif lsm_run.lsm_method == 'MLR_BEV':
                self.plot_drift_mlr_bev_legacy(lsm_run, surveys=selected_survey_names, stations=selected_station_names)
            elif lsm_run.lsm_method == 'LSM_non_diff':
                self.plot_drift_lsm_nondiff(lsm_run, surveys=selected_survey_names,
                                            stations=selected_station_names)
            elif lsm_run.lsm_method == 'VG_LSM_nondiff':
                self.plot_drift_lsm_nondiff(lsm_run, surveys=selected_survey_names,
                                            stations=selected_station_names)
            else:
                self.drift_plot.clear()  # Clear drift plot

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def plot_drift_lsm_nondiff(self, lsm_run, surveys=None, stations=None):
        """Create a drift plot for LSM estimation based on non-differential observations.

        The drift plot shows the evaluated drift polynomial (estimated). Additionally the post-fit residuals are plotted
        w.r.t. the drift polynomial, i.e. estimated drift at observation epochs - residuals at corresponding epochs.

        Notes
        -----
        This method is applicable for the LSM methods VG_LSM_nondiff and LSM_non_diff.

        Parameters
        ----------
        lsm_run : LSMNonDiff object.
            LSM object for VG estimation based non-differential observations.
        surveys : `None` (default) or list of survey names (str)
            To filter for surveys that will be displayed.
        stations : `None` (default) or list of station names (str)
            To filter for stations that will be displayed.
        """
        self.drift_plot.clear()
        self.drift_plot.legend.clear()

        drift_pol_df = lsm_run.drift_pol_df

        plotted_stations = []

        # Loop over surveys (setup data) in the selected lsm run object and plot data:
        for survey_name, setup_data in lsm_run.setups.items():
            setup_df_orig = setup_data['setup_df']
            # Filter for surveys:
            if surveys is not None:
                if survey_name not in surveys:
                    continue

            # Prep data:
            drift_pol_df_short = drift_pol_df.loc[drift_pol_df['survey_name'] == survey_name]
            setup_df = setup_df_orig.copy(deep=True)  # Make hard copy to protect original data!
            # stat_obs_df_short = stat_obs_df.loc[:, ['station_name', 'g_est_mugal', 'sd_g_est_mugal']]
            # setup_df = pd.merge(setup_df, stat_obs_df_short, on='station_name')
            # setup_df['g_plot_mugal'] = setup_df['g_mugal'] - setup_df['g_est_mugal']
            setup_df.sort_values(by='delta_t_campaign_h', inplace=True)
            # Merge df => Colum with post-fit residuals ("v_obs_est_mugal") added to "setup_df":
            if not hasattr(setup_data, 'setup_calc_method'):  # to grand downward compatibility < v0.2.0
                setup_calc_method = 'variance_weighted_mean'
            else:
                setup_calc_method = setup_data['setup_calc_method']
            if setup_calc_method != 'individual_obs':
                # Merge on the reference epochs, because the setup-IDs are not unique!
                setup_obs_df_short = lsm_run.setup_obs_df.loc[
                    lsm_run.setup_obs_df['survey_name'] == survey_name, ['ref_epoch_dt', 'v_obs_est_mugal']].copy(deep=True)
                setup_df = pd.merge(setup_df, setup_obs_df_short, how='left', left_on='epoch_dt', right_on='ref_epoch_dt')  # WEG!
            else:
                setup_obs_df_short = lsm_run.setup_obs_df.loc[
                    lsm_run.setup_obs_df['survey_name'] == survey_name, ['ref_epoch_dt', 'v_obs_est_mugal']].copy(
                    deep=True)
                setup_df = pd.merge(setup_df, setup_obs_df_short, how='left', left_on='epoch_dt',
                                    right_on='ref_epoch_dt')
            # Evaluate drift polynomial:
            coeff_list = drift_pol_df_short['coefficient'].to_list()
            coeff_list.reverse()
            if lsm_run.drift_ref_epoch_type == 'survey':
                delta_t_min_h = setup_df['delta_t_h'].min()  # = 0
                delta_t_max_h = setup_df['delta_t_h'].max()
            elif lsm_run.drift_ref_epoch_type == 'campaign':
                delta_t_min_h = setup_df['delta_t_campaign_h'].min()  # = 0
                delta_t_max_h = setup_df['delta_t_campaign_h'].max()
            delta_t_h = np.linspace(delta_t_min_h, delta_t_max_h, settings.DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION)
            drift_polynomial_mugal = np.polyval(coeff_list, delta_t_h)

            # Evaluate Drift for observation epochs:
            if lsm_run.drift_ref_epoch_type == 'survey':
                delta_t_h = setup_df['delta_t_h']
            elif lsm_run.drift_ref_epoch_type == 'campaign':
                delta_t_h = setup_df['delta_t_campaign_h']
            drift_at_obs_epochs_mugal = np.polyval(coeff_list, delta_t_h)
            setup_df['drift_at_obs_epochs_mugal'] = drift_at_obs_epochs_mugal  # Add to df

            # Drift function time reference as UNIX time (needed for plots):
            epoch_unix_min = setup_df['epoch_unix'].min()
            epoch_unix_max = setup_df['epoch_unix'].max()
            delta_t_epoch_unix = np.linspace(epoch_unix_min, epoch_unix_max,
                                             settings.DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION)

            yy_mugal = drift_polynomial_mugal

            # Constant offset to be subtracted from y-axis for better readability:
            subtr_const_mugal = round(yy_mugal.mean()/ 1000) * 1000

            # Plot drift function:
            pen = pg.mkPen(color='k', width=2)
            self.drift_plot.plot(delta_t_epoch_unix, yy_mugal - subtr_const_mugal,
                                 name=f'drift: {survey_name}',
                                 pen=pen, symbol='o', symbolSize=4, symbolBrush='k')

            # Plot observation data (setup observations):
            # - Example: https://www.geeksforgeeks.org/pyqtgraph-different-colored-spots-on-scatter-plot-graph/
            scatter = pg.ScatterPlotItem()
            spots = []
            # - prep. data for scatterplot:
            for index, row in setup_df.iterrows():
                if stations is None:
                    brush_color = self.station_colors_dict_results[row['station_name']]
                else:
                    if row['station_name'] not in stations:
                        brush_color = 'w'
                    else:
                        brush_color = self.station_colors_dict_results[row['station_name']]
                spot_dic = {'pos': (row['epoch_unix'], row['drift_at_obs_epochs_mugal'] - row['v_obs_est_mugal'] - subtr_const_mugal),
                            'size': settings.DRIFT_PLOT_SCATTER_PLOT_SYMBOL_SIZE,
                            'pen': {'color': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_COLOR,
                                    'width': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_WIDTH},
                            'brush': brush_color}
                spots.append(spot_dic)
                plotted_stations.append(row['station_name'])

            scatter.addPoints(spots)
            self.drift_plot.addItem(scatter)

        # Add station items to legend:
        plotted_stations = unique_ordered_list(plotted_stations)
        for station in plotted_stations:
            color = self.station_colors_dict_results[station]
            s_item_tmp = pg.ScatterPlotItem()
            s_item_tmp.setBrush(color)
            s_item_tmp.setPen({'color': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_COLOR,
                               'width': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_WIDTH})
            s_item_tmp.setSize(settings.DRIFT_PLOT_SCATTER_PLOT_SYMBOL_SIZE)
            self.drift_plot.legend.addItem(s_item_tmp, station)

        # Adjust plot window:
        self.drift_plot.showGrid(x=True, y=True)
        self.drift_plot.setLabel(axis='left', text=f'g [µGal] + {subtr_const_mugal / 1000:.1f} mGal')
        self.drift_plot.setTitle(f'Drift function w.r.t. setup observations')
        self.drift_plot.autoRange()

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def plot_drift_lsm_diff(self, lsm_run, surveys=None, stations=None, offset_user_defined_mugal=0):
        """Create a drift plot for LSM runs based on differential observations (method: LSM_diff)

        Parameters
        -----------
        lsm_run : LSMDiff object
            LSM object for differential observations.
        surveys : `None` (default) or list of survey names (str)
            To filter for surveys that will be displayed.
        stations : `None` (default) or list of station names (str)
            To filter for stations that will be displayed.
        offset_user_defined_mugal : int, optional (default=0)
            User defined offset between the drift polynomial function and the observed data points.
        """
        self.drift_plot.clear()
        self.drift_plot.legend.clear()

        plotted_stations = []

        stat_obs_df = lsm_run.stat_obs_df
        drift_pol_df = lsm_run.drift_pol_df
        if stat_obs_df is not None and drift_pol_df is not None:
            # Loop over surveys (setup data) in the selected lsm run object and plot data:
            for survey_name, setup_data in lsm_run.setups.items():
                setup_df_orig = setup_data['setup_df']
                # Filter for surveys:
                if surveys is not None:
                    if survey_name not in surveys:
                        continue

                # Prep data:
                drift_pol_df_short = drift_pol_df.loc[drift_pol_df['survey_name'] == survey_name]
                setup_df = setup_df_orig.copy(deep=True)  # Make hard copy to protect original data!
                stat_obs_df_short = stat_obs_df.loc[:, ['station_name', 'g_est_mugal', 'sd_g_est_mugal']]
                setup_df = pd.merge(setup_df, stat_obs_df_short, on='station_name')
                setup_df['g_plot_mugal'] = setup_df['g_mugal'] - setup_df['g_est_mugal']
                setup_df.sort_values(by='delta_t_campaign_h', inplace=True)

                # Evaluate drift polynomial:
                coeff_list = drift_pol_df_short['coefficient'].to_list()
                coeff_list.reverse()
                coeff_list.append(0)
                if lsm_run.drift_ref_epoch_type == 'survey':
                    delta_t_min_h = setup_df['delta_t_h'].min()  # = 0
                    delta_t_max_h = setup_df['delta_t_h'].max()
                elif lsm_run.drift_ref_epoch_type == 'campaign':
                    delta_t_min_h = setup_df['delta_t_campaign_h'].min()  # = 0
                    delta_t_max_h = setup_df['delta_t_campaign_h'].max()
                delta_t_h = np.linspace(delta_t_min_h, delta_t_max_h, settings.DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION)
                drift_polynomial_mugal = np.polyval(coeff_list, delta_t_h)

                # Drift function time reference as UNIX time (needed for plots):
                epoch_unix_min = setup_df['epoch_unix'].min()
                epoch_unix_max = setup_df['epoch_unix'].max()
                delta_t_epoch_unix = np.linspace(epoch_unix_min, epoch_unix_max,
                                                 settings.DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION)

                # !!! Due to the differential observations, the constant bias (N0) of the gravity reading cannot be estimated!
                # In order to draw the drift polynomial function w.r.t. the gravity meter observations (for the sake of visual
                # assessment of the drift function), the const. bias N0 is approximated, see below.
                offset_mugal = setup_df['g_plot_mugal'].mean() - drift_polynomial_mugal.mean()
                offset_mugal = offset_mugal + offset_user_defined_mugal
                yy_mugal = drift_polynomial_mugal + offset_mugal

                # Constant to be subtracted from y-axis:
                subtr_const_mugal = round(setup_df['g_plot_mugal'].mean() / 1000) * 1000

                # Plot drift function:
                pen = pg.mkPen(color='k', width=2)
                self.drift_plot.plot(delta_t_epoch_unix, yy_mugal - subtr_const_mugal,
                                     name=f'drift: {survey_name} (offset: {offset_mugal:.1f} µGal)',
                                     pen=pen, symbol='o', symbolSize=4, symbolBrush='k')

                # plot observation data (setup observations):
                # - Example: https://www.geeksforgeeks.org/pyqtgraph-different-colored-spots-on-scatter-plot-graph/
                scatter = pg.ScatterPlotItem()
                spots = []
                # - prep. data for scatterplot:
                for index, row in setup_df.iterrows():
                    if stations is None:
                        brush_color = self.station_colors_dict_results[row['station_name']]
                    else:
                        if row['station_name'] not in stations:
                            brush_color = 'w'
                        else:
                            brush_color = self.station_colors_dict_results[row['station_name']]
                    spot_dic = {'pos': (row['epoch_unix'], row['g_plot_mugal'] - subtr_const_mugal),
                                'size': settings.DRIFT_PLOT_SCATTER_PLOT_SYMBOL_SIZE,
                                'pen': {'color': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_COLOR,
                                        'width': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_WIDTH},
                                'brush': brush_color}
                    spots.append(spot_dic)
                    plotted_stations.append(row['station_name'])

                scatter.addPoints(spots)
                self.drift_plot.addItem(scatter)

            # Add station items to legend:
            # - https://pyqtgraph.readthedocs.io/en/latest/graphicsItems/legenditem.html
            plotted_stations = unique_ordered_list(plotted_stations)
            for station in plotted_stations:
                color = self.station_colors_dict_results[station]
                s_item_tmp = pg.ScatterPlotItem()
                s_item_tmp.setBrush(color)
                s_item_tmp.setPen({'color': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_COLOR,
                                   'width': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_WIDTH})
                s_item_tmp.setSize(settings.DRIFT_PLOT_SCATTER_PLOT_SYMBOL_SIZE)
                self.drift_plot.legend.addItem(s_item_tmp, station)

            # Adjust plot window:
            self.drift_plot.showGrid(x=True, y=True)
            self.drift_plot.setLabel(axis='left', text=f'g [µGal] + {subtr_const_mugal / 1000:.1f} mGal')
            self.drift_plot.setTitle(f'Drift function w.r.t. setup observations (with arbitrary offset!)')
            self.drift_plot.autoRange()

    def plot_drift_mlr_bev_legacy(self, lsm_run, surveys=None, stations=None):
        """Create a drift plot for LSM runs using multiple linear regression (method: MLR BEV legacy)

        Parameters
        -----------
        lsm_run : BEVLegacyProcessing object
            Parameter adjustment object.
        surveys : `None` (default) or list of survey names (str)
            To filter for surveys that will be displayed.
        stations : `None` (default) or list of station names (str)
            To filter for stations that will be displayed.
        """
        self.drift_plot.clear()
        self.drift_plot.legend.clear()

        plotted_stations = []

        stat_obs_df = lsm_run.stat_obs_df
        drift_pol_df = lsm_run.drift_pol_df

        # Loop over surveys (setup data) in the selected lsm run object and plot data:
        flag_first_survey = True
        for survey_name, setup_data in lsm_run.setups.items():
            setup_df_orig = setup_data['setup_df']
            if flag_first_survey:
                flag_first_survey = False
                plot_setup_df = pd.DataFrame(columns=setup_df_orig.columns)
            # Filter for surveys:
            if surveys is not None:
                if survey_name not in surveys:
                    continue
            # Prep data:
            setup_df = setup_df_orig.copy(deep=True)  # Make hard copy to protect original data!
            stat_obs_df_short = stat_obs_df.loc[:, ['station_name', 'g_drift_est_mugal']]
            setup_df = pd.merge(setup_df, stat_obs_df_short, on='station_name')
            setup_df['g_plot_mugal'] = setup_df['g_mugal'] - setup_df['g_drift_est_mugal']
            setup_df.sort_values(by='delta_t_campaign_h', inplace=True)
            plot_setup_df = pd.concat([plot_setup_df, setup_df])

        plot_setup_df.sort_values(by='delta_t_campaign_h', inplace=True)

        # Evaluate drift function:
        coeff_list = drift_pol_df['coefficient'].to_list()
        coeff_list.reverse()
        coeff_list.append(0)
        delta_t_min_h = plot_setup_df['delta_t_campaign_h'].min()  # = 0
        delta_t_max_h = plot_setup_df['delta_t_campaign_h'].max()
        delta_t_h = np.linspace(delta_t_min_h, delta_t_max_h, settings.DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION)
        drift_polynomial_mugal = np.polyval(coeff_list, delta_t_h)

        # Drift function time reference as UNIX time (needed for plots):
        epoch_unix_min = plot_setup_df['epoch_unix'].min()
        epoch_unix_max = plot_setup_df['epoch_unix'].max()
        delta_t_epoch_unix = np.linspace(epoch_unix_min, epoch_unix_max,
                                         settings.DRIFT_PLOT_NUM_ITEMS_IN_DRIFT_FUNCTION)

        # Plot drift function:
        pen = pg.mkPen(color='k', width=2)
        self.drift_plot.plot(delta_t_epoch_unix, drift_polynomial_mugal,
                             name=f'drift polynomial',
                             pen=pen, symbol='o', symbolSize=4, symbolBrush='k')

        # plot observation data (setup observations):
        # - Example: https://www.geeksforgeeks.org/pyqtgraph-different-colored-spots-on-scatter-plot-graph/
        scatter = pg.ScatterPlotItem()
        spots = []
        # - prep. data for scatterplot:
        for index, row in plot_setup_df.iterrows():
            if stations is None:
                brush_color = self.station_colors_dict_results[row['station_name']]
            else:
                if row['station_name'] not in stations:
                    brush_color = 'w'
                else:
                    brush_color = self.station_colors_dict_results[row['station_name']]
            spot_dic = {'pos': (row['epoch_unix'], row['g_plot_mugal']),
                        'size': settings.DRIFT_PLOT_SCATTER_PLOT_SYMBOL_SIZE,
                        'pen': {'color': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_COLOR,
                                'width': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_WIDTH},
                        'brush': brush_color}
            spots.append(spot_dic)
            plotted_stations.append(row['station_name'])

        scatter.addPoints(spots)
        self.drift_plot.addItem(scatter)

        # Add station items to legend:
        # - https://pyqtgraph.readthedocs.io/en/latest/graphicsItems/legenditem.html
        plotted_stations = unique_ordered_list(plotted_stations)
        for station in plotted_stations:
            color = self.station_colors_dict_results[station]
            s_item_tmp = pg.ScatterPlotItem()
            s_item_tmp.setBrush(color)
            s_item_tmp.setPen({'color': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_COLOR,
                               'width': settings.DRIFT_PLOT_SCATTER_PLOT_PEN_WIDTH})
            s_item_tmp.setSize(settings.DRIFT_PLOT_SCATTER_PLOT_SYMBOL_SIZE)
            self.drift_plot.legend.addItem(s_item_tmp, station)

        # Adjust plot window:
        self.drift_plot.showGrid(x=True, y=True)
        self.drift_plot.setLabel(axis='left', text=f'g [µGal]')
        self.drift_plot.setTitle(f'Drift function w.r.t. setup observations')
        self.drift_plot.autoRange()

    def set_up_obseration_results_plots_hist_method_comboBox(self):
        """Set up the histogram method combo box."""
        self.comboBox_results_obs_plot_hist_method.addItems(list(settings.NUMPY_HISTOGRAM_BIN_EDGES_OPTIONS.keys()))
        # self.histogram_bin_method_selection(0)  # default: i. item in dict

    def set_up_obseration_results_plots_widget(self):
        """Set up `self.graphicsLayoutWidget_results_observations_plots`."""
        self.glw_obs_results = self.graphicsLayoutWidget_results_observations_plots
        self.glw_obs_results.setBackground('w')  # white background color

        # Create sub-plots:
        # - Initialize time series or histogram:
        if self.radioButton_results_obs_plot_timeseries.isChecked():
            self.init_observation_results_plots_timerseries()
        elif self.radioButton_results_obs_plot_histogram.isChecked():
            self.init_observation_results_plots_histogram()
        else:
            pass

    def init_observation_results_plots_timerseries(self):
        """Initialize the time series plot for observation results."""
        if self.glw_obs_results.getItem(0,0) is not None:
            self.glw_obs_results.removeItem(self.glw_obs_results.getItem(0, 0))
        self.plot_obs_results = self.glw_obs_results.addPlot(0, 0, name='obs_results',
                                                             axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_results.setLabel(axis='left', text='')
        self.plot_obs_results.addLegend()

    def init_observation_results_plots_histogram(self):
        """Initialize the histogram plot for observation results."""
        if self.glw_obs_results.getItem(0,0) is not None:
            self.glw_obs_results.removeItem(self.glw_obs_results.getItem(0, 0))
        self.plot_obs_results = self.glw_obs_results.addPlot(0, 0, name='obs_results')
        self.plot_obs_results.setLabel(axis='left', text='')

    def get_hist_bin_method(self):
        """Get selected bin method from GUI."""
        bins = self.comboBox_results_obs_plot_hist_method.currentText()
        if bins == 'Num. of bins':
            bins = self.spinBox_results_obs_plot_number_bins.value()
        return bins

    def on_histogram_bin_method_currentIndexChanged(self):
        """Select bin determination method."""
        hist_bin_method = self.comboBox_results_obs_plot_hist_method.currentText()
        self.comboBox_results_obs_plot_hist_method.setToolTip(
            settings.NUMPY_HISTOGRAM_BIN_EDGES_OPTIONS[hist_bin_method])
        if hist_bin_method == 'Num. of bins':
            self.label_results_obs_plot_number_bins.setEnabled(True)
            self.spinBox_results_obs_plot_number_bins.setEnabled(True)
        else:
            self.label_results_obs_plot_number_bins.setEnabled(False)
            self.spinBox_results_obs_plot_number_bins.setEnabled(False)

    def plot_observation_results(self, results_obs_df=None, column_name='', type='timeseries'):
        """Plots observation data to the GraphicsLayoutWidget.

        Notes
        -----
        If input parameters `` or/and `` is/are `None` or '', the plot content is deleted and the plot is reseted.
        """
        # Clear plot in any case:
        self.plot_obs_results.clear()
        # Plot time series:
        if self.radioButton_results_obs_plot_timeseries.isChecked():
            self.init_observation_results_plots_timerseries()
            if results_obs_df is not None:  # Data available for plotting
                # Get data:
                data = results_obs_df[column_name].values
                if isinstance(data, object):
                    data = data.astype(float)
                obs_epoch_timestamps = (results_obs_df['ref_epoch_dt'].values - np.datetime64(
                    '1970-01-01T00:00:00')) / np.timedelta64(1, 's')
                plot_name = self.results_observation_model.get_short_column_description(column_name)
                self.plot_xy_data(self.plot_obs_results, obs_epoch_timestamps, data, plot_name=plot_name, color='b',
                                  symbol='o', symbol_size=10)
                self.plot_obs_results.showGrid(x=True, y=True)
                column_description = self.results_observation_model.get_plotable_columns()[column_name]
                self.plot_obs_results.setLabel(axis='left', text=column_description)
                self.plot_obs_results.autoRange()

        # Plot histogram:
        elif self.radioButton_results_obs_plot_histogram.isChecked():
            self.init_observation_results_plots_histogram()
            if results_obs_df is not None:  # Data available for plotting
                # Get bin method:
                bins = self.get_hist_bin_method()
                # Get data:
                data = results_obs_df[column_name].values
                if isinstance(data, object):
                    data = data.astype(float)
                y, x = np.histogram(data, bins=bins)
                self.plot_obs_results.plot(x, y, stepMode=True, fillLevel=0, brush=(0, 0, 255, 80))
                self.plot_obs_results.showGrid(x=True, y=True)

                # Add legend with statistics:
                mean = np.mean(data)
                median = np.median(data)
                std = np.std(data)
                q25 = np.quantile(data, 0.25)
                q75 = np.quantile(data, 0.75)
                iqr = q75 - q25
                legend = self.plot_obs_results.addLegend()
                legend.setLabelTextColor('k')
                style1 = pg.PlotDataItem(pen='w')
                style2 = pg.PlotDataItem(pen='w')
                style3 = pg.PlotDataItem(pen='w')
                style4 = pg.PlotDataItem(pen='w')
                legend.addItem(style1, f'Mean = {mean:0.3f}')
                legend.addItem(style2, f'Std. = {std:0.3f}')
                legend.addItem(style3, f'Median = {median:0.3f}')
                legend.addItem(style4, f'IQR = {iqr:0.3f}')

                self.plot_obs_results.autoRange()
        else:
            pass

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_results_obs_plots(self, idx=0):
        """Update the observation results plots in the results tab.

        Notes
        -----
        The input argument `idx` was added in order to enable logging execution times with the `@conditional_decorator(time_it, settings.DEBUG_TIME_IT)` decorator. As
        this method is also used as pyqt slot, it would otherwise raise an error.
        """
        # Get data from selected column
        col_idx, column_name = self.get_selected_obs_data_column()
        if col_idx != -1:
            filtered_results_obs_df = self.results_observation_model.get_model_data_df
        else:  # Invalid indices
            filtered_results_obs_df = None
        try:
            self.plot_observation_results(filtered_results_obs_df, column_name)
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))

    def set_up_results_vg_view_model(self):
        """Set up the view model for the VG results table view."""
        try:
            self.results_vg_model = ResultsVGModel(self.campaign.lsm_runs)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No LSM-adjustment results available!')
            self.statusBar().showMessage(f"No LSM-adjustment results available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.tableView_results_vg_table.setModel(self.results_vg_model)
            self.tableView_results_vg_table.resizeColumnsToContents()

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_results_vg_table_view(self, lsm_run_index: int):
        """Update the VG results table view after changing the table model."""
        self.results_vg_model.load_lsm_runs(self.campaign.lsm_runs)
        self.results_vg_model.update_view_model(lsm_run_index)
        self.results_vg_model.layoutChanged.emit()  # Show changes in table view
        self.tableView_results_vg_table.resizeColumnsToContents()

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

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_results_drift_table_view(self, lsm_run_index: int, survey_name=None):
        """Update the drift results table view after changing the table model."""
        self.results_drift_model.load_lsm_runs(self.campaign.lsm_runs)
        self.results_drift_model.update_view_model(lsm_run_index, survey_name=survey_name)
        self.results_drift_model.layoutChanged.emit()  # Show changes in table view
        # self.tableView_results_drift.resizeColumnsToContents()
        resize_table_view_columns(table_view=self.tableView_results_drift,
                                  n=1000,
                                  add_pixel=10)

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

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_results_observation_table_view(self, lsm_run_index: int, station_name=None, survey_name=None):
        """Update the observation results table view after changing the table model."""
        self.results_observation_model.load_lsm_runs(self.campaign.lsm_runs)
        self.results_observation_model.update_view_model(lsm_run_index,
                                                         station_name=station_name,
                                                         survey_name=survey_name,
                                                         gui_simple_mode=self.dlg_options.gui_simple_mode)
        self.results_observation_model.layoutChanged.emit()  # Show changes in table view
        # self.tableView_results_observations.resizeColumnsToContents()
        resize_table_view_columns(table_view=self.tableView_results_observations,
                                  n=100,
                                  add_pixel=10)
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
            self.results_station_model.layoutChanged.connect(self.results_station_model_layoutChanged)

    def results_station_model_layoutChanged(self):
        """Execute whenever the layout of the station results model changes."""
        self.update_comboBox_results_stations_statistics_select_col_based_on_table_view()
        self.display_station_results_statistics()

    def update_comboBox_results_stations_statistics_select_col_based_on_table_view(self):
        """Update items in the combobox for selecting the stations results table column for presenting statistics."""
        self.comboBox_results_stations_statistics_select_col.blockSignals(True)
        data_columns_dict = self.results_station_model.get_columns_for_descriptive_statistics()
        # Get current item:
        idx, current_column_name = self.get_selected_station_results_stats_column()
        self.comboBox_results_stations_statistics_select_col.clear()
        # Add items (text and data items):
        for col_name, short_description in data_columns_dict.items():
            self.comboBox_results_stations_statistics_select_col.addItem(short_description, userData=col_name)
        # Try to select the previous item again:
        if idx != -1:  # Previous selection available
            try:
                current_short_description = self.results_station_model.get_short_column_description(
                    current_column_name)
                self.comboBox_results_stations_statistics_select_col.setCurrentText(current_short_description)
            except:
                pass
        self.comboBox_results_stations_statistics_select_col.blockSignals(False)

    def get_selected_station_results_stats_column(self):
        """Get the selected column for displaying station results statistics in the GUI."""
        idx = self.comboBox_results_stations_statistics_select_col.currentIndex()
        column_name = self.comboBox_results_stations_statistics_select_col.currentData()
        return idx, column_name

    def display_station_results_statistics(self):
        """Display station results statistics in the GUI."""
        # Get current selection:
        idx, col_name = self.get_selected_station_results_stats_column()
        model_data_df = self.results_station_model.get_model_data_df()
        if (idx == -1) or (model_data_df is None):
            self.label_results_stations_statistics_mean.clear()
            self.label_results_stations_statistics_std.clear()
            self.label_results_stations_statistics_min.clear()
            self.label_results_stations_statistics_max.clear()
            self.label__results_stations_statistics_median.clear()
            self.label__results_stations_statistics_iqr.clear()
            return
        data_col =  model_data_df[col_name]
        iqr = data_col.quantile(0.75) - data_col.quantile(0.25)
        self.label_results_stations_statistics_mean.setText(f'{data_col.mean():.3f}')
        self.label_results_stations_statistics_std.setText(f'{data_col.std():.3f}')
        self.label_results_stations_statistics_min.setText(f'{data_col.min():.3f}')
        self.label_results_stations_statistics_max.setText(f'{data_col.max():.3f}')
        self.label__results_stations_statistics_median.setText(f'{data_col.median():.3f}')
        self.label__results_stations_statistics_iqr.setText(f'{iqr:.3f}')

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_results_station_table_view(self, lsm_run_index: int, station_name=None, survey_name=None):
        """Update the station results table view after changing the table model."""
        self.results_station_model.load_lsm_runs(self.campaign.lsm_runs)
        self.results_station_model.update_view_model(lsm_run_index,
                                                     station_name=station_name,
                                                     survey_name=survey_name,
                                                     surveys=self.campaign.surveys)
        self.results_station_model.layoutChanged.emit()  # Show changes in table view
        # self.tableView_results_stations.resizeColumnsToContents()
        resize_table_view_columns(table_view=self.tableView_results_stations,
                                  n=10,
                                  add_pixel=10)

    def set_up_results_correlation_matrix_view_model(self):
        """Set up the view model for the correlation matrix table view."""
        try:
            self.results_correlation_matrix_model = ResultsCorrelationMatrixModel(self.campaign.lsm_runs)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No LSM-adjustment results available!')
            self.statusBar().showMessage(f"No LSM-adjustment results available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.tableView_results_correlation_matrix.setModel(self.results_correlation_matrix_model)
            self.tableView_results_correlation_matrix.resizeColumnsToContents()

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_results_correlation_matrix_table_view(self, lsm_run_index: int, station_name=None, survey_name=None):
        """Update the correlation matrix table view after changing the table model.

        Notes
        -----
        `idx = -1` indicates that no valid lsm run was selected. In this case the table view should be cleared.
        """
        self.results_correlation_matrix_model.load_lsm_runs(self.campaign.lsm_runs)
        self.results_correlation_matrix_model.update_view_model(lsm_run_index,
                                                                station_name=station_name,
                                                                survey_name=survey_name)
        self.results_correlation_matrix_model.layoutChanged.emit()  # Show changes in table view
        # self.tableView_results_correlation_matrix.resizeColumnsToContents()  # Extremely slow!
        resize_table_view_columns(table_view=self.tableView_results_correlation_matrix,
                                  n=1,
                                  add_pixel=4)

    @pyqtSlot(int)
    def on_comboBox_results_lsm_run_selection_current_index_changed(self, index: int):
        """Invoked whenever the index of the selected item in the combobox changed."""
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

    @conditional_decorator(time_it, settings.DEBUG_TIME_IT)
    def update_results_tab(self, select_latest_item=False):
        """Update the results tab of the main tab widget with the content of the selected lsm run."""
        # Update combo box content for lsm run selection:
        self.update_comboBox_lsm_run_selection(select_latest_item=select_latest_item)

        # Get the currently selected lsm run object:
        idx, time_str = self.get_selected_lsm_run()
        if idx != -1:  # Valid index
            lsm_run = self.campaign.lsm_runs[idx]

            # Update the comment string in the line edit below the lsm run selection combo box:
            self.label_results_lsm_run_comment_display.setText(lsm_run.comment)

            # Update the list in the combobox:
            self.update_comboBox_results_selection_station(observed_stations=lsm_run.observed_stations)
            self.update_comboBox_results_selection_surrvey(survey_names=list(lsm_run.setups.keys()))

            # Get data from LSM object an populate GUI widgets:
            # - Info tab:
            self.label_results_comment.setText(lsm_run.comment)
            self.label_results_adjustment_method.setText(settings.ADJUSTMENT_METHODS[lsm_run.lsm_method])
            self.label_results_time_and_date.setText(lsm_run.init_time.strftime("%Y-%m-%d, %H:%M:%S"))
            try:
                number_of_iterations = lsm_run.number_of_iterations
                if number_of_iterations == 0:
                    number_of_iterations_str = 'No iteration'
                else:
                    number_of_iterations_str = str(number_of_iterations)
            except:
                number_of_iterations_str = 'No iteration'
            self.label_results_number_of_iterations.setText(number_of_iterations_str)
            if lsm_run.s02_a_posteriori is not None:
                self.label_results_sig0.setText(f'{lsm_run.s02_a_posteriori:1.3f}')
            else:
                self.label_results_sig0.clear()
            if lsm_run.write_log:
                self.plainTextEdit_results_log.setPlainText(lsm_run.log_str)

            try:
                if lsm_run.number_of_outliers is not None:
                    self.label_results_number_of_outliers.setText(str(lsm_run.number_of_outliers))
                else:
                    self.label_results_number_of_outliers.clear()
            except:
                self.label_results_number_of_outliers.clear()
            try:
                self.label_results_goodness_of_fit_test_status.setText(lsm_run.global_model_test_status)
            except:
                self.label_results_goodness_of_fit_test_status.clear()

            # Get station and/or survey names for filtering the displayed data:
            stat_idx, current_station_name = self.get_selected_station()
            if current_station_name == 'All stations':
                station_name = None
            else:
                station_name = current_station_name
            survey_idx, current_survey_name = self.get_selected_survey_results()
            if current_survey_name == 'All surveys':
                survey_name = None
            else:
                survey_name = current_survey_name

            # Get unique station colors for plotting data:
            self.station_colors_dict_results = get_station_color_dict(lsm_run.stat_obs_df['station_name'].to_list(),
                                                                      randomize=True)

            # Update widgets:
            self.update_results_station_table_view(idx, station_name=station_name, survey_name=survey_name)
            self.update_results_correlation_matrix_table_view(idx)
            self.update_results_observation_table_view(idx, station_name=station_name, survey_name=survey_name)
            self.update_results_drift_table_view(idx, survey_name=survey_name)
            self.update_results_vg_table_view(idx)
            self.update_results_obs_plots()
            self.update_drift_plot()
            self.update_vg_plot()
        else:  # invalid index => Reset results views
            self.station_colors_dict_results = {}
            self.label_results_comment.clear()
            self.label_results_adjustment_method.clear()
            self.label_results_number_of_iterations.clear()
            self.label_results_time_and_date.clear()
            self.label_results_sig0.clear()
            self.label_results_lsm_run_comment_display.clear()
            self.label_results_goodness_of_fit_test_status.clear()
            self.label_results_number_of_outliers.clear()
            self.plainTextEdit_results_log.clear()
            self.update_results_station_table_view(idx, station_name=None, survey_name=None)  # Can handle idx=-1
            self.update_results_correlation_matrix_table_view(idx)  # Can handle idx=-1
            self.update_results_observation_table_view(idx, station_name=None, survey_name=None)  # Can handle idx=-1
            self.update_results_drift_table_view(idx, survey_name=None)
            self.update_results_vg_table_view(idx)
            self.update_comboBox_results_selection_station(observed_stations=[])
            self.update_comboBox_results_selection_surrvey(survey_names=[])
            self.update_results_obs_plots()
            self.update_drift_plot()
            self.update_vg_plot()

    def update_comboBox_results_obs_plot_select_data_column_based_on_table_view(self):
        """Update the observaterion results data column selection combo box in the results tab."""
        self.comboBox_results_obs_plot_select_data_column.blockSignals(True)
        # Get data columns with data that is plottable from the observations results table view model:
        data_columns_dict = self.results_observation_model.get_plotable_columns()
        data_columns = list(data_columns_dict.keys())
        # Get current item:
        idx, current_column_name = self.get_selected_obs_data_column()
        self.comboBox_results_obs_plot_select_data_column.clear()
        # Add items (text and data items):
        for col_name in data_columns:
            # data_columns_short_description.append(self.results_observation_model.get_short_column_description(name))
            short_description = self.results_observation_model.get_short_column_description(col_name)
            self.comboBox_results_obs_plot_select_data_column.addItem(short_description, userData=col_name)

        # Try to select the previous item again:
        if idx != -1:  # Previous selection available
            try:
                current_short_description = self.results_observation_model.get_short_column_description(
                    current_column_name)
                self.comboBox_results_obs_plot_select_data_column.setCurrentText(current_short_description)
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
        """Update the survey selection combobox in the results tab, based on the current lsm run."""
        self.comboBox_results_selection_survey.blockSignals(True)
        # Get current item:
        survey_idx, current_survey_name = self.get_selected_survey_results()
        self.comboBox_results_selection_survey.clear()
        self.comboBox_results_selection_survey.addItems(['All surveys'] + survey_names)
        # Try to select the previous item again:
        if survey_idx != -1:  # Previous selection available
            self.comboBox_results_selection_survey.setCurrentText(current_survey_name)
        self.comboBox_results_selection_survey.blockSignals(False)

    def update_comboBox_stations_selection_surrvey(self, survey_names: list):
        """Update the survey selection combobox in the stations tab, based on the current surveys in the campaign."""
        self.comboBox_stations_plot_obs_map_surveys.blockSignals(True)
        # Get current item:
        survey_idx, current_survey_name = self.get_selected_survey_stations()
        self.comboBox_stations_plot_obs_map_surveys.clear()
        self.comboBox_stations_plot_obs_map_surveys.addItems(['All surveys'] + survey_names)
        # Try to select the previous item again:
        if survey_idx != -1:  # Previous selection available
            self.comboBox_stations_plot_obs_map_surveys.setCurrentText(current_survey_name)
        self.comboBox_stations_plot_obs_map_surveys.blockSignals(False)

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
            time_str = self.campaign.lsm_runs[idx].init_time.strftime("%Y-%m-%d, %H:%M:%S")
            comment_str = self.campaign.lsm_runs[idx].comment
            msg_text = f'Date and time: {time_str}\nComment: {comment_str}'
            reply = QMessageBox.question(self,
                                         'Delete LSM run?',
                                         msg_text,
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                try:
                    self.campaign.delete_lsm_run(idx)
                except Exception as e:
                    QMessageBox.critical(self, 'Error!', str(e))
                    self.statusBar().showMessage(f'LSM run "{time_str}" not deleted.')
                else:
                    self.statusBar().showMessage(f'LSM run "{time_str}" deleted.')
                    self.update_results_tab()
            else:
                pass

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

    def get_selected_survey_results(self):
        """Get the selected survey in the results tab."""
        survey_name = self.comboBox_results_selection_survey.currentText()
        idx = self.comboBox_results_selection_survey.currentIndex()
        return idx, survey_name

    def get_selected_survey_stations(self):
        """Get the selected survey in the stations tab for displaying observations on the map."""
        survey_name = self.comboBox_stations_plot_obs_map_surveys.currentText()
        idx = self.comboBox_stations_plot_obs_map_surveys.currentIndex()
        return idx, survey_name

    def get_selected_obs_data_column(self):
        """Get the selected observation results data column (dataframe column name)."""
        # short_description = self.comboBox_results_obs_plot_select_data_column.currentText()
        idx = self.comboBox_results_obs_plot_select_data_column.currentIndex()
        # Convert comboBox text to dataframe column name:
        # column_name = self.results_observation_model.get_column_name_from_short_description(short_description)
        column_name = self.comboBox_results_obs_plot_select_data_column.currentData()
        return idx, column_name

    def on_checkBox_obs_plot_setup_data_state_changed(self):
        """Invoke, whenever the state of the checkbox changes."""
        self.refresh_observation_plot()

    def refresh_observation_plot(self):
        """Refresh the observation plot."""
        survey_name, setup_id = self.get_obs_tree_widget_selected_item()
        self.plot_observations(survey_name)

    def on_pushbutton_obs_comp_setup_data(self):
        """Invoked when pushing the button 'pushbutton_obs_comp_setup_data'."""
        self.compute_setup_data_for_campaign()
        survey_name, setup_id = self.get_obs_tree_widget_selected_item()
        self.observation_model.update_view_model(survey_name,
                                                 setup_id,
                                                 gui_simple_mode=self.dlg_options.gui_simple_mode)
        self.refresh_observation_plot()
        self.update_setup_table_view(survey_name, setup_id)

    def update_setup_table_view(self, survey_name, setup_id):
        """Update the setups table view after changing the table model."""
        self.setup_model.update_view_model(survey_name,
                                           setup_id,
                                           gui_simple_mode=self.dlg_options.gui_simple_mode)
        self.setup_model.layoutChanged.emit()  # Show changes in table view
        self.tableView_observations_setups.resizeColumnsToContents()
        # Show additional infos on the currently visualized setup data in the GUI:
        try:
            if self.setup_model.get_ref_heigth_type in settings.REFERENCE_HEIGHT_TYPE.keys():
                self.label_obs_setups_ref_height.setText(f'{self.setup_model.get_ref_heigth_type} ({settings.REFERENCE_HEIGHT_TYPE[self.setup_model.get_ref_heigth_type]})')
            else:
                self.label_obs_setups_ref_height.setText(f'{self.setup_model.get_ref_heigth_type}')
            if self.setup_model.get_tidal_corr_type in settings.TIDE_CORRECTION_TYPES.keys():
                self.label_obs_setups_tidal_corr.setText(f'{self.setup_model.get_tidal_corr_type} ({settings.TIDE_CORRECTION_TYPES[self.setup_model.get_tidal_corr_type]})')
            else:
                self.label_obs_setups_tidal_corr.setText(f'{self.setup_model.get_tidal_corr_type}')
            if self.setup_model.get_atm_pres_corr_type in settings.ATM_PRES_CORRECTION_TYPES.keys():
                self.label_obs_setups_atm_pres_corr.setText(f'{self.setup_model.get_atm_pres_corr_type} ({settings.ATM_PRES_CORRECTION_TYPES[self.setup_model.get_atm_pres_corr_type]})')
            else:
                self.label_obs_setups_atm_pres_corr.setText(f'{self.setup_model.get_atm_pres_corr_type}')

        except KeyError:
            self.label_obs_setups_ref_height.setText('')
            self.label_obs_setups_tidal_corr.setText('')
            self.label_obs_setups_atm_pres_corr.setText('')

    def compute_setup_data_for_campaign(self):
        """Compute setup data for the campaign."""

        # Get GUI settings:
        active_obs_only_for_ref_epoch = self.dlg_setup_data.checkBox_drift_ref_epoch_active_obs_only.checkState() == Qt.Checked
        if self.dlg_setup_data.radioButton_variance_weighted_mean.isChecked():
            method = 'variance_weighted_mean'
        elif self.dlg_setup_data.radioButton_individual_obs.isChecked():
            method = 'individual_obs'
        else:
            method = 'variance_weighted_mean' # Use default
        if self.dlg_setup_data.radioButton_sd_from_obsfile.isChecked():
            method_sd = 'sd_from_obs_file'
        elif self.dlg_setup_data.radioButton_sd_default_per_obs.isChecked():
            method_sd = 'sd_default_per_obs'
        elif self.dlg_setup_data.radioButton_sd_defaul_per_setup.isChecked():
            method_sd = 'sd_default_per_setup'
        else:
            method_sd = 'sd_from_obs_file'
        default_sd_mugal = self.dlg_setup_data.spinBox_sd_default.value()

        try:
            self.campaign.calculate_setup_data(obs_type='reduced',
                                               active_obs_only_for_ref_epoch=active_obs_only_for_ref_epoch,
                                               method=method,
                                               method_sd=method_sd,
                                               default_sd_mugal=default_sd_mugal,
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
        num_of_lsm_runs_in_camapaign_before_adjustment = len(self.campaign.lsm_runs)
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
            add_const_to_sd_mugal = self.dlg_estimation_settings.doubleSpinBox_add_const_sd.value()
            scaling_factor_obs_sd = self.dlg_estimation_settings.doubleSpinBox_mult_factor_sd.value()
            self.dlg_estimation_settings.radioButton_drift_ref_epoch_first_ob_campaign.isChecked()
            if self.dlg_estimation_settings.radioButton_drift_ref_epoch_first_obs_survey.isChecked() and not self.dlg_estimation_settings.radioButton_drift_ref_epoch_first_ob_campaign.isChecked():
                drift_ref_epoch_type = 'survey'
            elif not self.dlg_estimation_settings.radioButton_drift_ref_epoch_first_obs_survey.isChecked() and self.dlg_estimation_settings.radioButton_drift_ref_epoch_first_ob_campaign.isChecked():
                drift_ref_epoch_type = 'campaign'
            else:
                raise AssertionError('Invalid definition of the drift reference epoch in the GUI!')

            # Autoscale settings:
            autoscale_s0_a_posteriori = self.dlg_estimation_settings.checkBox_iterative_s0_scaling.checkState() == Qt.Checked

            iteration_approach = self.dlg_estimation_settings.comboBox_iteration_approach.currentText()
            s02_target = self.dlg_estimation_settings.doubleSpinBox_target_s02.value()
            s02_target_delta = self.dlg_estimation_settings.doubleSpinBox_delta_target_s02.value()
            max_number_iterations = self.dlg_estimation_settings.spinBox_max_number_of_iterations.value()
            add_const_to_sd_of_observations_step_size_mugal = self.dlg_estimation_settings.doubleSpinBox_initial_step_size.value()
            max_total_additive_const_to_sd_mugal = self.dlg_estimation_settings.doubleSpinBox_max_additive_const_to_sd.value()
            max_multiplicative_factor_to_sd_percent = self.dlg_estimation_settings.doubleSpinBox_max_multiplicative_factor_to_sd_percent.value()
            min_multiplicative_factor_to_sd_percent = self.dlg_estimation_settings.doubleSpinBox_min_multiplicative_factor_to_sd_percent.value()
            initial_step_size_percent = self.dlg_estimation_settings.doubleSpinBox_initial_step_size_percent.value()
            noise_floor_mugal = self.dlg_estimation_settings.doubleSpinBox_gravity_noise_floor_mugal.value()

            # Initialize LSM object and add it to the campaign object:
            self.campaign.initialize_and_add_lsm_run(lsm_method=lsm_method, comment=comment, write_log=True)

            # Run the estimation:
            if lsm_method == 'LSM_diff':
                if autoscale_s0_a_posteriori:
                    self.campaign.lsm_runs[-1].adjust_autoscale_s0(
                        iteration_approach=iteration_approach,
                        s02_target=s02_target,
                        s02_target_delta=s02_target_delta,
                        max_number_iterations=max_number_iterations,
                        add_const_to_sd_of_observations_step_size_mugal=add_const_to_sd_of_observations_step_size_mugal,
                        max_total_additive_const_to_sd_mugal=max_total_additive_const_to_sd_mugal,
                        multiplicative_factor_step_size_percent=initial_step_size_percent,
                        max_multiplicative_factor_to_sd_percent=max_multiplicative_factor_to_sd_percent,
                        min_multiplicative_factor_to_sd_percent=min_multiplicative_factor_to_sd_percent,
                        drift_pol_degree=degree_drift_polynomial,
                        sig0_mugal=sig0,
                        scaling_factor_datum_observations=weight_factor_datum,
                        add_const_to_sd_of_observations_mugal=add_const_to_sd_mugal,
                        scaling_factor_for_sd_of_observations=scaling_factor_obs_sd,
                        confidence_level_chi_test=confidence_level_chi_test,
                        confidence_level_tau_test=confidence_level_tau_test,
                        drift_ref_epoch_type=drift_ref_epoch_type,
                        noise_floor_mugal=noise_floor_mugal,
                        verbose=IS_VERBOSE,
                    )
                else:  # no autoscale
                    self.campaign.lsm_runs[-1].adjust(
                        drift_pol_degree=degree_drift_polynomial,
                        sig0_mugal=sig0,
                        scaling_factor_datum_observations=weight_factor_datum,
                        add_const_to_sd_of_observations_mugal=add_const_to_sd_mugal,
                        scaling_factor_for_sd_of_observations=scaling_factor_obs_sd,
                        confidence_level_chi_test=confidence_level_chi_test,
                        confidence_level_tau_test=confidence_level_tau_test,
                        drift_ref_epoch_type=drift_ref_epoch_type,
                        noise_floor_mugal=noise_floor_mugal,
                        verbose=IS_VERBOSE
                    )
            elif lsm_method == 'LSM_non_diff':
                if autoscale_s0_a_posteriori:
                    self.campaign.lsm_runs[-1].adjust_autoscale_s0(
                        iteration_approach=iteration_approach,
                        s02_target=s02_target,
                        s02_target_delta=s02_target_delta,
                        max_number_iterations=max_number_iterations,
                        add_const_to_sd_of_observations_step_size_mugal=add_const_to_sd_of_observations_step_size_mugal,
                        max_total_additive_const_to_sd_mugal=max_total_additive_const_to_sd_mugal,
                        multiplicative_factor_step_size_percent=initial_step_size_percent,
                        max_multiplicative_factor_to_sd_percent=max_multiplicative_factor_to_sd_percent,
                        min_multiplicative_factor_to_sd_percent=min_multiplicative_factor_to_sd_percent,
                        drift_pol_degree=degree_drift_polynomial,
                        sig0_mugal=sig0,
                        scaling_factor_datum_observations=weight_factor_datum,
                        add_const_to_sd_of_observations_mugal=add_const_to_sd_mugal,
                        scaling_factor_for_sd_of_observations=scaling_factor_obs_sd,
                        confidence_level_chi_test=confidence_level_chi_test,
                        confidence_level_tau_test=confidence_level_tau_test,
                        drift_ref_epoch_type=drift_ref_epoch_type,
                        noise_floor_mugal=noise_floor_mugal,
                        verbose=IS_VERBOSE,
                    )
                else:
                    self.campaign.lsm_runs[-1].adjust(drift_pol_degree=degree_drift_polynomial,
                                                      sig0_mugal=sig0,
                                                      scaling_factor_datum_observations=weight_factor_datum,
                                                      add_const_to_sd_of_observations_mugal=add_const_to_sd_mugal,
                                                      scaling_factor_for_sd_of_observations=scaling_factor_obs_sd,
                                                      confidence_level_chi_test=confidence_level_chi_test,
                                                      confidence_level_tau_test=confidence_level_tau_test,
                                                      drift_ref_epoch_type=drift_ref_epoch_type,
                                                      noise_floor_mugal=noise_floor_mugal,
                                                      verbose=IS_VERBOSE)
            elif lsm_method == 'MLR_BEV':
                self.campaign.lsm_runs[-1].adjust(drift_pol_degree=degree_drift_polynomial,
                                                  verbose=IS_VERBOSE)
            elif lsm_method == 'VG_LSM_nondiff':
                # Get all required parameters from the GUI specific for the VG estimation:
                vg_polynomial_ref_height_offset_m = self.dlg_estimation_settings.doubleSpinBox_vg_polynomial_ref_height_offset_m.value()
                vg_polynomial_degree = self.dlg_estimation_settings.spinBox_vg_polynomial_degree.value()
                self.campaign.lsm_runs[-1].adjust(drift_pol_degree=degree_drift_polynomial,
                                                  vg_polynomial_degree=vg_polynomial_degree,
                                                  vg_polynomial_ref_height_offset_m=vg_polynomial_ref_height_offset_m,
                                                  sig0_mugal=sig0,
                                                  confidence_level_chi_test=confidence_level_chi_test,
                                                  confidence_level_tau_test=confidence_level_tau_test,
                                                  verbose=IS_VERBOSE)

        except AssertionError as e:
            QMessageBox.critical(self, 'Error!', str(e))
            # Delete failed lsm run object
            self.campaign.lsm_runs = self.campaign.lsm_runs[0:num_of_lsm_runs_in_camapaign_before_adjustment]
            self.statusBar().showMessage(f"Error! No parameters estimated.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
            # Delete failed lsm run object
            self.campaign.lsm_runs = self.campaign.lsm_runs[0:num_of_lsm_runs_in_camapaign_before_adjustment]
            self.statusBar().showMessage(f"Error! No parameters estimated.")
        else:
            # No errors when computing the setup data:
            self.statusBar().showMessage(f"Parameters estimated successfully!")
            # Update list of lsm runs in results tab of the GUI:
            self.update_results_tab(select_latest_item=True)
        pass

    def on_apply_autoselection(self):
        """Apply autoselection on the currently selected setup or survey according to the predefined settings."""

        # Get autoselect parameters from settings dialog
        flag_apply_tilt = self.dlg_autoselect_settings.checkBox_tilt.isChecked()
        flag_apply_g_sd = self.dlg_autoselect_settings.checkBox_sd.isChecked()
        flag_apply_delta_g = self.dlg_autoselect_settings.checkBox_delta_g.isChecked()
        flag_apply_duration = self.dlg_autoselect_settings.checkBox_duration.isChecked()

        treshold_g_sd_mugal = self.dlg_autoselect_settings.spinBox_sd.value()
        treshold_tilt_arcsec = self.dlg_autoselect_settings.spinBox_tilt.value()
        treshold_delta_sd_mugal = self.dlg_autoselect_settings.spinBox_delta_g.value()
        delta_g_number_of_points = self.dlg_autoselect_settings.spinBox_n.value()
        treshold_duration_sec = self.dlg_autoselect_settings.spinBox_duration.value()

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

        if flag_apply_duration:
            surv.autselect_duration(threshold_sec=treshold_duration_sec, setup_id=setup_id, verbose=IS_VERBOSE)
        if flag_apply_tilt:
            surv.autselect_tilt(threshold_arcsec=treshold_tilt_arcsec, setup_id=setup_id, verbose=IS_VERBOSE)
        if flag_apply_g_sd:
            surv.autselect_g_sd(threshold_mugal=treshold_g_sd_mugal, obs_type=reference_data, setup_id=setup_id,
                                verbose=IS_VERBOSE)
        if flag_apply_delta_g:
            surv.autselect_delta_g(threshold_mugal=treshold_delta_sd_mugal, n_obs=delta_g_number_of_points,
                                   obs_type=reference_data, setup_id=setup_id, verbose=IS_VERBOSE)

        # Update data visualization in GUI
        self.update_obs_table_view(survey_name, setup_id)
        self.plot_observations(survey_name)
        self.update_obs_tree_widgget_from_observation_model()

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
            # self.update_obs_table_view(survey_name, setup_id)  # TODO: Why here? Needed?
            # self.plot_observations(survey_name)  # TODO: Why here? Needed?
        else:
            if IS_VERBOSE:
                print('No item or multiple items selected!')
            survey_name = None
            setup_id = None
        return survey_name, setup_id

    def on_menu_observations_autoselection_settings(self):
        """Launch dialog for defining the autoselection settings."""
        return_value = self.dlg_autoselect_settings.exec()

    def on_menu_observations_setup_data(self):
        """Launch dialog for defining setup data options."""
        return_value = self.dlg_setup_data.exec()

    def on_menu_estimation_settings(self):
        """Launch dialog for defining the estimation settings."""
        return_value = self.dlg_estimation_settings.exec()

    def on_menu_gis_export_settings(self):
        """Launch dialog for defining gis export settings."""
        return_value = self.dlg_gis_export_settings.exec()

    @pyqtSlot()
    def on_menu_file_export_results(self):
        """Launch dialog for exporting results of an LSM run."""
        _OBS_FILE_EXPORT_TYPES = {'all observations': 'all_obs',
                                  'only active observations': 'active_only',
                                  'only inactive observations': 'inactive_only'}
        dlg = DialogExportResults(campaign=self.campaign)
        return_value = dlg.exec()
        flag_export_successful = True
        if return_value == QDialog.Accepted:
            # Get GUI settings and data:
            append_lsm_run_comment_to_filenames = False
            output_path = dlg.label_export_path_show.text()
            filename = self.campaign.campaign_name

            # Check if the output path exists:
            if not os.path.exists(output_path):
                QMessageBox.critical(self, 'Error!',
                                     f"The campaign's output path doe not exist ({output_path})! Change the path.")
                flag_export_successful = False
            else:
                # Get LSM run index and check data availability:
                lsm_run_idx = dlg.comboBox_select_lsm_run.currentIndex() - 1
                flag_lsm_run_selected = dlg.flag_lsm_runs_available and (lsm_run_idx > -1)

                try:
                    # Append lsm run comment to filename, if available:
                    if flag_lsm_run_selected:
                        lsm_run = self.campaign.lsm_runs[lsm_run_idx]
                        # Append LSM comment to filename?
                        if dlg.checkBox_add_lsm_comment_to_filename.checkState() == Qt.Checked:
                            append_lsm_run_comment_to_filenames = True
                            if not lsm_run.comment:  # Empty string
                                append_lsm_run_comment_to_filenames = False
                        else:
                            append_lsm_run_comment_to_filenames = False
                        if append_lsm_run_comment_to_filenames:
                            filename = self.campaign.campaign_name + '_' + lsm_run.comment

                except Exception as e:
                    QMessageBox.critical(self, 'Error!', str(e))
                    flag_export_successful = False

                # Write observation list (CSV file):
                if dlg.checkBox_write_observation_list.checkState() == Qt.Checked:
                    try:
                        filename_obs_list = filename + '_obs.csv'
                        export_type = _OBS_FILE_EXPORT_TYPES[
                            dlg.comboBox_observation_list_export_options.currentText()]
                        if flag_lsm_run_selected:
                            self.campaign.write_obs_list_of_lsm_run_csv(
                                filename_csv=os.path.join(output_path, filename_obs_list),
                                lsm_run_index=lsm_run_idx,
                                export_type=export_type,
                                verbose=IS_VERBOSE)
                        else:
                            # if no lsm run selected:
                            if dlg.flag_observation_data_available:
                                self.campaign.write_obs_list_csv(
                                    filename_csv=os.path.join(output_path, filename_obs_list),
                                    export_type=export_type,
                                    verbose=IS_VERBOSE)
                    except Exception as e:
                        QMessageBox.critical(self, 'Error!', str(e))
                        flag_export_successful = False

                # If LSM run selected:
                if flag_lsm_run_selected:
                    # Write nsb file:
                    if lsm_run.lsm_method in settings.LSM_METHODS_NSD_FILE_EXPORT:
                        try:
                            if dlg.checkBox_write_nsb_file.checkState() == Qt.Checked:
                                filename_nsb = filename + '.nsb'
                                if dlg.radioButton_mean_dhb_dhf.isChecked():
                                    vertical_offset_mode = 'mean'
                                elif dlg.radioButton_first_dhb_dhf.isChecked():
                                    vertical_offset_mode = 'first'
                                else:
                                    raise AssertionError(f'Undefined vertical offset mode!')
                                if dlg.radioButton_export_se.isChecked():
                                    formal_error_type = 'se'
                                    if IS_VERBOSE:
                                        print(f' - Write gravity standard errors (SE) to the nsb file.')
                                elif dlg.radioButton_export_sd.isChecked():
                                    if IS_VERBOSE:
                                        print(f' - Write gravity standard deviations (SD) to the nsb file.')
                                    formal_error_type = 'sd'
                                else:
                                    raise AssertionError(f'Invalid formal error type! Valid: "sd" or "se".')

                                if dlg.checkBox_nsb_remove_datum_stations.checkState() == Qt.Checked:
                                    exclude_datum_stations = True
                                else:
                                    exclude_datum_stations = False
                                self.campaign.write_nsb_file(filename=os.path.join(output_path, filename_nsb),
                                                             lsm_run_index=lsm_run_idx,
                                                             vertical_offset_mode=vertical_offset_mode,
                                                             exclude_datum_stations=exclude_datum_stations,
                                                             formal_error_type=formal_error_type,
                                                             verbose=IS_VERBOSE)
                        except Exception as e:
                            QMessageBox.critical(self, 'Error!', str(e))
                            flag_export_successful = False

                    # Write log file:
                    try:
                        if dlg.checkBox_write_log_file.checkState() == Qt.Checked:
                            filename_log = filename + '.log'
                            self.campaign.write_log_file(filename=os.path.join(output_path, filename_log),
                                                         lsm_run_index=lsm_run_idx,
                                                         verbose=IS_VERBOSE)
                    except Exception as e:
                        QMessageBox.critical(self, 'Error!', str(e))
                        flag_export_successful = False

                    # Save drift plot to PNG file:
                    # Reference: https://pyqtgraph.readthedocs.io/en/latest/exporting.html
                    # - WARNING: Sometimes a linAlg error occurrs ('MAtrix singular') when trying to export the
                    #            drift plot. It is not possible at this point in the code to cath this error! Therefore,
                    #            the export of the drift plot is the last step in the file export procedure! This was
                    #            at leat all other files can be written.
                    try:
                        if dlg.checkBox_save_drift_plot_png.checkState() == Qt.Checked:
                            # To prevent an error the plot has to be shown first by bringing the drift plot tab to front:
                            # - Main tab widget:
                            for idx in range(self.tabWidget_Main.count()):
                                if self.tabWidget_Main.tabText(idx) == "Results":
                                    self.tabWidget_Main.setCurrentIndex(idx)
                            # - Results tab widget:
                            for idx in range(self.tabWidget_results.count()):
                                if self.tabWidget_results.tabText(idx) == "Drift Plot":
                                    self.tabWidget_results.setCurrentIndex(idx)
                            filename_png = filename + '_drift_plot.png'
                            # The linAlg error is raised in the following line of code:
                            exporter = pg.exporters.ImageExporter(self.graphicsLayoutWidget_results_drift_plot.scene())
                            flag_export_successful = exporter.export(os.path.join(output_path, filename_png))
                    except Exception as e:
                        QMessageBox.critical(self, 'Error!', str(e))
                        flag_export_successful = False

                    # Save VG plot to PNG file:
                    if lsm_run.lsm_method in settings.LSM_METHODS_VG_PLOT:
                        try:
                            if dlg.checkBox_save_vg_plot_png.checkState() == Qt.Checked:
                                # To prevent an error the plot has to be shown first by bringing the drift plot tab to front:
                                # - Main tab widget:
                                for idx in range(self.tabWidget_Main.count()):
                                    if self.tabWidget_Main.tabText(idx) == "Results":
                                        self.tabWidget_Main.setCurrentIndex(idx)
                                # - Results tab widget:
                                for idx in range(self.tabWidget_results.count()):
                                    if self.tabWidget_results.tabText(idx) == "VG Plot":
                                        self.tabWidget_results.setCurrentIndex(idx)
                                filename_png = filename + '_vg_plot.png'
                                # The linAlg error is raised in the following line of code:
                                exporter = pg.exporters.ImageExporter(self.graphicsLayoutWidget_results_drift_plot.scene())
                                flag_export_successful = exporter.export(os.path.join(output_path, filename_png))
                        except Exception as e:
                            QMessageBox.critical(self, 'Error!', str(e))
                            flag_export_successful = False

                    # Save shapefile:
                    if lsm_run.lsm_method in settings.LSM_METHODS_GIS_EXPORT:
                        if dlg.checkBox_gis_write_shapefile.isChecked():
                            if self.dlg_gis_export_settings.radioButton_campaign_output_dir.isChecked():
                                gis_output_dir = self.campaign.output_directory
                            else:
                                gis_output_dir = self.dlg_gis_export_settings.lineEdit_gis_output_dir.text()
                            if not os.path.isdir(gis_output_dir):
                                QMessageBox.critical(self, 'Error!',
                                                     f'Invalid output directory for GIS files: {gis_output_dir}')
                            else:
                                try:
                                    epsg_code = int(self.dlg_gis_export_settings.lineEdit_stat_coord_epsg.text())
                                except ValueError:
                                    QMessageBox.critical(self, 'Error!', 'Invalid EPSG code. Need to be an integer value.')
                                else:
                                    try:
                                        filename_shp = os.path.join(gis_output_dir, filename + '_obs.shp')
                                        lsm_run.export_obs_results_shapefile(filename=filename_shp, epsg_code=epsg_code)
                                    except Exception as e:
                                        QMessageBox.critical(self, 'Error!', str(e))
                                        flag_export_successful = False
                                    try:
                                        filename_shp = os.path.join(gis_output_dir, filename + '_stat.shp')
                                        lsm_run.export_obs_results_shapefile(filename=filename_shp, epsg_code=epsg_code)
                                    except Exception as e:
                                        QMessageBox.critical(self, 'Error!', str(e))
                                        flag_export_successful = False

            if flag_export_successful:
                self.statusBar().showMessage(f"Export to {output_path} successful!")
            else:
                self.statusBar().showMessage(f"Problems at data export!")
        else:
            self.statusBar().showMessage(f"No exports.")

    def set_up_obseration_plots_widget(self):
        """Set up `self.GraphicsLayoutWidget_observations`."""
        l = self.GraphicsLayoutWidget_observations
        l.setBackground('w')  # white background color

        # date_axis = TimeAxisItem(orientation='bottom')

        # Create sub-plots:
        # Gravity g [µGal]
        self.plot_obs_g = l.addPlot(0, 0, name='plot_obs_g', axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_g.setLabel(axis='left', text='g [µGal]')
        self.plot_obs_g.addLegend()

        # Standard deviation of gravity g [µGal]
        self.plot_obs_sd_g = l.addPlot(1, 0, name='plot_obs_sd_g',
                                       axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_sd_g.setLabel(axis='left', text='SD [µGal]')
        self.plot_obs_sd_g.addLegend()
        self.plot_obs_sd_g.setXLink(self.plot_obs_g)

        # Instrument tilt in X and Y directions [arcsec]
        self.plot_obs_tilt = l.addPlot(2, 0, name='plot_obs_tilt',
                                       axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_tilt.setLabel(axis='left', text='tilt [asec]')
        self.plot_obs_tilt.addLegend()
        self.plot_obs_tilt.setXLink(self.plot_obs_g)

        # Observation corrections [µGal]
        self.plot_obs_corrections = l.addPlot(3, 0, name='plot_obs_corrections',
                                              axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        self.plot_obs_corrections.setLabel(axis='left', text='Corrections [µGal]')
        self.plot_obs_corrections.addLegend()
        self.plot_obs_corrections.setXLink(self.plot_obs_g)

    def plot_observations(self, survey_name=None):
        """Plots observation data to the GraphicsLayoutWidget.

        Parameters
        ----------
        survey_name : str, optional (default=None)
            Specifies the survey for which the data should be plotted. `None` indicates, that the plots should be just
            cleared (e.g. new campaign and no observation data available).
        """
        obs_df = self.observation_model.get_data
        # Wipe plots, if survey_name is None:
        self.plot_obs_g.clear()
        self.plot_obs_sd_g.clear()
        self.plot_obs_corrections.clear()
        self.plot_obs_tilt.clear()

        if obs_df is not None and survey_name is not None:
            setup_df = self.observation_model.get_setup_data
            obs_epoch_timestamps = (obs_df['obs_epoch'].values - np.datetime64(
                '1970-01-01T00:00:00')) / np.timedelta64(1,
                                                          's')
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
                g_mugal = obs_df['g_red_mugal'].astype(float).values
                sd_g_mugal = obs_df['sd_g_red_mugal'].values
                corr_tide = obs_df['corr_tide_red_mugal'].values
                corr_tide_name = self.campaign.surveys[survey_name].red_tide_correction_type

                # try-except scheme to provide downward compatibility for versions < 0.2.3
                try:
                    corr_atm_pres_name = self.campaign.surveys[survey_name].red_atm_pres_correction_type
                except AttributeError:
                    corr_atm_pres_name = '' # ... won't be plotted in GUI.
                corr_atm_pres = obs_df['corr_atm_pres_red_mugal'].values
                if (corr_atm_pres == None).all():
                    corr_atm_pres_name = ''  # ... won't be plotted in GUI.
                elif np.isnan(corr_atm_pres).all():
                    corr_atm_pres_name = ''  # ... won't be plotted in GUI.
                try:
                    corr_description_str = self.campaign.surveys[survey_name].red_tide_correction_description
                except AttributeError:
                    pass
                else:
                    if corr_description_str:
                        corr_tide_name = corr_tide_name + ': ' + corr_description_str
                ref_height_name = self.campaign.surveys[survey_name].red_reference_height_type
            else:
                g_mugal = obs_df['g_obs_mugal'].values
                sd_g_mugal = obs_df['sd_g_obs_mugal'].values
                corr_tide = obs_df['corr_tide_mugal'].values
                corr_tide_name = self.campaign.surveys[survey_name].obs_tide_correction_type
                ref_height_name = self.campaign.surveys[survey_name].obs_reference_height_type
                corr_atm_pres_name = ''  # ... won't be plotted in GUI.

            # Gravity g [µGal]
            # - Plot with marker symbols according to their 'keep_obs' states and connect the 'sigPointsClicked' event.
            self.plot_obs_g.clear()
            self.plot_obs_g.setLabel(axis='left', text=f'g [µGal]')
            self.plot_obs_g.setTitle(f'Observed gravity [µGal]')
            pen = pg.mkPen(color='b')
            flags_keep_obs = obs_df['keep_obs'].values
            symbol_brushes = []
            for flag in flags_keep_obs:
                if flag:
                    symbol_brushes.append(self.BRUSH_ACTIVE_OBS)
                else:
                    symbol_brushes.append(self.BRUSH_INACTIVE_OBS)

            plot_offset_mgal = round(g_mugal.mean() / 1000)
            plot_offset_mugal = plot_offset_mgal * 1000

            # setup data: g
            if setup_df is not None and self.checkBox_obs_plot_setup_data.isChecked():
                self.plot_xy_data(self.plot_obs_g, setup_df['epoch_unix'].values,
                                  setup_df['g_mugal'].values - plot_offset_mugal,
                                  plot_name='setup', color='k', symbol='x', symbol_size=25)

            # Type of 'self.plot_obs_g_data_item': PlotDataItem
            self.plot_obs_g_data_item = self.plot_obs_g.plot(obs_epoch_timestamps, g_mugal - plot_offset_mugal,
                                                             name=f'Ref.: {ref_height_name}',
                                                             pen=pen, symbol='o', symbolSize=10,
                                                             symbolBrush=symbol_brushes)
            self.plot_obs_g_data_item.sigPointsClicked.connect(self.on_observation_plot_data_item_clicked)
            self.plot_obs_g.showGrid(x=True, y=True)
            self.plot_obs_g.setTitle(f'Observed gravity [µGal] + {plot_offset_mgal:.1f} mGal')
            self.plot_obs_g.setLabel(axis='left', text=f'g [µGal] + {plot_offset_mgal:.0f} mGal')
            self.plot_obs_g.autoRange()

            # Standard deviation of gravity g [µGal]
            self.plot_obs_sd_g.clear()
            # setup data: sd_g
            if setup_df is not None and self.checkBox_obs_plot_setup_data.isChecked():
                self.plot_xy_data(self.plot_obs_sd_g, setup_df['epoch_unix'].values, setup_df['sd_g_mugal'].values,
                                  plot_name='setup', color='k', symbol='x', symbol_size=25)
            self.plot_xy_data(self.plot_obs_sd_g, obs_epoch_timestamps, sd_g_mugal, plot_name='SD [µGal]',
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
            # - Tidal corrections:
            self.plot_xy_data(self.plot_obs_corrections, obs_epoch_timestamps, corr_tide,
                              plot_name=f'tides ({corr_tide_name})', color='b', symbol='o', symbol_size=10)
            # - Atm. pressure corrections (if available):
            if corr_atm_pres_name != '' and corr_atm_pres_name != 'no_atm_pres_corr':
                self.plot_xy_data(self.plot_obs_corrections, obs_epoch_timestamps, corr_atm_pres,
                                  plot_name=f'atm. pres ({corr_atm_pres_name})', color='r', symbol='o', symbol_size=10)

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
        try:
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
            index = self.tableView_observations.model().index(spot_item._index,
                                                              self.observation_model._data_column_names.index(
                                                                  'keep_obs'))
            self.observation_model.dataChanged.emit(index, index,
                                                    [9999])  # Call connected method "on_observation_model_data_changed"
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))

    def on_observation_model_data_changed(self, topLeft, bottomRight, role):
        """Invoked whenever data in the observation table view changed."""
        if IS_VERBOSE:
            # print('TableModel: dataChanged: ', topLeft, bottomRight, role)
            pass
        if topLeft == bottomRight and 9999 in role:  # Only one item selected and role=9999 => keep obs flag changed
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
            pass  # More than one item selected/changed.

    def on_tree_widget_item_changed(self, item, column):
        """Invoked whenever an item in th observation tree view is changed, e.g. if the check-state changes."""
        self.treeWidget_observations.blockSignals(True)  # To avoid recursive effects
        # if IS_VERBOSE:
        #     print('TreeView: itemChanged: ', item, column)
        #     print(f' - CheckState of item "{item.text(0)}": {item.checkState(0)}')
        flag_checked_state = checked_state_to_bool(item.checkState(0))
        # Is parent (survey) or child (setup):
        if item.parent() is None:  # Is survey item
            self.on_obs_tree_widget_item_selected()
        else:  # Is a setup item
            # Update table view and data in dataframe:
            survey_name = item.parent().text(0)
            setup_id = int(item.text(0))
            self.campaign.surveys[survey_name].activate_setup(setup_id, flag_checked_state)
            self.survey_model.emit_data_changed_survey(survey_name)
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
            self.plot_observations(survey_name)
        else:
            if IS_VERBOSE:
                print('No item or multiple items selected!')

    def update_obs_table_view(self, survey_name: str, setup_id: int):
        """Update the observation table view according to the selected survey and instrument setup."""
        # if IS_VERBOSE:
        #    print(f'survey name: {survey_name}; setup ID: {setup_id}')
        # Update the observation table view model according to the selected
        self.observation_model.update_view_model(survey_name,
                                                 setup_id,
                                                 gui_simple_mode=self.dlg_options.gui_simple_mode)  # Show added survey in table
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

    def survey_tree_widget_collapse_all(self):
        """Collapse all items in the survey tree widget"""
        self.treeWidget_observations.collapseAll()

    def survey_tree_widget_expand_all(self):
        """Expand all items in the survey tree widget"""
        self.treeWidget_observations.expandAll()

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
                        child.setCheckState(0, Qt.Unchecked)
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
        self.update_stations_map(auto_range=False)

    @pyqtSlot(int)  # Required, because 2 signals are emitted and one (int) has to be selected!
    def on_checkBox_filter_observed_stat_only_toggled(self, state, auto_range_stations_plot=False):
        """Event handler for the filter observed stations only checkbox."""
        self.proxy_station_model.setFilterKeyColumn(self.campaign.stations._STAT_DF_COLUMNS.index('is_observed'))
        if state == Qt.Checked:
            self.proxy_station_model.setFilterFixedString('True')
        else:
            self.proxy_station_model.setFilterFixedString('')
        self.update_stations_map(auto_range=auto_range_stations_plot)

    @pyqtSlot()
    def exit_application(self):
        """Exit the application."""
        sys.exit()

    def on_menu_file_options(self):
        """Launch dialog with general program options."""
        return_value = self.dlg_options.exec()
        if return_value == QDialog.Accepted:
            try:
                self.apply_options()
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
            else:
                pass
        else:
            pass  # Do nothing
            # self.statusBar().showMessage(f"Changes in options not applied.")

    def apply_options(self):
        """Apply options that are set in the options dialog."""
        if self.campaign is None:
            return
        # Simple or advanced GUI appearance:
        # - Update the results tab in order to change the appearance.
        self.update_results_tab()
        # - Stations Tab:
        self.set_up_station_view_model()
        self.enable_station_view_options_based_on_model()
        self.set_up_proxy_station_model()
        # - Observations tab:
        survey_name, setup_id = self.get_obs_tree_widget_selected_item()
        self.update_obs_table_view(survey_name, setup_id)
        self.update_setup_table_view(survey_name, setup_id)
        self.set_up_survey_view_model()

    def on_menu_help_about(self):
        """Launch the about dialog."""
        _ = self.dlg_about.exec()

    def on_menu_observations_corrections(self):
        """Launch diaglog to select and apply observation corrections."""
        return_value = self.dlg_corrections.exec()
        if return_value == QDialog.Accepted and self.campaign:
            try:
                self.apply_observation_corrections()
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f"Error: No observation corrections applied.")
            else:
                # Load survey from campaign data to observations vie model:
                self.observation_model.load_surveys(self.campaign.surveys)
                self.on_obs_tree_widget_item_selected()
                self.statusBar().showMessage(f"Observation corrections applied.")
        else:
            self.statusBar().showMessage(f"No observation corrections applied.")

    def apply_observation_corrections(self):
        """Apply observation corrections according to the selected settings."""
        flag_selection_ok = True
        error_msg = ''
        tide_corr_timeseries_interpol_method = ''
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
            QMessageBox.critical('Error!', error_msg)

        if self.dlg_corrections.radioButton_corr_tides_no_correction.isChecked():
            target_tide_corr = 'no_tide_corr'
        elif self.dlg_corrections.radioButton_corr_tides_cg5_model.isChecked():
            target_tide_corr = 'cg5_longman1959'
        elif self.dlg_corrections.radioButton_corr_tides_longman1959.isChecked():
            target_tide_corr = 'longman1959'
        elif self.dlg_corrections.radioButton_corr_tides_time_series.isChecked():
            target_tide_corr = 'from_time_series'
            tide_corr_timeseries_interpol_method = self.dlg_corrections.comboBox_tides_interpolation_method.currentText()
        else:
            flag_selection_ok = False
            error_msg = f'Invalid selection of tidal correction in GUI (observation corrections dialog).'
            QMessageBox.critical('Error!', error_msg)

        if self.dlg_corrections.checkBox_corrections_atm_pressure.isChecked():
            target_atm_pres_corr = 'iso_2533_1975'
        else:
            target_atm_pres_corr = 'no_atm_pres_corr'

        atm_pres_admittance = self.dlg_corrections.doubleSpinBox_atm_pres_admittance.value()

        if flag_selection_ok:
            try:
                self.campaign.reduce_observations_in_all_surveys(
                    target_ref_height=target_ref_height,
                    target_tide_corr=target_tide_corr,
                    target_atm_pres_corr=target_atm_pres_corr,
                    atm_pres_admittance=atm_pres_admittance,
                    tide_corr_timeseries_interpol_method=tide_corr_timeseries_interpol_method,
                    verbose=IS_VERBOSE)
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))


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
            self.menu_Add_Stations.setEnabled(True)
            self.action_Export_Results.setEnabled(True)
            self.action_Save_Campaign.setEnabled(True)
            self.action_Change_output_directory.setEnabled(True)
            self.action_Change_Campaign_name.setEnabled(True)

            # Set up GUI (models and widgets):
            # - Stations Tab:
            self.set_up_station_view_model()
            self.enable_station_view_options_based_on_model()
            self.set_up_proxy_station_model()
            self.update_comboBox_stations_selection_surrvey(survey_names=self.campaign.survey_names)
            self.update_stations_map(auto_range=True)
            # - Observations Tab:
            self.set_up_observation_view_model()
            self.enable_menu_observations_based_on_campaign_data()
            self.populate_survey_tree_widget()
            self.set_up_setup_view_model()
            self.plot_observations(survey_name=None)  # wipe observations plot
            self.set_up_survey_view_model()
            # - Results tab:
            self.set_up_results_stations_view_model()
            self.set_up_results_correlation_matrix_view_model()
            self.set_up_results_observations_view_model()
            self.set_up_results_drift_view_model()
            self.set_up_results_vg_view_model()
            self.update_results_tab(select_latest_item=True)
            # - Corrections time series dialog:
            self.dlg_correction_time_series.reset_update_gui()

            self.statusBar().showMessage(f"New Campaign created (name: {self.campaign.campaign_name}, "
                                         f"output directory: {self.campaign.output_directory})")
            self.setWindowTitle('GravTools - Campaign: ' + self.campaign.campaign_name)

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

    def on_menu_file_load_stations_from_oesgn_table(self):
        """Load stations from OESGN table file."""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        oesgn_filename, _ = QFileDialog.getOpenFileName(self,
                                                               'Select OESGN table file',
                                                               self.campaign.output_directory,
                                                               "OESGN table file (*.tab)",
                                                               options=options)
        if oesgn_filename:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            oesgn_filename = QDir.toNativeSeparators(oesgn_filename)
            self.load_stations(file_type='oesgn', filename=oesgn_filename)

    def on_menu_file_load_stations_from_csv_file(self):
        """Load stations from a CSV file."""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        csv_filename, _ = QFileDialog.getOpenFileName(self,
                                                        'Select station CSV file',
                                                        self.campaign.output_directory,
                                                        "Station CSV file (*.csv)",
                                                        options=options)
        if csv_filename:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            csv_filename = QDir.toNativeSeparators(csv_filename)
            self.load_stations(file_type='csv', filename=csv_filename)

    def load_stations(self, file_type, filename):
        """Load stations from stations file to the campaign and update the GUI."""

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
            if file_type == 'oesgn':
                self.campaign.add_stations_from_oesgn_table_file(filename,
                                                                 is_datum=settings.INIT_OESGN_STATION_AS_DATUM,
                                                                 verbose=IS_VERBOSE)
            elif file_type == 'csv':
                self.campaign.add_stations_from_csv_file(filename, verbose=IS_VERBOSE)
            else:
                raise AssertionError(f'Unknown station file type: {file_type}!')
            self.campaign.synchronize_stations_and_surveys(verbose=IS_VERBOSE)
            self.refresh_stations_table_model_and_view()
            self.set_up_proxy_station_model()
            self.on_checkBox_filter_observed_stat_only_toggled(
                state=self.checkBox_filter_observed_stat_only.checkState())

        except FileNotFoundError:
            QMessageBox.critical(self, 'File not found error', f'"{filename}" not found.')
            self.statusBar().showMessage(f"No stations added.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
            self.statusBar().showMessage(f"No stations added.")
        else:
            self.enable_station_view_options_based_on_model()
            # Show observed stations only based on Checkbox state:
            self.on_checkBox_filter_observed_stat_only_toggled(self.checkBox_filter_observed_stat_only.checkState(),
                                                               auto_range_stations_plot=True)
            # self.update_stations_map() => Called in self.on_checkBox_filter_observed_stat_only_toggled() above!
            number_of_stations_added = self.campaign.stations.get_number_of_stations - number_of_stations_old

            # Re-calculate observation corrections, based on the new station data (VG is relevant for height
            # reduction!):
            try:
                self.apply_observation_corrections()
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f"Error: No observation corrections applied.")
            else:
                self.statusBar().showMessage(f"Observation corrections applied.")
            finally:
                # Update the observations table and plot (in case the reduced obs. changed):
                survey_name, setup_id = self.get_obs_tree_widget_selected_item()
                self.update_obs_table_view(survey_name, setup_id)
                self.plot_observations(survey_name)

                self.statusBar().showMessage(f"{number_of_stations_added} stations added.")

    def enable_station_view_options_based_on_model(self):
        """Enable or disable the station view options based on the number of stations in the model."""
        if len(self.station_model._data) > 0:
            self.groupBox_filter_options.setEnabled(True)
            self.groupBox_stations_map_view_options.setEnabled(True)
        else:
            self.groupBox_filter_options.setEnabled(False)
            self.groupBox_stations_map_view_options.setEnabled(False)

    def refresh_stations_table_model_and_view(self):
        """Refresh the station table model and the respective table view."""
        self.station_model.load_stat_df(self.campaign.stations.stat_df)
        self.station_model.layoutChanged.emit()  # Refresh the table view, after changing the model's size.
        self.tableView_Stations.resizeColumnsToContents()

    @pyqtSlot()
    def set_up_station_view_model(self):
        """Set up station data view model and show station data table view."""
        try:
            self.station_model = StationTableModel(self.campaign.stations.stat_df,
                                                   gui_simple_mode=self.dlg_options.gui_simple_mode)
            self.station_model.dataChanged.connect(self.on_station_model_data_changed)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No stations available!')
            self.statusBar().showMessage(f"No stations available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.connect_station_model_to_table_view()  # Set view
            self.tableView_Stations.resizeColumnsToContents()  # TODO This line takes a lot of time!!!!!!!!!!!!
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

    def set_up_survey_view_model(self):
        """Set up the survey view model."""
        # Set model:
        try:
            self.survey_model = SurveyTableModel(self.campaign.surveys,
                                                   gui_simple_mode=self.dlg_options.gui_simple_mode)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No surveys available!')
            self.statusBar().showMessage(f"No surveys available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.tableView_surveys.setModel(self.survey_model)
            self.tableView_surveys.resizeColumnsToContents()

    @pyqtSlot()
    def on_menu_file_load_survey_from_cg5_observation_file(self):
        """Launch file selection dialog to select a CG5 observation file and load the data to the campaign."""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        cg5_obs_file_filenames, _ = QFileDialog.getOpenFileNames(self,
                                                               'Select CG5 observation file',
                                                               self.campaign.output_directory,
                                                               "CG5 observation file (*.TXT)",
                                                               options=options)
        if not cg5_obs_file_filenames:
            self.statusBar().showMessage(f"No survey data added.")
            return

        added_surveys_list = []
        # Add surveys to the campaign:
        for cg5_obs_file_filename in cg5_obs_file_filenames:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the
            # underlying operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            cg5_obs_file_filename = QDir.toNativeSeparators(cg5_obs_file_filename)
            # Add survey data to Campaign:
            try:
                new_cg5_survey = Survey.from_cg5_obs_file(cg5_obs_file_filename)
                if new_cg5_survey.obs_tide_correction_type == 'unknown':
                    raise RuntimeError('Type of tidal correction is unknown!')
                self.campaign.add_survey(survey_add=new_cg5_survey, verbose=IS_VERBOSE)
            except Exception as e:
                QMessageBox.critical(self, 'Error!', f'Error while loading {cg5_obs_file_filename}: ' + str(e))
                # self.statusBar().showMessage(f"No survey data added.")
                continue
            else:
                added_surveys_list.append(new_cg5_survey.name + f' ({new_cg5_survey.get_number_of_observations()} obs.)')

        if  not added_surveys_list:
            self.statusBar().showMessage(f"No survey data added.")
            return

        self.update_comboBox_stations_selection_surrvey(survey_names=self.campaign.survey_names)

        # Calculate reduced observation data:
        if settings.CALCULATE_REDUCED_OBS_WHEN_LOADING_DATA:
            try:
                self.apply_observation_corrections()
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f"Error: No observation corrections applied.")
            else:
                self.statusBar().showMessage(f"Observation corrections applied.")

        self.connect_station_model_to_table_view()
        self.campaign.synchronize_stations_and_surveys(verbose=IS_VERBOSE)
        self.refresh_stations_table_model_and_view()
        self.populate_survey_tree_widget()
        # Select the added survey in the tree view:
        for tree_item_idx in range(self.treeWidget_observations.topLevelItemCount()):
            if self.treeWidget_observations.topLevelItem(tree_item_idx).text(0) == new_cg5_survey.name:
                self.treeWidget_observations.topLevelItem(tree_item_idx).setSelected(True)
        self.enable_menu_observations_based_on_campaign_data()

        #
        self.set_up_proxy_station_model()  # Re-connect the sort & filter proxy model to the station view.
        self.set_up_survey_view_model()
        self.on_checkBox_filter_observed_stat_only_toggled(
            state=self.checkBox_filter_observed_stat_only.checkState())
        self.enable_station_view_options_based_on_model()
        # Show observed stations only based on Checkbox state:
        self.on_checkBox_filter_observed_stat_only_toggled(self.checkBox_filter_observed_stat_only.checkState(),
                                                           auto_range_stations_plot=True)
        self.statusBar().showMessage(
            f'{len(added_surveys_list)} survey added to campaign: ' + ','.join(added_surveys_list))


    @pyqtSlot()
    def on_manu_observations_flag_observations(self):
        """Load observation list file to flag observations in current campaign accordingly.

        The observation list file contains a list of gravity meter observations in CSV format with the following columns:
          - survey_name: Name of the survey
          - obs_epoch: Observation reference epoch (datetime string)
          - station_name: Name of the station
          - keep_obs: flag that indicates whether this particular observation is active or not.

        The observations in the file are matched with those in the campaign based on `survey_name`, `obs_epoch` and
        `station_name`. All observations with `keep_obs = False` in the loaded file are flagged on the campaign.
        """
        # Get the filename of the observation list csv file:
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        obs_list_filename, _ = QFileDialog.getOpenFileName(self,
                                                           'Load observation list file',
                                                           self.campaign.output_directory,
                                                           "Gravtools observation list file(*.csv)",
                                                           options=options)
        if obs_list_filename:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            obs_list_filename = QDir.toNativeSeparators(obs_list_filename)
            # Flag observations in current campaign:
            try:
                flag_log_str = self.campaign.flag_observations_based_on_obs_list_csv_file(
                    obs_list_filename=obs_list_filename,
                    update_type='all_obs',
                    verbose=IS_VERBOSE)
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
                self.statusBar().showMessage(f"Error while flagging observations.")
            else:
                # Update data visualization in GUI:
                for tree_idx in range(0, self.treeWidget_observations.topLevelItemCount()):
                    survey_name_tree = self.treeWidget_observations.topLevelItem(tree_idx).text(0)
                    # Update observation view model:
                    self.observation_model.update_view_model(survey_name_tree,
                                                             setup_id=None,
                                                             gui_simple_mode=self.dlg_options.gui_simple_mode)
                    self.update_obs_tree_widgget_from_observation_model()
                survey_name, setup_id = self.get_obs_tree_widget_selected_item()
                self.update_obs_table_view(survey_name, setup_id)
                self.plot_observations(survey_name)
                # print flag log sting to Message Box:

                # QMessageBox.information(self, 'Update Information!', flag_log_str)
                res = ScrollMessageBox(QMessageBox.Information, 'Update Information!', flag_log_str)

    def on_results_vg_plot_type_radiobuttons_changed(self, checked):
        """Invoked whenever a radiobutton for selecting the VG plot type in the results tab is changed."""
        if not checked:
            return
        if self.radioButton_results_vg_plot_details.isChecked():
            self.checkBox_results_vg_plot_show_residuals.setEnabled(True)
            self.update_vg_plot()
        if self.radioButton_results_vg_plot_full_polynomial.isChecked():
            self.checkBox_results_vg_plot_show_residuals.setEnabled(False)
            self.update_vg_plot()

    def checkBox_results_vg_plot_show_residuals_state_changed(self):
        """Invoked whenever the state of the checkbox changes."""
        self.update_vg_plot()

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
        output_dir_name = QFileDialog.getExistingDirectory(self, 'Select a directory', initial_folder_path,
                                                           QFileDialog.ShowDirsOnly)

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

    def __init__(self, campaign_output_dir, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)
        # connect signals and slots:
        self.pushButton_select_oesgn_table.clicked.connect(self.select_oesgn_table_file_dialog)
        self.campaign_output_dir = campaign_output_dir

    def select_oesgn_table_file_dialog(self):
        """Open file selection dialog."""
        # initial_oesgn_file_path = self.lineEdit_oesgn_table_file_path.text()
        if self.lineEdit_oesgn_table_file_path.text():
            initial_oesgn_file_path = os.path.dirname(os.path.abspath(self.lineEdit_oesgn_table_file_path.text()))
        else:
            initial_oesgn_file_path = self.campaign_output_dir
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


class DialogGisExportSettings(QDialog, Ui_Dialog_gis_settings):
    """Dialog to define the gis export settings."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)

        # Connect signals and slots:
        self.pushButton_select_gis_output_dir.clicked.connect(self.get_output_directory_dialog)
        # Set up GUI widgets:
        self.lineEdit_output_subdir.setText(settings.GIS_RESULTS_OUTPUT_SUBDIR)

    def get_output_directory_dialog(self):
        """Open dialog to get the output directory."""
        if self.lineEdit_gis_output_dir.text():
            initial_folder_path = os.path.dirname(os.path.abspath(self.lineEdit_gis_output_dir.text()))
        else:
            try:
                initial_folder_path = self.campaign_output_dir
            except:
                initial_folder_path = os.getcwd()
        output_dir_name = QFileDialog.getExistingDirectory(self, 'Select a directory', initial_folder_path,
                                                           QFileDialog.ShowDirsOnly)

        if output_dir_name:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            output_dir_name = QDir.toNativeSeparators(output_dir_name)
            if IS_VERBOSE:
                print(f'GIS output directory selectzed: {output_dir_name}')

        # Check, if path exists:
        if os.path.isdir(output_dir_name):
            self.lineEdit_gis_output_dir.setText(output_dir_name)
        else:
            QMessageBox.critical(self, "ERROR", f"The output directory does not exist: {output_dir_name}")
            self.lineEdit_gis_output_dir.setText('')


class DialogSetupData(QDialog, Ui_Dialog_setup_data):
    """Dialog to define setup data options."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)


class DialogAbout(QDialog, Ui_Dialog_about):
    """About dialog."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)


class DialogOptions(QDialog, Ui_Dialog_options):
    """Dialog for setting general program options."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)

    @property
    def gui_simple_mode(self) -> bool:
        """Returns flag whether to use the simple GUI mode."""
        if self.radioButton_gui_mode_simple.isChecked():
            return True
        elif self.radioButton_gui_mode_advanced.isChecked():
            return False
        else:
            raise AssertionError('Invalid GUI mode selected!')


class DialogExportResults(QDialog, Ui_Dialog_export_results):
    """Dialog to define the estimation settings."""

    def __init__(self, campaign, parent=None):
        super().__init__(parent)
        # Get data:
        self._lsm_runs = campaign.lsm_runs
        self.flag_observation_data_available = len(campaign.surveys) > 0
        self.flag_lsm_runs_available = len(campaign.lsm_run_times) > 0
        # Run the .setupUi() method to show the GUI
        self.setupUi(self)
        # Populate the combo box to select an LSM run (enable/disable groupBoxes accordingly):
        self.comboBox_select_lsm_run.clear()
        self.comboBox_select_lsm_run.addItems(['Current observation selection (no LSM run)'] + campaign.lsm_run_times)

        # Set the lineEdit with the export path:
        self.label_export_path_show.setText(campaign.output_directory)

        # Set initial item in comboBox (last entry/LSM run):
        idx = self.comboBox_select_lsm_run.count() - 1  # Index of last lsm run
        self.comboBox_select_lsm_run.setCurrentIndex(idx)

        # Select LSM run, etc.
        if self.flag_lsm_runs_available:
            lsm_run_idx = idx - 1
        else:
            lsm_run_idx = -1

        # Enable/disable GUI items and add lsm run comment:
        self.write_lsm_run_comment_to_gui(lsm_run_idx)
        self.enable_gui_widgets_based_on_lsm_run_selection(lsm_run_idx)

        # connect signals and slots:
        self.comboBox_select_lsm_run.currentIndexChanged.connect(self.on_comboBox_select_lsm_run_current_index_changed)

        # Optional dependency for GIS data export:
        if _has_geopandas:
            self.groupBox_gis.setEnabled(True)
        else:
            self.groupBox_gis.setEnabled(False)
            self.checkBox_gis_write_shapefile.setChecked(Qt.Unchecked)

    @pyqtSlot(int)
    def on_comboBox_select_lsm_run_current_index_changed(self, index: int):
        """Invoked whenever the index of the selected item in the combobox changed."""
        lsm_run_idx = index - 1
        self.write_lsm_run_comment_to_gui(lsm_run_idx)
        self.enable_gui_widgets_based_on_lsm_run_selection(lsm_run_idx)

    def write_lsm_run_comment_to_gui(self, lsm_run_idx):
        """Writes the lsm run comment of the selected lsm run to the GUI line edit."""
        if lsm_run_idx > 1:
            self.label_export_comment_show.setText(self.get_lsm_run_comment(lsm_run_idx))
        else:
            self.label_export_comment_show.setText('')

    def get_lsm_run_comment(self, lsm_run_idx: int):
        """Returns the lsm run comment of the run with the specified index."""
        try:
            return self._lsm_runs[lsm_run_idx].comment
        except:
            return ''

    def enable_gui_widgets_based_on_lsm_run_selection(self, lsm_run_idx: int):
        """Enable or disable GUI elements based on the lsm run selection."""

        flag_lsm_run_selected = self.flag_lsm_runs_available and (lsm_run_idx > -1)

        if flag_lsm_run_selected:
            self.groupBox_other_files.setEnabled(True)
            self.groupBox_nsb_file.setEnabled(True)
            self.groupBox_observation_list.setEnabled(True)
            if _has_geopandas:
                self.groupBox_gis.setEnabled(True)
        else:
            self.groupBox_gis.setEnabled(False)
            self.groupBox_other_files.setEnabled(False)
            self.groupBox_nsb_file.setEnabled(False)
            self.groupBox_observation_list.setEnabled(False)
            if self.flag_observation_data_available:
                self.groupBox_observation_list.setEnabled(True)
        if flag_lsm_run_selected or self.flag_observation_data_available:
            self.buttonBox.buttons()[0].setEnabled(True)  # OK button in buttonBox
        else:
            self.buttonBox.buttons()[0].setEnabled(False)  # OK button in buttonBox

        if flag_lsm_run_selected:  # => lsm_run_idx >= 0
            lsm_method = self._lsm_runs[lsm_run_idx].lsm_method
            if lsm_method == 'LSM_diff' or lsm_method == 'LSM_non_diff' or lsm_method == 'MLR_BEV':  # network adjust.
                self.groupBox_nsb_file.setEnabled(True)
                self.checkBox_save_vg_plot_png.setEnabled(False)
            elif lsm_method == 'VG_LSM_nondiff':  # VG estimation
                self.groupBox_nsb_file.setEnabled(False)
                self.checkBox_save_vg_plot_png.setEnabled(True)
            else:
                raise AssertionError(f'Invalid LSM method: {lsm_method}!')

def main():
    """Main program to start the GUI."""
    # Create the application
    app = QApplication(sys.argv)

    # Create and show the application's main window
    main_window = MainWindow()
    main_window.show()
    # Run the application's main loop:
    sys.exit(
        app.exec())  # exit or error code of Qt (app.exec_) is passed to sys.exit. Terminates pgm with standard python method


if __name__ == "__main__":
    """Main Program."""
    main()
