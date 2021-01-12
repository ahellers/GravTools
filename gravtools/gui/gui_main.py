import sys
import os
from PyQt5.QtWidgets import QApplication, QDialog, QMainWindow, QFileDialog, QMessageBox, QTreeWidgetItem, QHeaderView
from PyQt5.QtCore import QDir, QAbstractTableModel, Qt, QSortFilterProxyModel, pyqtSlot, QRegExp
from PyQt5 import QtGui

from MainWindow import Ui_MainWindow
from dialog_new_campaign import Ui_Dialog_new_Campaign
from dialog_load_stations import Ui_Dialog_load_stations

from gravtools.models.survey import Campaign, Survey, Station
from gui_models import StationModel

DEFAULT_OUTPUT_DIR = os.path.abspath(os.getcwd())  # Current working directory
DEFAULT_CG5_OBS_FILE_PATH = os.path.abspath(os.getcwd())  # Current working directory
IS_VERBOSE = True  # Define, whether screen output is enabled.


class MainWindow(QMainWindow, Ui_MainWindow):
    """Main Window of the application."""

    def __init__(self):
        """Initializer."""

        # Instance Attributes:
        self.campaign = None

        # GUI:
        super().__init__()
        self.setupUi(self)

        # Overwrite setting from ui file, if necessary:
        # ...

        # Connect signals and slots
        self.action_Exit.triggered.connect(self.exit_application)
        self.action_New_Campaign.triggered.connect(self.on_menu_file_new_campaign)
        self.action_Add_Stations.triggered.connect(self.on_menu_file_load_stations)
        # self.actionShow_Stations.triggered.connect(self.show_station_data)
        self.action_from_CG5_observation_file.triggered.connect(self.on_menu_file_load_survey_from_cg5_observation_file)
        self.lineEdit_filter_stat_name.textChanged.connect(self.on_lineEdit_filter_stat_name_textChanged)
        self.checkBox_filter_observed_stat_only.stateChanged.connect(self.on_checkBox_filter_observed_stat_only_toggled)

        # Set up GUI items and widgets:
        self.set_up_survey_tree_widget()
        # self.observations_splitter.setSizes([1000, 10])

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

    def populate_survey_tree_widget(self):
        """Populate the survey tree widget."""
        # Delete existing items:
        self.delete_all_items_from_survey_tree_widget()

        # Add new items:
        # Parent items (surveys):
        for survey_name, survey in self.campaign.surveys.items():
            # print (survey_name, survey)
            # survey.keep_survey
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
                # print(setup_id)
                # get all observations
                setup_keep_obs_flags = obs_df.loc[obs_df['setup_id'] == setup_id, 'keep_obs']
                setup_station_names = obs_df.loc[obs_df['setup_id'] == setup_id, 'station_name']
                num_of_obs_in_setup = len(setup_keep_obs_flags)

                # get station name (has to be unique within one setup!):
                if len(setup_station_names.unique()) > 1:
                    QMessageBox.critical(self, 'Error!', 'Setup {} includes observations to more than one '
                                                         'station ({})!'.format(setup_id,
                                                         ', '.join(setup_station_names.unique().tolist())))
                    self.delete_all_items_from_survey_tree_widget()
                    return
                else:
                    setup_station_name = setup_station_names.unique()[0]
                    child = QTreeWidgetItem(parent)
                    child.setFlags(child.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsTristate )
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

    def create_new_campaign(self, campaign_name=''):
        """Create a new and empty campaign."""
        print(f'New campaing: {campaign_name}')

    @pyqtSlot()
    def on_menu_file_new_campaign(self):
        """Launching dialog to create a new campaign."""
        dlg = DialogNewCampaign(old_campaign=self.campaign)
        return_value = dlg.exec()
        if return_value == QDialog.Accepted:
            # print('Accepted')

            # Create Campaign object:
            self.campaign = Campaign(
                campaign_name=dlg.lineEdit_campaign_name.text(),
                output_directory=dlg.lineEdit_output_directory.text(),
                surveys=None,  # Always use non-mutable default arguments!
                stations=None,  # Always use non-mutable default arguments!
            )
            # Activate main menu items:
            self.menuAdd_Survey.setEnabled(True)
            self.action_Add_Stations.setEnabled(True)

            self.statusBar().showMessage(f"New Campaign created (name: {self.campaign.campaign_name}, "
                                         f"output directory: {self.campaign.output_directory})")
            self.setWindowTitle('GravTools - Campaign: ' + self.campaign.campaign_name)

            self.set_up_station_view_model()
            # self.set_up_proxy_station_model()

        elif return_value == QDialog.Rejected:
            self.statusBar().showMessage(f"Canceled creating new campaign.")

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
            self.station_model = StationModel(self.campaign.stations.stat_df)
        except AttributeError:
            QMessageBox.warning(self, 'Warning!', 'No stations available!')
            self.statusBar().showMessage(f"No stations available.")
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.connect_station_model_to_table_view()  # Set view
            self.tableView_Stations.resizeColumnsToContents()
            self.statusBar().showMessage(f"{self.campaign.number_of_stations} stations in current project.")

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
                                                               'Select OESGN table file',
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
                    self.connect_station_model_to_table_view()  # Disconnect sort & filter proxy model from station view.
                    self.campaign.synchronize_stations_and_surveys(verbose=IS_VERBOSE)
                    self.refresh_stations_table_model_and_view()
                    self.populate_survey_tree_widget()
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
        ## print(initial_folder_path)
        output_dir_name = QFileDialog.getExistingDirectory(self, 'Select a directory', initial_folder_path)

        if output_dir_name:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            output_dir_name = QDir.toNativeSeparators(output_dir_name)
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
        # print(initial_oesgn_file_path)
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        oesgn_filename, _ = QFileDialog.getOpenFileName(self, 'Select OESGN table file', initial_oesgn_file_path,
                                                        "OESGN table (*.TAB)", options=options)
        if oesgn_filename:
            # Returns pathName with the '/' separators converted to separators that are appropriate for the underlying
            # operating system.
            # On Windows, toNativeSeparators("c:/winnt/system32") returns "c:\winnt\system32".
            oesgn_filename = QDir.toNativeSeparators(oesgn_filename)
            # print(oesgn_filename)
            self.lineEdit_oesgn_table_file_path.setText(oesgn_filename)


if __name__ == "__main__":
    """Main Program."""
    # Create the application
    app = QApplication(sys.argv)

    # Create and show the application's main window
    main_window = MainWindow()
    main_window.show()
    # Run the application's main loop
    sys.exit(
        app.exec())  # exit or error code of Qt (app.exec_) is passed to sys.exit. Terminates pgm with standard python method
