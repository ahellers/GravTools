import sys
import os
from PyQt5.QtWidgets import QApplication, QDialog, QMainWindow, QFileDialog, QMessageBox
from PyQt5.QtCore import QDir

from MainWindow import Ui_MainWindow
from dialog_new_campaign import Ui_Dialog_new_Campaign

from gravtools.models.survey import Campaign, Survey, Station


DEFAULT_OUTPUT_DIR = os.path.abspath(os.getcwd())  # Current working directory


class MainWindow(QMainWindow, Ui_MainWindow):
    """Main Window of the application."""

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

    def exit_application(self):
        """Exit the application."""
        sys.exit()

    def create_new_campaign(self, campaign_name=''):
        """Create a new and empty campaign."""
        print(f'New campaing: {campaign_name}')

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
        elif return_value == QDialog.Rejected:
            self.statusBar().showMessage(f"Canceled creating new campaign.")

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
        self.pushButton_change_output_directory.clicked.connect(self.get_output_directory)
        # self.accepted.connect(self.test_fun)

        # Initialize the name and the output directory from an existing campaign (is it exists):
        # self.old_campaign = old_campaign
        if isinstance(old_campaign, Campaign):  # not None!
            self.lineEdit_output_directory.setText(old_campaign.output_directory)
            self.lineEdit_campaign_name.setText(old_campaign.campaign_name)
        else:
            self.lineEdit_output_directory.setText(DEFAULT_OUTPUT_DIR)
            self.lineEdit_campaign_name.setText('')

    # def test_fun(self):
    #     print('test_fun')

    def get_output_directory(self):
        """Open dialog to get the output directory."""

        initial_folder_path = self.lineEdit_output_directory.text()
        print(initial_folder_path)
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


if __name__ == "__main__":
    """Main Program."""
    # Create the application
    app = QApplication(sys.argv)

    # Create and show the application's main window
    main_window = MainWindow()
    main_window.show()
    # Run the application's main loop
    sys.exit(app.exec())
