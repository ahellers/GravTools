"""Dialog for managing time series data for observation corrections.

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

import os
import pytz
import datetime as dt

from PyQt5.QtWidgets import QDialog, QMessageBox, QTreeWidgetItem
import pyqtgraph as pg

from gravtools.gui.dialog_correction_time_series import Ui_DialogCorrectionTimeSeries
from gravtools.gui.dlg_load_tsf_file import DialogLoadTsfFile
from gravtools import __version__
from gravtools import settings

# Set splitter stretch:
_SPLITTER_STRETCH_LEFT = 1
_SPLITTER_STRETCH_RIGHT = 10


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


class DialogCorrectionTimeSeries(QDialog, Ui_DialogCorrectionTimeSeries):
    """Dialog for managing time series data for observation corrections."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # Run the .setupUi() method to show the GUI
        self.setupUi(self)

        # Set the relative size of the left and right layout of the splitter:
        self.splitter.setStretchFactor(_SPLITTER_STRETCH_LEFT, _SPLITTER_STRETCH_RIGHT)

        # Connect signals and slots:
        self.pushButton_load_tsf_file.pressed.connect(self.launch_dlg_load_tsf_file)
        self.pushButton_collapse.pressed.connect(self.collapse_all_tree_widget)
        self.pushButton_expand.pressed.connect(self.expand_all_tree_widget)
        self.treeWidget_selection.currentItemChanged.connect(self.tree_widget_current_item_changed)
        self.pushButton_delete_selection.pressed.connect(self.delete_selected_item)

        # Attributes
        self.dlg_load_tsf_file = DialogLoadTsfFile(self)

        # Other initial settings:
        self.reset_stacked_widget()
        self.set_up_time_series_plot_widget()

    def launch_dlg_load_tsf_file(self):
        """Launch dialog to load TSF file."""
        if self.parent().campaign is None:
            QMessageBox.critical(self, 'Error!', 'No Campaign available! Please create a campaign first.')
            return

        self.check_correction_time_series_object()

        # Launch file selection dialog:
        return_value = self.dlg_load_tsf_file.exec()
        if return_value != QDialog.Accepted:
            return  # Do not load any data

        # Get filename and other parameters from the file selection dialog:
        filename_tsf = self.dlg_load_tsf_file.lineEdit_filename.text()
        survey_name = self.dlg_load_tsf_file.lineEdit_survey_name.text()
        filter_location = self.dlg_load_tsf_file.lineEdit_filter_location.text()
        filter_instrument = self.dlg_load_tsf_file.lineEdit_filter_instrument.text()
        filter_data_type = self.dlg_load_tsf_file.lineEdit_filter_data_type.text()
        overwrite_channel = self.dlg_load_tsf_file.checkBox_overwrite_channel.isChecked()

        is_effect = self.dlg_load_tsf_file.radioButton_effect.isChecked()
        is_correction = self.dlg_load_tsf_file.radioButton_correction.isChecked()

        if not os.path.exists(filename_tsf):
            QMessageBox.critical(self, 'Error!',
                                 f"The selected TSF file ({filename_tsf}) does not exist!")
            return

        # Load the TSF file:
        try:
            self.parent().campaign.correction_time_series.load_tfs_file(survey_name=survey_name,
                                                                        filename_tsf=filename_tsf,
                                                                        location=filter_location,
                                                                        instrument=filter_instrument,
                                                                        data_type=filter_data_type,
                                                                        overwrite_channel=overwrite_channel,
                                                                        is_correction=is_correction)
        except Exception as e:
            QMessageBox.critical(self, 'Error!', str(e))

        self.reset_update_gui()

    def check_correction_time_series_object(self, parent=None):
        """Adds a correction time series object to the campaign, if missing."""
        if not hasattr(self.parent().campaign, 'correction_time_series'):
            self.parent().campaign.add_empty_correction_time_series()
            if parent is None:
                parent = self
            QMessageBox.warning(parent, 'Warning!', f'The campaign "{self.parent().campaign.campaign_name}" had no '
                                                  f'"correction_time_series" attribute. This '
                                                  f'may happens when loading campaign data from previous GravTools '
                                                  f'versions. Current version: {__version__}. '
                                                  f'Campaign data version: {self.parent().campaign.gravtools_version}. '
                                                  f'However, "correction_time_series" was added to the campaign. Please'
                                                  f'save the campaign to keep the changes!')

    def reset_update_gui(self):
        """Update the GUI after changing data, e.g. loading from TSF files or deleting data."""
        self.reset_stacked_widget()
        self.populate_survey_tree_widget()
        self.update_timer_series_plot()

    def populate_survey_tree_widget(self):
        """Populate the survey tree widget."""

        correction_time_series = self.parent().campaign.correction_time_series

        # Delete existing items:
        self.treeWidget_selection.blockSignals(True)
        self.treeWidget_selection.clear()

        # Add new items:
        # Parent items (surveys):
        for survey_name, survey_corrections in correction_time_series.surveys.items():
            parent = QTreeWidgetItem(self.treeWidget_selection)
            parent.setText(0, survey_name)
            for station_name, station_correction in survey_corrections.stations.items():
                child = QTreeWidgetItem(parent)
                child.setText(0, station_name)

        self.treeWidget_selection.show()
        self.treeWidget_selection.blockSignals(False)
        self.treeWidget_selection.expandToDepth(0)  # Expand items to a certain depth

    def tree_widget_current_item_changed(self, tree_item, column):
        """Whenever the current item changes in the selection tree widget."""
        # Show item details in stacked widget according to selection:
        depth = get_depth(tree_item)
        if depth == 1:
            survey_name = tree_item.text(0)
            self.stacked_widget_show_survey_details(survey_name=survey_name)
        elif depth == 2:
            station_name = tree_item.text(0)
            survey_name = tree_item.parent().text(0)
            self.stacked_widget_show_station_details(survey_name=survey_name, station_name=station_name)
            self.update_timer_series_plot()
        else:
            pass

    def stacked_widget_show_station_details(self, survey_name, station_name):
        """Show station details in the stacked widget."""
        self.stackedWidget_item_details.setCurrentIndex(2)
        correction_time_series = self.parent().campaign.correction_time_series
        tidal_correction = correction_time_series.surveys[survey_name].stations[station_name].tidal_correction
        self.label_station_station_name.setText(station_name)
        self.label_station_data_source.setText(tidal_correction.data_source)
        self.label_station_description.setText(tidal_correction.description)
        self.label_station_unit.setText(tidal_correction.unit)
        self.label_station_created.setText(tidal_correction.created_datetime_utc_str)
        self.label_station_model_type.setText(tidal_correction.model_type)
        self.label_station_start.setText(tidal_correction.start_datetime_str)
        self.label_station_end.setText(tidal_correction.end_datetime_str)
        self.label_station_duration.setText(tidal_correction.duration_dhms_str)
        self.label_station_nuber_of_datapoints.setText(str(tidal_correction.number_of_datapoints))

    def stacked_widget_show_survey_details(self, survey_name):
        """Show survey details in the stacked widget."""

        correction_time_series = self.parent().campaign.correction_time_series
        self.stackedWidget_item_details.setCurrentIndex(1)
        self.label_survey_name.setText(survey_name)
        self.label_num_of_stations.setText(str(correction_time_series.surveys[survey_name].number_of_stations))

    def collapse_all_tree_widget(self):
        """Collapse all items in selection tree widget."""
        self.treeWidget_selection.collapseAll()

    def expand_all_tree_widget(self):
        """Expand all items in tree widget."""
        self.treeWidget_selection.expandAll()

    def reset_stacked_widget(self):
        """Reset the stacked widget by displaying page with index 0"""
        self.stackedWidget_item_details.setCurrentIndex(0)

    def delete_selected_item(self):
        """Delete the selected item in the tree widget"""
        correction_time_series = self.parent().campaign.correction_time_series
        tree_item = self.treeWidget_selection.currentItem()
        depth = get_depth(tree_item)
        if depth == 1:  # Survey
            survey_name = tree_item.text(0)
            msg_text = f'Survey {survey_name}'
            reply = QMessageBox.question(self,
                                         'Delete survey corrections?',
                                         msg_text,
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                try:
                    correction_time_series.delete_survey_correction(survey_name)
                except KeyError:
                    QMessageBox.critical(self, 'Error!', f'No corrections for survey "{survey_name}" available!')
                except Exception as e:
                    QMessageBox.critical(self, 'Error!', str(e))
                else:
                    self.reset_update_gui()
            else:
                return
        elif depth == 2:  # Station
            station_name = tree_item.text(0)
            survey_name = tree_item.parent().text(0)
            msg_text = f'Station {survey_name}'
            reply = QMessageBox.question(self,
                                         'Delete station corrections?',
                                         msg_text,
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                try:
                    correction_time_series.surveys[survey_name].delete_station_correction(station_name)
                except KeyError:
                    QMessageBox.critical(self, 'Error!', f'No corrections for station "{station_name}" in survey "{survey_name}" available!')
                except Exception as e:
                    QMessageBox.critical(self, 'Error!', str(e))
                else:
                    self.reset_update_gui()
            else:
                return
        else:
            return

    def set_up_time_series_plot_widget(self):
        """Set up `self.graphicsLayoutWidget_results_vg_plot` widget."""
        self.glw_ts_plot = self.graphicsLayoutWidget_time_series
        self.glw_ts_plot.setBackground('w')  # white background color
        self.ts_plot = self.glw_ts_plot.addPlot(0, 0, name='ts_plot',
                                                      axisItems={'bottom': TimeAxisItem(orientation='bottom')})
        # self.ts_plot.setLabel(axis='bottom', text='Time')  # Takes too much space
        self.ts_plot.addLegend()
        self.ts_plot.showGrid(x=True, y=True)
        # self.ts_plot.setTitle('')  # Takes too much space

    def plot_time_series(self, data, timestamps, plot_name: str ='', y_label: str = ''):
        # Clear plot in any case:
        self.ts_plot.clear()
        self.plot_xy_data(self.ts_plot, timestamps, data, plot_name=plot_name, color='b', symbol='o', symbol_size=8)
        self.ts_plot.setLabel(axis='left', text=y_label)
        self.ts_plot.autoRange()

    @staticmethod
    def plot_xy_data(plot_item, x, y, plot_name, color='k', symbol='o', symbol_size=10):
        """Plot XY-data."""
        pen = pg.mkPen(color=color)
        plot_item.plot(x, y, name=plot_name, pen=pen, symbol=symbol, symbolSize=symbol_size, symbolBrush=(color))

    def update_timer_series_plot(self):
        """Update the time series plots according to the selection in the tree widget."""
        tree_item = self.treeWidget_selection.currentItem()
        if tree_item is None:
            return
        depth = get_depth(tree_item)
        if depth == 2:  # Time series data
            # Get the data:
            station_name = tree_item.text(0)
            survey_name = tree_item.parent().text(0)
            try:
                time_series = self.parent().campaign.correction_time_series.surveys[survey_name].stations[
                    station_name].tidal_correction
            except KeyError:
                return
            if time_series is None:
                return
            data = time_series.data
            timestamps = time_series.ref_time_unix
            # Plot the data:
            try:
                self.plot_time_series(data, timestamps, plot_name=time_series.description, y_label=time_series.unit)
            except Exception as e:
                QMessageBox.critical(self, 'Error!', str(e))
        else:
            self.ts_plot.clear()
            return


def get_subtree_nodes(tree_widget_item) -> list:
    """Returns al QTreeWidgetItems in a subtree rooted at a given node (recursive)."""
    nodes = []
    for i in range(tree_widget_item.childCount()):
        nodes.extend(get_subtree_nodes(tree_widget_item.child(i)))
    return nodes


def get_all_items(tree_widget) -> list:
    """Returns all QTreeWidgetItems in a given QTreeWidget."""
    all_items = []
    for i in range(tree_widget.topLevelItemCount()):
        top_item = tree_widget.topLevelItem(i)
        all_items.extend(get_subtree_nodes(top_item))
    return all_items


def get_depth(tree_widget_item) -> int:
    """Returns the depth of a QTreeWidgetItem within the QTreeWidget."""
    depth = 0
    while tree_widget_item is not None:
        depth += 1
        tree_widget_item = tree_widget_item.parent()
    return depth





