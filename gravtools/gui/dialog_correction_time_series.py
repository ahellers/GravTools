# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gravtools/gui/dialog_correction_time_series.ui'
#
# Created by: PyQt5 UI code generator 5.15.6
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DialogCorrectionTimeSeries(object):
    def setupUi(self, DialogCorrectionTimeSeries):
        DialogCorrectionTimeSeries.setObjectName("DialogCorrectionTimeSeries")
        DialogCorrectionTimeSeries.resize(1214, 600)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(DialogCorrectionTimeSeries)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.splitter = QtWidgets.QSplitter(DialogCorrectionTimeSeries)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setHandleWidth(4)
        self.splitter.setChildrenCollapsible(False)
        self.splitter.setObjectName("splitter")
        self.layoutWidget = QtWidgets.QWidget(self.splitter)
        self.layoutWidget.setObjectName("layoutWidget")
        self.verticalLayout_left = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout_left.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout_left.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_left.setObjectName("verticalLayout_left")
        self.treeWidget_selection = QtWidgets.QTreeWidget(self.layoutWidget)
        self.treeWidget_selection.setHeaderHidden(True)
        self.treeWidget_selection.setColumnCount(1)
        self.treeWidget_selection.setObjectName("treeWidget_selection")
        self.treeWidget_selection.headerItem().setText(0, "1")
        self.treeWidget_selection.header().setVisible(False)
        self.verticalLayout_left.addWidget(self.treeWidget_selection)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.pushButton_delete_selection = QtWidgets.QPushButton(self.layoutWidget)
        self.pushButton_delete_selection.setObjectName("pushButton_delete_selection")
        self.horizontalLayout_2.addWidget(self.pushButton_delete_selection)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.verticalLayout_left.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.pushButton_collapse = QtWidgets.QPushButton(self.layoutWidget)
        self.pushButton_collapse.setObjectName("pushButton_collapse")
        self.horizontalLayout_3.addWidget(self.pushButton_collapse)
        self.pushButton_expand = QtWidgets.QPushButton(self.layoutWidget)
        self.pushButton_expand.setObjectName("pushButton_expand")
        self.horizontalLayout_3.addWidget(self.pushButton_expand)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.verticalLayout_left.addLayout(self.horizontalLayout_3)
        self.layoutWidget1 = QtWidgets.QWidget(self.splitter)
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.verticalLayout_right = QtWidgets.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout_right.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout_right.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_right.setObjectName("verticalLayout_right")
        self.stackedWidget_item_details = QtWidgets.QStackedWidget(self.layoutWidget1)
        self.stackedWidget_item_details.setObjectName("stackedWidget_item_details")
        self.page_empty = QtWidgets.QWidget()
        self.page_empty.setObjectName("page_empty")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.page_empty)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_6 = QtWidgets.QLabel(self.page_empty)
        self.label_6.setObjectName("label_6")
        self.verticalLayout.addWidget(self.label_6)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.stackedWidget_item_details.addWidget(self.page_empty)
        self.page_survey = QtWidgets.QWidget()
        self.page_survey.setObjectName("page_survey")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.page_survey)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label = QtWidgets.QLabel(self.page_survey)
        self.label.setObjectName("label")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.label_survey_name = QtWidgets.QLabel(self.page_survey)
        self.label_survey_name.setText("")
        self.label_survey_name.setObjectName("label_survey_name")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.label_survey_name)
        self.label_2 = QtWidgets.QLabel(self.page_survey)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.label_num_of_stations = QtWidgets.QLabel(self.page_survey)
        self.label_num_of_stations.setText("")
        self.label_num_of_stations.setObjectName("label_num_of_stations")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.label_num_of_stations)
        self.verticalLayout_5.addLayout(self.formLayout)
        self.stackedWidget_item_details.addWidget(self.page_survey)
        self.page_station = QtWidgets.QWidget()
        self.page_station.setObjectName("page_station")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.page_station)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.label_3 = QtWidgets.QLabel(self.page_station)
        self.label_3.setObjectName("label_3")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.label_station_data_source = QtWidgets.QLabel(self.page_station)
        self.label_station_data_source.setText("")
        self.label_station_data_source.setObjectName("label_station_data_source")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.label_station_data_source)
        self.label_4 = QtWidgets.QLabel(self.page_station)
        self.label_4.setObjectName("label_4")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_4)
        self.label_station_description = QtWidgets.QLabel(self.page_station)
        self.label_station_description.setText("")
        self.label_station_description.setObjectName("label_station_description")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.label_station_description)
        self.label_5 = QtWidgets.QLabel(self.page_station)
        self.label_5.setObjectName("label_5")
        self.formLayout_2.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.label_5)
        self.label_station_unit = QtWidgets.QLabel(self.page_station)
        self.label_station_unit.setText("")
        self.label_station_unit.setObjectName("label_station_unit")
        self.formLayout_2.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.label_station_unit)
        self.label_10 = QtWidgets.QLabel(self.page_station)
        self.label_10.setObjectName("label_10")
        self.formLayout_2.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.label_10)
        self.label_station_model_type = QtWidgets.QLabel(self.page_station)
        self.label_station_model_type.setText("")
        self.label_station_model_type.setObjectName("label_station_model_type")
        self.formLayout_2.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.label_station_model_type)
        self.label_7 = QtWidgets.QLabel(self.page_station)
        self.label_7.setObjectName("label_7")
        self.formLayout_2.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label_7)
        self.label_station_start = QtWidgets.QLabel(self.page_station)
        self.label_station_start.setText("")
        self.label_station_start.setObjectName("label_station_start")
        self.formLayout_2.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.label_station_start)
        self.label_8 = QtWidgets.QLabel(self.page_station)
        self.label_8.setObjectName("label_8")
        self.formLayout_2.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.label_8)
        self.label_station_end = QtWidgets.QLabel(self.page_station)
        self.label_station_end.setText("")
        self.label_station_end.setObjectName("label_station_end")
        self.formLayout_2.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.label_station_end)
        self.label_9 = QtWidgets.QLabel(self.page_station)
        self.label_9.setObjectName("label_9")
        self.formLayout_2.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.label_9)
        self.label_station_duration = QtWidgets.QLabel(self.page_station)
        self.label_station_duration.setText("")
        self.label_station_duration.setObjectName("label_station_duration")
        self.formLayout_2.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.label_station_duration)
        self.label_12 = QtWidgets.QLabel(self.page_station)
        self.label_12.setObjectName("label_12")
        self.formLayout_2.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.label_12)
        self.label_station_nuber_of_datapoints = QtWidgets.QLabel(self.page_station)
        self.label_station_nuber_of_datapoints.setText("")
        self.label_station_nuber_of_datapoints.setObjectName("label_station_nuber_of_datapoints")
        self.formLayout_2.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.label_station_nuber_of_datapoints)
        self.label_11 = QtWidgets.QLabel(self.page_station)
        self.label_11.setObjectName("label_11")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_11)
        self.label_station_station_name = QtWidgets.QLabel(self.page_station)
        self.label_station_station_name.setText("")
        self.label_station_station_name.setObjectName("label_station_station_name")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.label_station_station_name)
        self.label_13 = QtWidgets.QLabel(self.page_station)
        self.label_13.setObjectName("label_13")
        self.formLayout_2.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_13)
        self.label_station_created = QtWidgets.QLabel(self.page_station)
        self.label_station_created.setText("")
        self.label_station_created.setObjectName("label_station_created")
        self.formLayout_2.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.label_station_created)
        self.verticalLayout_2.addLayout(self.formLayout_2)
        self.graphicsLayoutWidget_time_series = GraphicsLayoutWidget(self.page_station)
        self.graphicsLayoutWidget_time_series.setObjectName("graphicsLayoutWidget_time_series")
        self.verticalLayout_2.addWidget(self.graphicsLayoutWidget_time_series)
        self.stackedWidget_item_details.addWidget(self.page_station)
        self.verticalLayout_right.addWidget(self.stackedWidget_item_details)
        self.frame_actions = QtWidgets.QFrame(self.layoutWidget1)
        self.frame_actions.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_actions.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_actions.setObjectName("frame_actions")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.frame_actions)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.pushButton_load_tsf_file = QtWidgets.QPushButton(self.frame_actions)
        self.pushButton_load_tsf_file.setObjectName("pushButton_load_tsf_file")
        self.horizontalLayout.addWidget(self.pushButton_load_tsf_file)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem3)
        self.verticalLayout_4.addLayout(self.horizontalLayout)
        self.verticalLayout_right.addWidget(self.frame_actions)
        self.verticalLayout_3.addWidget(self.splitter)

        self.retranslateUi(DialogCorrectionTimeSeries)
        self.stackedWidget_item_details.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(DialogCorrectionTimeSeries)

    def retranslateUi(self, DialogCorrectionTimeSeries):
        _translate = QtCore.QCoreApplication.translate
        DialogCorrectionTimeSeries.setWindowTitle(_translate("DialogCorrectionTimeSeries", "Time series data"))
        self.pushButton_delete_selection.setToolTip(_translate("DialogCorrectionTimeSeries", "Delete selected item in the tree view."))
        self.pushButton_delete_selection.setText(_translate("DialogCorrectionTimeSeries", "Delete selected item"))
        self.pushButton_collapse.setToolTip(_translate("DialogCorrectionTimeSeries", "Toggle collapse tree view."))
        self.pushButton_collapse.setText(_translate("DialogCorrectionTimeSeries", "Collapse all"))
        self.pushButton_expand.setText(_translate("DialogCorrectionTimeSeries", "Expand all"))
        self.label_6.setText(_translate("DialogCorrectionTimeSeries", "<html><head/><body><p>← Select item in the tree widget to show details!</p></body></html>"))
        self.label.setToolTip(_translate("DialogCorrectionTimeSeries", "Survey name."))
        self.label.setText(_translate("DialogCorrectionTimeSeries", "Survey name:"))
        self.label_survey_name.setToolTip(_translate("DialogCorrectionTimeSeries", "Survey name."))
        self.label_2.setToolTip(_translate("DialogCorrectionTimeSeries", "Number of stations in the survey."))
        self.label_2.setText(_translate("DialogCorrectionTimeSeries", "Number of stations:"))
        self.label_num_of_stations.setToolTip(_translate("DialogCorrectionTimeSeries", "Number of stations in the survey."))
        self.label_3.setToolTip(_translate("DialogCorrectionTimeSeries", "Data source, e.g. file name."))
        self.label_3.setText(_translate("DialogCorrectionTimeSeries", "Data source:"))
        self.label_station_data_source.setToolTip(_translate("DialogCorrectionTimeSeries", "Data source, e.g. file name."))
        self.label_4.setToolTip(_translate("DialogCorrectionTimeSeries", "Description of the time series data."))
        self.label_4.setText(_translate("DialogCorrectionTimeSeries", "Description:"))
        self.label_station_description.setToolTip(_translate("DialogCorrectionTimeSeries", "Description of the time series data."))
        self.label_5.setToolTip(_translate("DialogCorrectionTimeSeries", "Unit."))
        self.label_5.setText(_translate("DialogCorrectionTimeSeries", "Unit:"))
        self.label_station_unit.setToolTip(_translate("DialogCorrectionTimeSeries", "Unit."))
        self.label_10.setToolTip(_translate("DialogCorrectionTimeSeries", "<html><head/><body><p>Indicates wheter the time series data represent the <span style=\" font-style:italic;\">effect</span> on the observed gravity, or <span style=\" font-style:italic;\">corrections</span> for the observed gravity (opposite signs!). </p></body></html>"))
        self.label_10.setText(_translate("DialogCorrectionTimeSeries", "Model type:"))
        self.label_station_model_type.setToolTip(_translate("DialogCorrectionTimeSeries", "<html><head/><body><p>Indicates wheter the time series data represent the <span style=\" font-style:italic;\">effect</span> on the observed gravity, or <span style=\" font-style:italic;\">corrections</span> for the observed gravity (opposite signs!). </p></body></html>"))
        self.label_7.setToolTip(_translate("DialogCorrectionTimeSeries", "Start date and time."))
        self.label_7.setText(_translate("DialogCorrectionTimeSeries", "Start:"))
        self.label_station_start.setToolTip(_translate("DialogCorrectionTimeSeries", "Start date and time."))
        self.label_8.setToolTip(_translate("DialogCorrectionTimeSeries", "<html><head/><body><p>End date and time.</p></body></html>"))
        self.label_8.setText(_translate("DialogCorrectionTimeSeries", "End:"))
        self.label_station_end.setToolTip(_translate("DialogCorrectionTimeSeries", "End date and time."))
        self.label_9.setToolTip(_translate("DialogCorrectionTimeSeries", "Duration"))
        self.label_9.setText(_translate("DialogCorrectionTimeSeries", "Duration:"))
        self.label_station_duration.setToolTip(_translate("DialogCorrectionTimeSeries", "Duration"))
        self.label_12.setToolTip(_translate("DialogCorrectionTimeSeries", "Total number of data points."))
        self.label_12.setText(_translate("DialogCorrectionTimeSeries", "No. datapoints:"))
        self.label_station_nuber_of_datapoints.setToolTip(_translate("DialogCorrectionTimeSeries", "Total number of data points."))
        self.label_11.setToolTip(_translate("DialogCorrectionTimeSeries", "Station name."))
        self.label_11.setText(_translate("DialogCorrectionTimeSeries", "Station name:"))
        self.label_station_station_name.setToolTip(_translate("DialogCorrectionTimeSeries", "Station name."))
        self.label_13.setToolTip(_translate("DialogCorrectionTimeSeries", "Date and time (UTC) when the time series data was loaded to GravTools."))
        self.label_13.setText(_translate("DialogCorrectionTimeSeries", "Created (UTC):"))
        self.label_station_created.setToolTip(_translate("DialogCorrectionTimeSeries", "Date and time (UTC) when the time series data was loaded to GravTools."))
        self.graphicsLayoutWidget_time_series.setToolTip(_translate("DialogCorrectionTimeSeries", "Plot of the time series."))
        self.pushButton_load_tsf_file.setToolTip(_translate("DialogCorrectionTimeSeries", "Load TSF file (e.g. createb by TSoft)"))
        self.pushButton_load_tsf_file.setText(_translate("DialogCorrectionTimeSeries", "Load TSF file"))
from pyqtgraph import GraphicsLayoutWidget
