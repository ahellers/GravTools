# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gravtools/gui/dialog_calc_drift.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog_calculate_drift(object):
    def setupUi(self, Dialog_calculate_drift):
        Dialog_calculate_drift.setObjectName("Dialog_calculate_drift")
        Dialog_calculate_drift.resize(1280, 697)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog_calculate_drift)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox_data_selection = QtWidgets.QGroupBox(Dialog_calculate_drift)
        self.groupBox_data_selection.setObjectName("groupBox_data_selection")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox_data_selection)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label_survey = QtWidgets.QLabel(self.groupBox_data_selection)
        self.label_survey.setObjectName("label_survey")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_survey)
        self.comboBox_survey = QtWidgets.QComboBox(self.groupBox_data_selection)
        self.comboBox_survey.setObjectName("comboBox_survey")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.comboBox_survey)
        self.label_station = QtWidgets.QLabel(self.groupBox_data_selection)
        self.label_station.setObjectName("label_station")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_station)
        self.comboBox_station = QtWidgets.QComboBox(self.groupBox_data_selection)
        self.comboBox_station.setObjectName("comboBox_station")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.comboBox_station)
        self.verticalLayout_4.addLayout(self.formLayout)
        self.checkBox_multiple_setups_only = QtWidgets.QCheckBox(self.groupBox_data_selection)
        self.checkBox_multiple_setups_only.setChecked(True)
        self.checkBox_multiple_setups_only.setObjectName("checkBox_multiple_setups_only")
        self.verticalLayout_4.addWidget(self.checkBox_multiple_setups_only)
        self.checkBox_active_obs_only = QtWidgets.QCheckBox(self.groupBox_data_selection)
        self.checkBox_active_obs_only.setChecked(True)
        self.checkBox_active_obs_only.setObjectName("checkBox_active_obs_only")
        self.verticalLayout_4.addWidget(self.checkBox_active_obs_only)
        self.verticalLayout.addWidget(self.groupBox_data_selection)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.groupBox_options = QtWidgets.QGroupBox(Dialog_calculate_drift)
        self.groupBox_options.setObjectName("groupBox_options")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_options)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.label_degree = QtWidgets.QLabel(self.groupBox_options)
        self.label_degree.setObjectName("label_degree")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_degree)
        self.spinBox_degree = QtWidgets.QSpinBox(self.groupBox_options)
        self.spinBox_degree.setMinimum(1)
        self.spinBox_degree.setMaximum(3)
        self.spinBox_degree.setObjectName("spinBox_degree")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.spinBox_degree)
        self.label_method = QtWidgets.QLabel(self.groupBox_options)
        self.label_method.setObjectName("label_method")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_method)
        self.comboBox_method = QtWidgets.QComboBox(self.groupBox_options)
        self.comboBox_method.setObjectName("comboBox_method")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.comboBox_method)
        self.label_min_no_obs = QtWidgets.QLabel(self.groupBox_options)
        self.label_min_no_obs.setObjectName("label_min_no_obs")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_min_no_obs)
        self.spinBox_min_no_obs = QtWidgets.QSpinBox(self.groupBox_options)
        self.spinBox_min_no_obs.setMinimum(1)
        self.spinBox_min_no_obs.setProperty("value", 2)
        self.spinBox_min_no_obs.setObjectName("spinBox_min_no_obs")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.spinBox_min_no_obs)
        self.verticalLayout_3.addLayout(self.formLayout_2)
        self.checkBox_remove_g_offset = QtWidgets.QCheckBox(self.groupBox_options)
        self.checkBox_remove_g_offset.setObjectName("checkBox_remove_g_offset")
        self.verticalLayout_3.addWidget(self.checkBox_remove_g_offset)
        self.horizontalLayout.addWidget(self.groupBox_options)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.graphicsLayoutWidget_drift_plot = GraphicsLayoutWidget(Dialog_calculate_drift)
        self.graphicsLayoutWidget_drift_plot.setObjectName("graphicsLayoutWidget_drift_plot")
        self.verticalLayout_2.addWidget(self.graphicsLayoutWidget_drift_plot)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog_calculate_drift)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Save)
        self.buttonBox.setCenterButtons(False)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_2.addWidget(self.buttonBox)

        self.retranslateUi(Dialog_calculate_drift)
        self.buttonBox.accepted.connect(Dialog_calculate_drift.accept) # type: ignore
        self.buttonBox.rejected.connect(Dialog_calculate_drift.reject) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(Dialog_calculate_drift)

    def retranslateUi(self, Dialog_calculate_drift):
        _translate = QtCore.QCoreApplication.translate
        Dialog_calculate_drift.setWindowTitle(_translate("Dialog_calculate_drift", "Calculate drift"))
        self.groupBox_data_selection.setToolTip(_translate("Dialog_calculate_drift", "Options for selection of observation data that is used for calculating the drift."))
        self.groupBox_data_selection.setTitle(_translate("Dialog_calculate_drift", "Data selection"))
        self.label_survey.setToolTip(_translate("Dialog_calculate_drift", "Select a survey."))
        self.label_survey.setText(_translate("Dialog_calculate_drift", "Survey"))
        self.comboBox_survey.setToolTip(_translate("Dialog_calculate_drift", "Select a survey."))
        self.label_station.setToolTip(_translate("Dialog_calculate_drift", "Either select all stations or one station for drift determination."))
        self.label_station.setText(_translate("Dialog_calculate_drift", "Station"))
        self.comboBox_station.setToolTip(_translate("Dialog_calculate_drift", "Either select all stations or one station for drift determination."))
        self.checkBox_multiple_setups_only.setToolTip(_translate("Dialog_calculate_drift", "Only consider stations with multiple instrument setups."))
        self.checkBox_multiple_setups_only.setText(_translate("Dialog_calculate_drift", "Stations with multiple setups only"))
        self.checkBox_active_obs_only.setToolTip(_translate("Dialog_calculate_drift", "Only consider active observations (based on the selection in the observations tab in the main window)."))
        self.checkBox_active_obs_only.setText(_translate("Dialog_calculate_drift", "Only consider active observations"))
        self.groupBox_options.setToolTip(_translate("Dialog_calculate_drift", "<html><head/><body><p>Options for the drift calculation and the visualization in the plotwindows below. The drift is modelled as time dependent polynomial of the specified degree.</p></body></html>"))
        self.groupBox_options.setTitle(_translate("Dialog_calculate_drift", "Options"))
        self.label_degree.setToolTip(_translate("Dialog_calculate_drift", "Degree of the drift polynomial."))
        self.label_degree.setText(_translate("Dialog_calculate_drift", "Poly. degree"))
        self.spinBox_degree.setToolTip(_translate("Dialog_calculate_drift", "Degree of the drift polynomial."))
        self.label_method.setToolTip(_translate("Dialog_calculate_drift", "Drift calculation method."))
        self.label_method.setText(_translate("Dialog_calculate_drift", "Drift calc. method"))
        self.comboBox_method.setToolTip(_translate("Dialog_calculate_drift", "Drift calculation method."))
        self.label_min_no_obs.setToolTip(_translate("Dialog_calculate_drift", "<html><head/><body><p>Minimum number of observations required to fit a drift polynomial.</p></body></html>"))
        self.label_min_no_obs.setText(_translate("Dialog_calculate_drift", "Min. n/o observations"))
        self.spinBox_min_no_obs.setToolTip(_translate("Dialog_calculate_drift", "<html><head/><body><p>Minimum number of observations required to fit a drift polynomial.</p></body></html>"))
        self.checkBox_remove_g_offset.setToolTip(_translate("Dialog_calculate_drift", "<html><head/><body><p>Remove gravity offsets between the drift functions and observations of different stations. The first station in the dropdown list is used as reference and remains unchanged. This allows for depicting drift functions of multiple stations as an overlay.</p></body></html>"))
        self.checkBox_remove_g_offset.setText(_translate("Dialog_calculate_drift", "Remove gravity offsets between stations in the plot"))
        self.graphicsLayoutWidget_drift_plot.setToolTip(_translate("Dialog_calculate_drift", "Drift plot depicting the observerd gravity values at stations and the fitted drift polynomials a functions of time."))
        self.buttonBox.setToolTip(_translate("Dialog_calculate_drift", "<html><head/><body><p>Click <span style=\" font-weight:600; font-style:italic;\">Save</span> to create a drift report file and a drift plot file in the campaign directory. Click <span style=\" font-weight:600; font-style:italic;\">Cancel</span> to close this windows.</p></body></html>"))
from pyqtgraph import GraphicsLayoutWidget