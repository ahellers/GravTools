# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gravtools/gui/dialog_export_results.ui'
#
# Created by: PyQt5 UI code generator 5.15.6
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog_export_results(object):
    def setupUi(self, Dialog_export_results):
        Dialog_export_results.setObjectName("Dialog_export_results")
        Dialog_export_results.resize(453, 745)
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(Dialog_export_results)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.groupBox_general_settings = QtWidgets.QGroupBox(Dialog_export_results)
        self.groupBox_general_settings.setObjectName("groupBox_general_settings")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.groupBox_general_settings)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label_select_lsm_run = QtWidgets.QLabel(self.groupBox_general_settings)
        self.label_select_lsm_run.setObjectName("label_select_lsm_run")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_select_lsm_run)
        self.comboBox_select_lsm_run = QtWidgets.QComboBox(self.groupBox_general_settings)
        self.comboBox_select_lsm_run.setObjectName("comboBox_select_lsm_run")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.comboBox_select_lsm_run)
        self.label_export_path = QtWidgets.QLabel(self.groupBox_general_settings)
        self.label_export_path.setObjectName("label_export_path")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_export_path)
        self.label_export_comment = QtWidgets.QLabel(self.groupBox_general_settings)
        self.label_export_comment.setObjectName("label_export_comment")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_export_comment)
        self.label_export_path_show = QtWidgets.QLabel(self.groupBox_general_settings)
        self.label_export_path_show.setText("")
        self.label_export_path_show.setObjectName("label_export_path_show")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.label_export_path_show)
        self.label_export_comment_show = QtWidgets.QLabel(self.groupBox_general_settings)
        self.label_export_comment_show.setText("")
        self.label_export_comment_show.setObjectName("label_export_comment_show")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.label_export_comment_show)
        self.verticalLayout_5.addLayout(self.formLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(30, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.checkBox_add_lsm_comment_to_filename = QtWidgets.QCheckBox(self.groupBox_general_settings)
        self.checkBox_add_lsm_comment_to_filename.setChecked(True)
        self.checkBox_add_lsm_comment_to_filename.setObjectName("checkBox_add_lsm_comment_to_filename")
        self.horizontalLayout_2.addWidget(self.checkBox_add_lsm_comment_to_filename)
        self.verticalLayout_5.addLayout(self.horizontalLayout_2)
        self.verticalLayout_4.addWidget(self.groupBox_general_settings)
        self.groupBox_nsb_file = QtWidgets.QGroupBox(Dialog_export_results)
        self.groupBox_nsb_file.setObjectName("groupBox_nsb_file")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_nsb_file)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.checkBox_write_nsb_file = QtWidgets.QCheckBox(self.groupBox_nsb_file)
        self.checkBox_write_nsb_file.setChecked(True)
        self.checkBox_write_nsb_file.setObjectName("checkBox_write_nsb_file")
        self.verticalLayout_3.addWidget(self.checkBox_write_nsb_file)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem1 = QtWidgets.QSpacerItem(30, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.formLayout_3 = QtWidgets.QFormLayout()
        self.formLayout_3.setObjectName("formLayout_3")
        self.radioButton_first_dhb_dhf = QtWidgets.QRadioButton(self.groupBox_nsb_file)
        self.radioButton_first_dhb_dhf.setChecked(False)
        self.radioButton_first_dhb_dhf.setObjectName("radioButton_first_dhb_dhf")
        self.buttonGroup_dhb_dhf = QtWidgets.QButtonGroup(Dialog_export_results)
        self.buttonGroup_dhb_dhf.setObjectName("buttonGroup_dhb_dhf")
        self.buttonGroup_dhb_dhf.addButton(self.radioButton_first_dhb_dhf)
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.radioButton_first_dhb_dhf)
        self.radioButton_mean_dhb_dhf = QtWidgets.QRadioButton(self.groupBox_nsb_file)
        self.radioButton_mean_dhb_dhf.setChecked(True)
        self.radioButton_mean_dhb_dhf.setObjectName("radioButton_mean_dhb_dhf")
        self.buttonGroup_dhb_dhf.addButton(self.radioButton_mean_dhb_dhf)
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.radioButton_mean_dhb_dhf)
        self.checkBox_nsb_remove_datum_stations = QtWidgets.QCheckBox(self.groupBox_nsb_file)
        self.checkBox_nsb_remove_datum_stations.setChecked(True)
        self.checkBox_nsb_remove_datum_stations.setObjectName("checkBox_nsb_remove_datum_stations")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.checkBox_nsb_remove_datum_stations)
        self.radioButton_export_sd = QtWidgets.QRadioButton(self.groupBox_nsb_file)
        self.radioButton_export_sd.setObjectName("radioButton_export_sd")
        self.buttonGroup_sd_se = QtWidgets.QButtonGroup(Dialog_export_results)
        self.buttonGroup_sd_se.setObjectName("buttonGroup_sd_se")
        self.buttonGroup_sd_se.addButton(self.radioButton_export_sd)
        self.formLayout_3.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.radioButton_export_sd)
        self.radioButton_export_se = QtWidgets.QRadioButton(self.groupBox_nsb_file)
        self.radioButton_export_se.setChecked(True)
        self.radioButton_export_se.setObjectName("radioButton_export_se")
        self.buttonGroup_sd_se.addButton(self.radioButton_export_se)
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.radioButton_export_se)
        self.horizontalLayout.addLayout(self.formLayout_3)
        self.verticalLayout_3.addLayout(self.horizontalLayout)
        self.verticalLayout_4.addWidget(self.groupBox_nsb_file)
        self.groupBox_other_files = QtWidgets.QGroupBox(Dialog_export_results)
        self.groupBox_other_files.setObjectName("groupBox_other_files")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_other_files)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.checkBox_write_log_file = QtWidgets.QCheckBox(self.groupBox_other_files)
        self.checkBox_write_log_file.setToolTip("")
        self.checkBox_write_log_file.setChecked(True)
        self.checkBox_write_log_file.setObjectName("checkBox_write_log_file")
        self.verticalLayout_2.addWidget(self.checkBox_write_log_file)
        self.checkBox_save_drift_plot_png = QtWidgets.QCheckBox(self.groupBox_other_files)
        self.checkBox_save_drift_plot_png.setEnabled(True)
        self.checkBox_save_drift_plot_png.setChecked(True)
        self.checkBox_save_drift_plot_png.setObjectName("checkBox_save_drift_plot_png")
        self.verticalLayout_2.addWidget(self.checkBox_save_drift_plot_png)
        self.checkBox_save_vg_plot_png = QtWidgets.QCheckBox(self.groupBox_other_files)
        self.checkBox_save_vg_plot_png.setChecked(True)
        self.checkBox_save_vg_plot_png.setObjectName("checkBox_save_vg_plot_png")
        self.verticalLayout_2.addWidget(self.checkBox_save_vg_plot_png)
        self.verticalLayout_4.addWidget(self.groupBox_other_files)
        self.groupBox_observation_list = QtWidgets.QGroupBox(Dialog_export_results)
        self.groupBox_observation_list.setObjectName("groupBox_observation_list")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_observation_list)
        self.verticalLayout.setObjectName("verticalLayout")
        self.checkBox_write_observation_list = QtWidgets.QCheckBox(self.groupBox_observation_list)
        self.checkBox_write_observation_list.setChecked(True)
        self.checkBox_write_observation_list.setObjectName("checkBox_write_observation_list")
        self.verticalLayout.addWidget(self.checkBox_write_observation_list)
        self.formLayout_5 = QtWidgets.QFormLayout()
        self.formLayout_5.setObjectName("formLayout_5")
        spacerItem2 = QtWidgets.QSpacerItem(30, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.formLayout_5.setItem(0, QtWidgets.QFormLayout.LabelRole, spacerItem2)
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.label_observation_list_export_options = QtWidgets.QLabel(self.groupBox_observation_list)
        self.label_observation_list_export_options.setObjectName("label_observation_list_export_options")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_observation_list_export_options)
        self.comboBox_observation_list_export_options = QtWidgets.QComboBox(self.groupBox_observation_list)
        self.comboBox_observation_list_export_options.setObjectName("comboBox_observation_list_export_options")
        self.comboBox_observation_list_export_options.addItem("")
        self.comboBox_observation_list_export_options.addItem("")
        self.comboBox_observation_list_export_options.addItem("")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.comboBox_observation_list_export_options)
        self.formLayout_5.setLayout(0, QtWidgets.QFormLayout.FieldRole, self.formLayout_2)
        self.verticalLayout.addLayout(self.formLayout_5)
        self.verticalLayout_4.addWidget(self.groupBox_observation_list)
        self.groupBox_gis = QtWidgets.QGroupBox(Dialog_export_results)
        self.groupBox_gis.setObjectName("groupBox_gis")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.groupBox_gis)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.checkBox_gis_write_shapefile = QtWidgets.QCheckBox(self.groupBox_gis)
        self.checkBox_gis_write_shapefile.setObjectName("checkBox_gis_write_shapefile")
        self.verticalLayout_6.addWidget(self.checkBox_gis_write_shapefile)
        self.verticalLayout_4.addWidget(self.groupBox_gis)
        spacerItem3 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem3)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog_export_results)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_4.addWidget(self.buttonBox)

        self.retranslateUi(Dialog_export_results)
        self.buttonBox.accepted.connect(Dialog_export_results.accept) # type: ignore
        self.buttonBox.rejected.connect(Dialog_export_results.reject) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(Dialog_export_results)

    def retranslateUi(self, Dialog_export_results):
        _translate = QtCore.QCoreApplication.translate
        Dialog_export_results.setWindowTitle(_translate("Dialog_export_results", "Export data"))
        self.groupBox_general_settings.setToolTip(_translate("Dialog_export_results", "General export settings."))
        self.groupBox_general_settings.setTitle(_translate("Dialog_export_results", "General settings"))
        self.label_select_lsm_run.setText(_translate("Dialog_export_results", "Select LSM run"))
        self.comboBox_select_lsm_run.setToolTip(_translate("Dialog_export_results", "Select one LSM run for data export. If no LSM run is selected, just the observation list csv file can be exported, containing the current observation selection drom the observations tab."))
        self.label_export_path.setToolTip(_translate("Dialog_export_results", "Export path for all files. The path is defined when creating a new campaign and can be changed via the menu \"File/Change output directory\"."))
        self.label_export_path.setText(_translate("Dialog_export_results", "Export_path:"))
        self.label_export_comment.setToolTip(_translate("Dialog_export_results", "Comment of the selected LSM run."))
        self.label_export_comment.setText(_translate("Dialog_export_results", "Comment:"))
        self.checkBox_add_lsm_comment_to_filename.setToolTip(_translate("Dialog_export_results", "Append the LSM run comment to the name of the output file(s), e.g. <campaign name>_<lsm run comment>.nsb"))
        self.checkBox_add_lsm_comment_to_filename.setText(_translate("Dialog_export_results", "Append comment to filenames"))
        self.groupBox_nsb_file.setToolTip(_translate("Dialog_export_results", "Export options for nsb files that contain the estimation results (gravity at stations) of a network adjustment. This option is only available if station gravity was estimated (e.g. not for estimation of gravity gradients)"))
        self.groupBox_nsb_file.setTitle(_translate("Dialog_export_results", "nsb file"))
        self.checkBox_write_nsb_file.setToolTip(_translate("Dialog_export_results", "Write a nsb file."))
        self.checkBox_write_nsb_file.setText(_translate("Dialog_export_results", "Write nsb file"))
        self.radioButton_first_dhb_dhf.setToolTip(_translate("Dialog_export_results", "Take vertical offsets between instrument top and ground (dhb) and instrument top and reference point (dhf) at a station from the first setup at this station."))
        self.radioButton_first_dhb_dhf.setText(_translate("Dialog_export_results", "Take dhb and dhf from first setup at a station"))
        self.radioButton_mean_dhb_dhf.setToolTip(_translate("Dialog_export_results", "Calculate mean values of the vertical offsets between instrument top and ground (dhb) and instrument top and reference point (dhf) from all setups at a station. Write these mean values to the nsb file."))
        self.radioButton_mean_dhb_dhf.setText(_translate("Dialog_export_results", "Use mean dhb and dhf at stations"))
        self.checkBox_nsb_remove_datum_stations.setToolTip(_translate("Dialog_export_results", "Remove all datum stations rom the nsb file."))
        self.checkBox_nsb_remove_datum_stations.setText(_translate("Dialog_export_results", "Remove datum stations from nsb file"))
        self.radioButton_export_sd.setToolTip(_translate("Dialog_export_results", "<html><head/><body><p>Export the standard deviations (SD) of the gravity estimates.</p></body></html>"))
        self.radioButton_export_sd.setText(_translate("Dialog_export_results", "Export std. deviation (SD)"))
        self.radioButton_export_se.setToolTip(_translate("Dialog_export_results", "<html><head/><body><p>Export the standard errors (SE) of the gravity estimates.</p></body></html>"))
        self.radioButton_export_se.setText(_translate("Dialog_export_results", "Export std. error (SE)"))
        self.groupBox_other_files.setTitle(_translate("Dialog_export_results", "Other files"))
        self.checkBox_write_log_file.setText(_translate("Dialog_export_results", "Write log file (text file)"))
        self.checkBox_save_drift_plot_png.setToolTip(_translate("Dialog_export_results", "Save the current drift plot to a PNG file."))
        self.checkBox_save_drift_plot_png.setText(_translate("Dialog_export_results", "Save current drift plot view (PNG file)"))
        self.checkBox_save_vg_plot_png.setToolTip(_translate("Dialog_export_results", "Save the current VG plot to a PNG file (only available if the VG was estimated)."))
        self.checkBox_save_vg_plot_png.setText(_translate("Dialog_export_results", "Save the current VG plot view (PNG file)"))
        self.groupBox_observation_list.setToolTip(_translate("Dialog_export_results", "The observation list file contains either a list of all observations, or of all active or inactive observations. If a LSM run is selected, the observation list that was used to calculate the setup observations for this particular run will be exported. If no LSM run is selected, the current observation list that is visualized in the observations tab (obs. plot and obs. table) will be exported. Observation list files (csv) can be loaded by GravTools at any time to restore a particular selection of observations. "))
        self.groupBox_observation_list.setTitle(_translate("Dialog_export_results", "Observation list (CSV file)"))
        self.checkBox_write_observation_list.setToolTip(_translate("Dialog_export_results", "Write an observation list file."))
        self.checkBox_write_observation_list.setText(_translate("Dialog_export_results", "Write Observation list"))
        self.label_observation_list_export_options.setToolTip(_translate("Dialog_export_results", "Select which observations will be written to the file depending in the activity status."))
        self.label_observation_list_export_options.setText(_translate("Dialog_export_results", "Export "))
        self.comboBox_observation_list_export_options.setToolTip(_translate("Dialog_export_results", "Select which observations will be written to the file depending in the activity status."))
        self.comboBox_observation_list_export_options.setItemText(0, _translate("Dialog_export_results", "all observations"))
        self.comboBox_observation_list_export_options.setItemText(1, _translate("Dialog_export_results", "only active observations"))
        self.comboBox_observation_list_export_options.setItemText(2, _translate("Dialog_export_results", "only inactive observations"))
        self.groupBox_gis.setTitle(_translate("Dialog_export_results", "GIS data"))
        self.checkBox_gis_write_shapefile.setToolTip(_translate("Dialog_export_results", "<html><head/><body><p>Write the stations and observation results to an ESRI shapefile (content of the tabs observation and station results table). The coordinate reference system (crs) of the station coordinates has to be specified by an EPSG code!</p></body></html>"))
        self.checkBox_gis_write_shapefile.setText(_translate("Dialog_export_results", "Write station and observation results to a shapefile"))
