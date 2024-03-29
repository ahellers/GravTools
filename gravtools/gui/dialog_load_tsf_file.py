# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gravtools/gui/dialog_load_tsf_file.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DialogLoadTsfFile(object):
    def setupUi(self, DialogLoadTsfFile):
        DialogLoadTsfFile.setObjectName("DialogLoadTsfFile")
        DialogLoadTsfFile.resize(623, 385)
        DialogLoadTsfFile.setToolTip("")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(DialogLoadTsfFile)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.pushButton_select_file = QtWidgets.QPushButton(DialogLoadTsfFile)
        self.pushButton_select_file.setMaximumSize(QtCore.QSize(100, 16777215))
        self.pushButton_select_file.setObjectName("pushButton_select_file")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.pushButton_select_file)
        self.lineEdit_filename = QtWidgets.QLineEdit(DialogLoadTsfFile)
        self.lineEdit_filename.setObjectName("lineEdit_filename")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.lineEdit_filename)
        self.label_survey_name = QtWidgets.QLabel(DialogLoadTsfFile)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_survey_name.sizePolicy().hasHeightForWidth())
        self.label_survey_name.setSizePolicy(sizePolicy)
        self.label_survey_name.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_survey_name.setObjectName("label_survey_name")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_survey_name)
        self.lineEdit_survey_name = QtWidgets.QLineEdit(DialogLoadTsfFile)
        self.lineEdit_survey_name.setObjectName("lineEdit_survey_name")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.lineEdit_survey_name)
        self.verticalLayout_3.addLayout(self.formLayout_2)
        self.groupBox_filter = QtWidgets.QGroupBox(DialogLoadTsfFile)
        self.groupBox_filter.setObjectName("groupBox_filter")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_filter)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label_2 = QtWidgets.QLabel(self.groupBox_filter)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.lineEdit_filter_location = QtWidgets.QLineEdit(self.groupBox_filter)
        self.lineEdit_filter_location.setObjectName("lineEdit_filter_location")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.lineEdit_filter_location)
        self.label_3 = QtWidgets.QLabel(self.groupBox_filter)
        self.label_3.setObjectName("label_3")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.lineEdit_filter_instrument = QtWidgets.QLineEdit(self.groupBox_filter)
        self.lineEdit_filter_instrument.setObjectName("lineEdit_filter_instrument")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.lineEdit_filter_instrument)
        self.label_4 = QtWidgets.QLabel(self.groupBox_filter)
        self.label_4.setObjectName("label_4")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_4)
        self.lineEdit_filter_data_type = QtWidgets.QLineEdit(self.groupBox_filter)
        self.lineEdit_filter_data_type.setObjectName("lineEdit_filter_data_type")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.lineEdit_filter_data_type)
        self.verticalLayout_2.addLayout(self.formLayout)
        self.verticalLayout_3.addWidget(self.groupBox_filter)
        self.groupBox_options = QtWidgets.QGroupBox(DialogLoadTsfFile)
        self.groupBox_options.setObjectName("groupBox_options")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_options)
        self.verticalLayout.setObjectName("verticalLayout")
        self.radioButton_effect = QtWidgets.QRadioButton(self.groupBox_options)
        self.radioButton_effect.setChecked(True)
        self.radioButton_effect.setObjectName("radioButton_effect")
        self.buttonGroup_effect_or_correction = QtWidgets.QButtonGroup(DialogLoadTsfFile)
        self.buttonGroup_effect_or_correction.setObjectName("buttonGroup_effect_or_correction")
        self.buttonGroup_effect_or_correction.addButton(self.radioButton_effect)
        self.verticalLayout.addWidget(self.radioButton_effect)
        self.radioButton_correction = QtWidgets.QRadioButton(self.groupBox_options)
        self.radioButton_correction.setObjectName("radioButton_correction")
        self.buttonGroup_effect_or_correction.addButton(self.radioButton_correction)
        self.verticalLayout.addWidget(self.radioButton_correction)
        self.checkBox_overwrite_channel = QtWidgets.QCheckBox(self.groupBox_options)
        self.checkBox_overwrite_channel.setChecked(True)
        self.checkBox_overwrite_channel.setObjectName("checkBox_overwrite_channel")
        self.verticalLayout.addWidget(self.checkBox_overwrite_channel)
        self.verticalLayout_3.addWidget(self.groupBox_options)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem)
        self.buttonBox = QtWidgets.QDialogButtonBox(DialogLoadTsfFile)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_3.addWidget(self.buttonBox)

        self.retranslateUi(DialogLoadTsfFile)
        self.buttonBox.accepted.connect(DialogLoadTsfFile.accept) # type: ignore
        self.buttonBox.rejected.connect(DialogLoadTsfFile.reject) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(DialogLoadTsfFile)

    def retranslateUi(self, DialogLoadTsfFile):
        _translate = QtCore.QCoreApplication.translate
        DialogLoadTsfFile.setWindowTitle(_translate("DialogLoadTsfFile", "Load TSF file"))
        self.pushButton_select_file.setToolTip(_translate("DialogLoadTsfFile", "Select TSF file"))
        self.pushButton_select_file.setText(_translate("DialogLoadTsfFile", "Select TSF file"))
        self.lineEdit_filename.setToolTip(_translate("DialogLoadTsfFile", "Filename and path of the TSF file"))
        self.label_survey_name.setToolTip(_translate("DialogLoadTsfFile", "Name of the survey for which the time series data is loaded from the TSF file"))
        self.label_survey_name.setText(_translate("DialogLoadTsfFile", "Survey name"))
        self.lineEdit_survey_name.setToolTip(_translate("DialogLoadTsfFile", "Name of the survey for which the time series data is loaded from the TSF file"))
        self.groupBox_filter.setToolTip(_translate("DialogLoadTsfFile", "<html><head/><body><p>Filter the channels that are loaded from the TSF file by defining restrictions for the location, the instrumt and/or the data type. The filters are applied on the channel names in the [CHANNELS] block that have to follow the convention <span style=\" font-style:italic;\">&lt;location&gt;:&lt;instrument&gt;:&lt;data type&gt;</span>.</p></body></html>"))
        self.groupBox_filter.setTitle(_translate("DialogLoadTsfFile", "Filter"))
        self.label_2.setToolTip(_translate("DialogLoadTsfFile", "Filter the channels that are loaded from the TSF file by defining restrictions for the location, the instrumt and/or the data type. The filters are applied on the channel names in the [CHANNELS] block that have to follow the convention <location>:<instrument>:<data type>."))
        self.label_2.setText(_translate("DialogLoadTsfFile", "Location:"))
        self.lineEdit_filter_location.setToolTip(_translate("DialogLoadTsfFile", "Filter the channels that are loaded from the TSF file by defining restrictions for the location, the instrumt and/or the data type. The filters are applied on the channel names in the [CHANNELS] block that have to follow the convention <location>:<instrument>:<data type>."))
        self.label_3.setToolTip(_translate("DialogLoadTsfFile", "Filter the channels that are loaded from the TSF file by defining restrictions for the location, the instrumt and/or the data type. The filters are applied on the channel names in the [CHANNELS] block that have to follow the convention <location>:<instrument>:<data type>."))
        self.label_3.setText(_translate("DialogLoadTsfFile", "Instrument:"))
        self.lineEdit_filter_instrument.setToolTip(_translate("DialogLoadTsfFile", "Filter the channels that are loaded from the TSF file by defining restrictions for the location, the instrumt and/or the data type. The filters are applied on the channel names in the [CHANNELS] block that have to follow the convention <location>:<instrument>:<data type>."))
        self.label_4.setToolTip(_translate("DialogLoadTsfFile", "Filter the channels that are loaded from the TSF file by defining restrictions for the location, the instrumt and/or the data type. The filters are applied on the channel names in the [CHANNELS] block that have to follow the convention <location>:<instrument>:<data type>."))
        self.label_4.setText(_translate("DialogLoadTsfFile", "Data type:"))
        self.lineEdit_filter_data_type.setToolTip(_translate("DialogLoadTsfFile", "Filter the channels that are loaded from the TSF file by defining restrictions for the location, the instrumt and/or the data type. The filters are applied on the channel names in the [CHANNELS] block that have to follow the convention <location>:<instrument>:<data type>."))
        self.groupBox_options.setToolTip(_translate("DialogLoadTsfFile", "General options"))
        self.groupBox_options.setTitle(_translate("DialogLoadTsfFile", "Options"))
        self.radioButton_effect.setToolTip(_translate("DialogLoadTsfFile", "<html><head/><body><p>The channel data represents the gravity effect of phenomena, e.g. the gravitational attraction of sun and moon. Effects are <span style=\" font-weight:600;\">subtracted</span> from gravity observations in order to correct them (default for synthetic gravity calculated in TSoft).</p><p>See <span style=\" font-style:italic;\">TSoft manual (version 2.2.4 Release date 2015-09-09), p. 15</span>: &quot;<span style=\" font-style:italic;\">The loading calculation by Tsoft will result in the correction, not the effect. On the other hand the prediction of the solid Earth tides using the WDD parameter set provided by Tsoft directly will result in the effect. Therefore both must be treated with different sign in order to reduce a gravity time series correctly.</span>&quot;</p></body></html>"))
        self.radioButton_effect.setText(_translate("DialogLoadTsfFile", "Time series describe the gravity effect of phenomena"))
        self.radioButton_correction.setToolTip(_translate("DialogLoadTsfFile", "<html><head/><body><p>The channel data represents corrections for gravity observations, i.e. these values are <span style=\" font-weight:600;\">added</span> to gravity observations in order to correct them.</p><p>See <span style=\" font-style:italic;\">TSoft manual (version 2.2.4 Release date 2015-09-09), p. 15</span>: &quot;<span style=\" font-style:italic;\">The loading calculation by Tsoft will result in the correction, not the effect. On the other hand the prediction of the solid Earth tides using the WDD parameter set provided by Tsoft directly will result in the effect. Therefore both must be treated with different sign in order to reduce a gravity time series correctly.</span>&quot;</p></body></html>"))
        self.radioButton_correction.setText(_translate("DialogLoadTsfFile", "Time series represent corrections for gravity observations"))
        self.checkBox_overwrite_channel.setToolTip(_translate("DialogLoadTsfFile", "Overwrite existing channels with the same location (= station name) when loading a TSF file."))
        self.checkBox_overwrite_channel.setText(_translate("DialogLoadTsfFile", "Overwrite existing channels with identical locations"))
