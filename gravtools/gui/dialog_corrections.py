# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gravtools/gui/dialog_corrections.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog_corrections(object):
    def setupUi(self, Dialog_corrections):
        Dialog_corrections.setObjectName("Dialog_corrections")
        Dialog_corrections.resize(274, 311)
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(Dialog_corrections)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.groupBox_corrections_tides = QtWidgets.QGroupBox(Dialog_corrections)
        self.groupBox_corrections_tides.setObjectName("groupBox_corrections_tides")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_corrections_tides)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.radioButton_corr_tides_no_correction = QtWidgets.QRadioButton(self.groupBox_corrections_tides)
        self.radioButton_corr_tides_no_correction.setObjectName("radioButton_corr_tides_no_correction")
        self.buttonGroup_corrections_tides = QtWidgets.QButtonGroup(Dialog_corrections)
        self.buttonGroup_corrections_tides.setObjectName("buttonGroup_corrections_tides")
        self.buttonGroup_corrections_tides.addButton(self.radioButton_corr_tides_no_correction)
        self.verticalLayout.addWidget(self.radioButton_corr_tides_no_correction)
        self.radioButton_corr_tides_cg5_model = QtWidgets.QRadioButton(self.groupBox_corrections_tides)
        self.radioButton_corr_tides_cg5_model.setChecked(True)
        self.radioButton_corr_tides_cg5_model.setObjectName("radioButton_corr_tides_cg5_model")
        self.buttonGroup_corrections_tides.addButton(self.radioButton_corr_tides_cg5_model)
        self.verticalLayout.addWidget(self.radioButton_corr_tides_cg5_model)
        self.verticalLayout_3.addLayout(self.verticalLayout)
        self.verticalLayout_5.addWidget(self.groupBox_corrections_tides)
        self.groupBox_corrections_ref_heights = QtWidgets.QGroupBox(Dialog_corrections)
        self.groupBox_corrections_ref_heights.setObjectName("groupBox_corrections_ref_heights")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox_corrections_ref_heights)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.radioButton_corr_ref_heights_sensor = QtWidgets.QRadioButton(self.groupBox_corrections_ref_heights)
        self.radioButton_corr_ref_heights_sensor.setObjectName("radioButton_corr_ref_heights_sensor")
        self.buttonGroup_corrections_ref_heights = QtWidgets.QButtonGroup(Dialog_corrections)
        self.buttonGroup_corrections_ref_heights.setObjectName("buttonGroup_corrections_ref_heights")
        self.buttonGroup_corrections_ref_heights.addButton(self.radioButton_corr_ref_heights_sensor)
        self.verticalLayout_2.addWidget(self.radioButton_corr_ref_heights_sensor)
        self.radioButton_corr_ref_heights_instrument_top = QtWidgets.QRadioButton(self.groupBox_corrections_ref_heights)
        self.radioButton_corr_ref_heights_instrument_top.setObjectName("radioButton_corr_ref_heights_instrument_top")
        self.buttonGroup_corrections_ref_heights.addButton(self.radioButton_corr_ref_heights_instrument_top)
        self.verticalLayout_2.addWidget(self.radioButton_corr_ref_heights_instrument_top)
        self.radioButton_corr_ref_heights_ground = QtWidgets.QRadioButton(self.groupBox_corrections_ref_heights)
        self.radioButton_corr_ref_heights_ground.setObjectName("radioButton_corr_ref_heights_ground")
        self.buttonGroup_corrections_ref_heights.addButton(self.radioButton_corr_ref_heights_ground)
        self.verticalLayout_2.addWidget(self.radioButton_corr_ref_heights_ground)
        self.radioButton_corr_ref_heights_control_point = QtWidgets.QRadioButton(self.groupBox_corrections_ref_heights)
        self.radioButton_corr_ref_heights_control_point.setChecked(True)
        self.radioButton_corr_ref_heights_control_point.setObjectName("radioButton_corr_ref_heights_control_point")
        self.buttonGroup_corrections_ref_heights.addButton(self.radioButton_corr_ref_heights_control_point)
        self.verticalLayout_2.addWidget(self.radioButton_corr_ref_heights_control_point)
        self.verticalLayout_4.addLayout(self.verticalLayout_2)
        self.verticalLayout_5.addWidget(self.groupBox_corrections_ref_heights)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog_corrections)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_5.addWidget(self.buttonBox)
        self.verticalLayout_6.addLayout(self.verticalLayout_5)

        self.retranslateUi(Dialog_corrections)
        self.buttonBox.accepted.connect(Dialog_corrections.accept)
        self.buttonBox.rejected.connect(Dialog_corrections.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog_corrections)

    def retranslateUi(self, Dialog_corrections):
        _translate = QtCore.QCoreApplication.translate
        Dialog_corrections.setWindowTitle(_translate("Dialog_corrections", "Corrections"))
        self.groupBox_corrections_tides.setTitle(_translate("Dialog_corrections", "Tidal correction"))
        self.radioButton_corr_tides_no_correction.setText(_translate("Dialog_corrections", "No correction"))
        self.radioButton_corr_tides_cg5_model.setText(_translate("Dialog_corrections", "CG-5 model (Longman, 1959)"))
        self.groupBox_corrections_ref_heights.setTitle(_translate("Dialog_corrections", "Reference height"))
        self.radioButton_corr_ref_heights_sensor.setText(_translate("Dialog_corrections", "Sensor"))
        self.radioButton_corr_ref_heights_instrument_top.setText(_translate("Dialog_corrections", "Instrument top"))
        self.radioButton_corr_ref_heights_ground.setText(_translate("Dialog_corrections", "Ground"))
        self.radioButton_corr_ref_heights_control_point.setText(_translate("Dialog_corrections", "Control point"))
        self.buttonBox.setToolTip(_translate("Dialog_corrections", "When clicking OK the corrections are applied on all observations in the current campaign."))
