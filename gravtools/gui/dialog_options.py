# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gravtools/gui/dialog_options.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog_options(object):
    def setupUi(self, Dialog_options):
        Dialog_options.setObjectName("Dialog_options")
        Dialog_options.resize(319, 192)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(Dialog_options)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.groupBox_gui_options = QtWidgets.QGroupBox(Dialog_options)
        self.groupBox_gui_options.setObjectName("groupBox_gui_options")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_gui_options)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label = QtWidgets.QLabel(self.groupBox_gui_options)
        self.label.setObjectName("label")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.radioButton_gui_mode_simple = QtWidgets.QRadioButton(self.groupBox_gui_options)
        self.radioButton_gui_mode_simple.setChecked(True)
        self.radioButton_gui_mode_simple.setObjectName("radioButton_gui_mode_simple")
        self.buttonGroup_gui_mode = QtWidgets.QButtonGroup(Dialog_options)
        self.buttonGroup_gui_mode.setObjectName("buttonGroup_gui_mode")
        self.buttonGroup_gui_mode.addButton(self.radioButton_gui_mode_simple)
        self.verticalLayout.addWidget(self.radioButton_gui_mode_simple)
        self.radioButton_gui_mode_advanced = QtWidgets.QRadioButton(self.groupBox_gui_options)
        self.radioButton_gui_mode_advanced.setObjectName("radioButton_gui_mode_advanced")
        self.buttonGroup_gui_mode.addButton(self.radioButton_gui_mode_advanced)
        self.verticalLayout.addWidget(self.radioButton_gui_mode_advanced)
        self.formLayout.setLayout(0, QtWidgets.QFormLayout.FieldRole, self.verticalLayout)
        self.verticalLayout_2.addLayout(self.formLayout)
        self.verticalLayout_3.addWidget(self.groupBox_gui_options)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog_options)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_3.addWidget(self.buttonBox)

        self.retranslateUi(Dialog_options)
        self.buttonBox.accepted.connect(Dialog_options.accept) # type: ignore
        self.buttonBox.rejected.connect(Dialog_options.reject) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(Dialog_options)

    def retranslateUi(self, Dialog_options):
        _translate = QtCore.QCoreApplication.translate
        Dialog_options.setWindowTitle(_translate("Dialog_options", "Options"))
        self.groupBox_gui_options.setTitle(_translate("Dialog_options", "Graphical user interface (GUI) options"))
        self.label.setToolTip(_translate("Dialog_options", "Select the level of details that is preented in the GUI. The advanced mode provides more options."))
        self.label.setText(_translate("Dialog_options", "GUI mode"))
        self.radioButton_gui_mode_simple.setToolTip(_translate("Dialog_options", "Simple GUI appearance."))
        self.radioButton_gui_mode_simple.setText(_translate("Dialog_options", "Simple mode"))
        self.radioButton_gui_mode_advanced.setToolTip(_translate("Dialog_options", "Advanced GUI appearance."))
        self.radioButton_gui_mode_advanced.setText(_translate("Dialog_options", "Advanced mode"))
