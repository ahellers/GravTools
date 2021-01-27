# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gravtools/gui/MainWindow.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1403, 729)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.tabWidget_Main = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget_Main.setTabPosition(QtWidgets.QTabWidget.East)
        self.tabWidget_Main.setMovable(True)
        self.tabWidget_Main.setObjectName("tabWidget_Main")
        self.tab_main_stations = QtWidgets.QWidget()
        self.tab_main_stations.setObjectName("tab_main_stations")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.tab_main_stations)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.groupBox_filter_options = QtWidgets.QGroupBox(self.tab_main_stations)
        self.groupBox_filter_options.setEnabled(False)
        self.groupBox_filter_options.setFlat(False)
        self.groupBox_filter_options.setCheckable(False)
        self.groupBox_filter_options.setObjectName("groupBox_filter_options")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_filter_options)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.checkBox_filter_observed_stat_only = QtWidgets.QCheckBox(self.groupBox_filter_options)
        self.checkBox_filter_observed_stat_only.setObjectName("checkBox_filter_observed_stat_only")
        self.verticalLayout_2.addWidget(self.checkBox_filter_observed_stat_only)
        self.label = QtWidgets.QLabel(self.groupBox_filter_options)
        self.label.setObjectName("label")
        self.verticalLayout_2.addWidget(self.label)
        self.lineEdit_filter_stat_name = QtWidgets.QLineEdit(self.groupBox_filter_options)
        self.lineEdit_filter_stat_name.setText("")
        self.lineEdit_filter_stat_name.setObjectName("lineEdit_filter_stat_name")
        self.verticalLayout_2.addWidget(self.lineEdit_filter_stat_name)
        self.verticalLayout_3.addWidget(self.groupBox_filter_options)
        self.groupBox_edit_options = QtWidgets.QGroupBox(self.tab_main_stations)
        self.groupBox_edit_options.setEnabled(False)
        self.groupBox_edit_options.setObjectName("groupBox_edit_options")
        self.verticalLayout_3.addWidget(self.groupBox_edit_options)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem)
        self.horizontalLayout_2.addLayout(self.verticalLayout_3)
        self.tab_Widget_Stations = QtWidgets.QTabWidget(self.tab_main_stations)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tab_Widget_Stations.sizePolicy().hasHeightForWidth())
        self.tab_Widget_Stations.setSizePolicy(sizePolicy)
        self.tab_Widget_Stations.setObjectName("tab_Widget_Stations")
        self.tab_Stations_Table = QtWidgets.QWidget()
        self.tab_Stations_Table.setObjectName("tab_Stations_Table")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.tab_Stations_Table)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tableView_Stations = QtWidgets.QTableView(self.tab_Stations_Table)
        self.tableView_Stations.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.tableView_Stations.setObjectName("tableView_Stations")
        self.tableView_Stations.horizontalHeader().setCascadingSectionResizes(True)
        self.verticalLayout.addWidget(self.tableView_Stations)
        self.tab_Widget_Stations.addTab(self.tab_Stations_Table, "")
        self.tab_Stations_Plot = QtWidgets.QWidget()
        self.tab_Stations_Plot.setObjectName("tab_Stations_Plot")
        self.tab_Widget_Stations.addTab(self.tab_Stations_Plot, "")
        self.horizontalLayout_2.addWidget(self.tab_Widget_Stations)
        self.horizontalLayout_2.setStretch(0, 1)
        self.horizontalLayout_2.setStretch(1, 4)
        self.verticalLayout_4.addLayout(self.horizontalLayout_2)
        self.tabWidget_Main.addTab(self.tab_main_stations, "")
        self.tab_main_observations = QtWidgets.QWidget()
        self.tab_main_observations.setObjectName("tab_main_observations")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.tab_main_observations)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.splitter_obs = QtWidgets.QSplitter(self.tab_main_observations)
        self.splitter_obs.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_obs.setObjectName("splitter_obs")
        self.layoutWidget = QtWidgets.QWidget(self.splitter_obs)
        self.layoutWidget.setObjectName("layoutWidget")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.treeWidget_observations = QtWidgets.QTreeWidget(self.layoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.treeWidget_observations.sizePolicy().hasHeightForWidth())
        self.treeWidget_observations.setSizePolicy(sizePolicy)
        self.treeWidget_observations.setMinimumSize(QtCore.QSize(280, 0))
        self.treeWidget_observations.setHeaderHidden(True)
        self.treeWidget_observations.setObjectName("treeWidget_observations")
        self.treeWidget_observations.header().setVisible(False)
        self.treeWidget_observations.header().setStretchLastSection(True)
        self.verticalLayout_5.addWidget(self.treeWidget_observations)
        self.groupBox_obs_view_options = QtWidgets.QGroupBox(self.layoutWidget)
        self.groupBox_obs_view_options.setObjectName("groupBox_obs_view_options")
        self.verticalLayout_9 = QtWidgets.QVBoxLayout(self.groupBox_obs_view_options)
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.checkBox_obs_show_all_columns = QtWidgets.QCheckBox(self.groupBox_obs_view_options)
        self.checkBox_obs_show_all_columns.setChecked(True)
        self.checkBox_obs_show_all_columns.setObjectName("checkBox_obs_show_all_columns")
        self.verticalLayout_9.addWidget(self.checkBox_obs_show_all_columns)
        self.checkBox_obs_plot_reduced_observations = QtWidgets.QCheckBox(self.groupBox_obs_view_options)
        self.checkBox_obs_plot_reduced_observations.setObjectName("checkBox_obs_plot_reduced_observations")
        self.verticalLayout_9.addWidget(self.checkBox_obs_plot_reduced_observations)
        self.verticalLayout_5.addWidget(self.groupBox_obs_view_options)
        self.tabWidget_observations = QtWidgets.QTabWidget(self.splitter_obs)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget_observations.sizePolicy().hasHeightForWidth())
        self.tabWidget_observations.setSizePolicy(sizePolicy)
        self.tabWidget_observations.setObjectName("tabWidget_observations")
        self.tab_observations_table = QtWidgets.QWidget()
        self.tab_observations_table.setObjectName("tab_observations_table")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.tab_observations_table)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.tableView_observations = QtWidgets.QTableView(self.tab_observations_table)
        self.tableView_observations.setObjectName("tableView_observations")
        self.verticalLayout_6.addWidget(self.tableView_observations)
        self.tabWidget_observations.addTab(self.tab_observations_table, "")
        self.tab_observations_plots = QtWidgets.QWidget()
        self.tab_observations_plots.setObjectName("tab_observations_plots")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.tab_observations_plots)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.GraphicsLayoutWidget_observations = GraphicsLayoutWidget(self.tab_observations_plots)
        self.GraphicsLayoutWidget_observations.setObjectName("GraphicsLayoutWidget_observations")
        self.verticalLayout_8.addWidget(self.GraphicsLayoutWidget_observations)
        self.tabWidget_observations.addTab(self.tab_observations_plots, "")
        self.verticalLayout_7.addWidget(self.splitter_obs)
        self.tabWidget_Main.addTab(self.tab_main_observations, "")
        self.tab_main_setups = QtWidgets.QWidget()
        self.tab_main_setups.setObjectName("tab_main_setups")
        self.tabWidget_Main.addTab(self.tab_main_setups, "")
        self.horizontalLayout_3.addWidget(self.tabWidget_Main)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1403, 22))
        self.menubar.setObjectName("menubar")
        self.menu_File = QtWidgets.QMenu(self.menubar)
        self.menu_File.setObjectName("menu_File")
        self.menuAdd_Survey = QtWidgets.QMenu(self.menu_File)
        self.menuAdd_Survey.setEnabled(False)
        self.menuAdd_Survey.setObjectName("menuAdd_Survey")
        self.menu_Observations = QtWidgets.QMenu(self.menubar)
        self.menu_Observations.setTearOffEnabled(False)
        self.menu_Observations.setObjectName("menu_Observations")
        self.menuEstimation_settings = QtWidgets.QMenu(self.menubar)
        self.menuEstimation_settings.setObjectName("menuEstimation_settings")
        self.menuStations = QtWidgets.QMenu(self.menubar)
        self.menuStations.setObjectName("menuStations")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setInputMethodHints(QtCore.Qt.ImhNone)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.action_New_Campaign = QtWidgets.QAction(MainWindow)
        self.action_New_Campaign.setObjectName("action_New_Campaign")
        self.action_Save_Campaign = QtWidgets.QAction(MainWindow)
        self.action_Save_Campaign.setEnabled(False)
        self.action_Save_Campaign.setObjectName("action_Save_Campaign")
        self.action_Load_Campaign = QtWidgets.QAction(MainWindow)
        self.action_Load_Campaign.setEnabled(False)
        self.action_Load_Campaign.setObjectName("action_Load_Campaign")
        self.action_Add_Stations = QtWidgets.QAction(MainWindow)
        self.action_Add_Stations.setEnabled(False)
        self.action_Add_Stations.setObjectName("action_Add_Stations")
        self.action_Exit = QtWidgets.QAction(MainWindow)
        self.action_Exit.setObjectName("action_Exit")
        self.action_Edit_Observations = QtWidgets.QAction(MainWindow)
        self.action_Edit_Observations.setObjectName("action_Edit_Observations")
        self.actionResiduals = QtWidgets.QAction(MainWindow)
        self.actionResiduals.setObjectName("actionResiduals")
        self.actionEstimates = QtWidgets.QAction(MainWindow)
        self.actionEstimates.setObjectName("actionEstimates")
        self.action_Corrections = QtWidgets.QAction(MainWindow)
        self.action_Corrections.setObjectName("action_Corrections")
        self.actionShow_Stations = QtWidgets.QAction(MainWindow)
        self.actionShow_Stations.setObjectName("actionShow_Stations")
        self.action_from_CG5_observation_file = QtWidgets.QAction(MainWindow)
        self.action_from_CG5_observation_file.setObjectName("action_from_CG5_observation_file")
        self.action_from_BEV_observation_file = QtWidgets.QAction(MainWindow)
        self.action_from_BEV_observation_file.setObjectName("action_from_BEV_observation_file")
        self.menuAdd_Survey.addAction(self.action_from_CG5_observation_file)
        self.menuAdd_Survey.addAction(self.action_from_BEV_observation_file)
        self.menu_File.addAction(self.action_New_Campaign)
        self.menu_File.addAction(self.action_Save_Campaign)
        self.menu_File.addAction(self.action_Load_Campaign)
        self.menu_File.addSeparator()
        self.menu_File.addAction(self.menuAdd_Survey.menuAction())
        self.menu_File.addAction(self.action_Add_Stations)
        self.menu_File.addSeparator()
        self.menu_File.addAction(self.action_Exit)
        self.menu_Observations.addAction(self.action_Edit_Observations)
        self.menu_Observations.addSeparator()
        self.menu_Observations.addAction(self.action_Corrections)
        self.menuStations.addAction(self.actionShow_Stations)
        self.menubar.addAction(self.menu_File.menuAction())
        self.menubar.addAction(self.menu_Observations.menuAction())
        self.menubar.addAction(self.menuEstimation_settings.menuAction())
        self.menubar.addAction(self.menuStations.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget_Main.setCurrentIndex(1)
        self.tab_Widget_Stations.setCurrentIndex(1)
        self.tabWidget_observations.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "GravTools"))
        self.groupBox_filter_options.setToolTip(_translate("MainWindow", "Display options for stations"))
        self.groupBox_filter_options.setTitle(_translate("MainWindow", "Filter Options"))
        self.checkBox_filter_observed_stat_only.setToolTip(_translate("MainWindow", "Display only stations that were observed in this campaign. "))
        self.checkBox_filter_observed_stat_only.setText(_translate("MainWindow", "Observed only"))
        self.label.setText(_translate("MainWindow", "Station name filter (regex)"))
        self.groupBox_edit_options.setTitle(_translate("MainWindow", "Edit Options"))
        self.tab_Widget_Stations.setTabText(self.tab_Widget_Stations.indexOf(self.tab_Stations_Table), _translate("MainWindow", "Table"))
        self.tab_Widget_Stations.setTabText(self.tab_Widget_Stations.indexOf(self.tab_Stations_Plot), _translate("MainWindow", "Plot"))
        self.tabWidget_Main.setTabText(self.tabWidget_Main.indexOf(self.tab_main_stations), _translate("MainWindow", "Stations"))
        self.treeWidget_observations.setSortingEnabled(False)
        self.treeWidget_observations.headerItem().setText(0, _translate("MainWindow", "Survey"))
        self.treeWidget_observations.headerItem().setText(1, _translate("MainWindow", "Station"))
        self.treeWidget_observations.headerItem().setText(2, _translate("MainWindow", "#Obs"))
        self.groupBox_obs_view_options.setTitle(_translate("MainWindow", "View Options"))
        self.checkBox_obs_show_all_columns.setText(_translate("MainWindow", "Show all columns"))
        self.checkBox_obs_plot_reduced_observations.setText(_translate("MainWindow", "Plot reduced observations"))
        self.tabWidget_observations.setTabText(self.tabWidget_observations.indexOf(self.tab_observations_table), _translate("MainWindow", "Table"))
        self.tabWidget_observations.setTabText(self.tabWidget_observations.indexOf(self.tab_observations_plots), _translate("MainWindow", "Plots"))
        self.tabWidget_Main.setTabText(self.tabWidget_Main.indexOf(self.tab_main_observations), _translate("MainWindow", "Observations"))
        self.tabWidget_Main.setTabText(self.tabWidget_Main.indexOf(self.tab_main_setups), _translate("MainWindow", "Setups"))
        self.tabWidget_Main.setTabToolTip(self.tabWidget_Main.indexOf(self.tab_main_setups), _translate("MainWindow", "Setup data: Corrected and reduces observation data accumulated for each instrument setup."))
        self.menu_File.setTitle(_translate("MainWindow", "File"))
        self.menuAdd_Survey.setTitle(_translate("MainWindow", "Add Survey"))
        self.menu_Observations.setTitle(_translate("MainWindow", "Observations"))
        self.menuEstimation_settings.setTitle(_translate("MainWindow", "Estimation settings"))
        self.menuStations.setTitle(_translate("MainWindow", "Stations"))
        self.action_New_Campaign.setText(_translate("MainWindow", "&New Campaign"))
        self.action_Save_Campaign.setText(_translate("MainWindow", "&Save Campaign"))
        self.action_Load_Campaign.setText(_translate("MainWindow", "&Load Campaign"))
        self.action_Add_Stations.setText(_translate("MainWindow", "Add Stations"))
        self.action_Exit.setText(_translate("MainWindow", "E&xit"))
        self.action_Exit.setToolTip(_translate("MainWindow", "Exit GravTools"))
        self.action_Edit_Observations.setText(_translate("MainWindow", "Edit Observations"))
        self.actionResiduals.setText(_translate("MainWindow", "Residuals"))
        self.actionEstimates.setText(_translate("MainWindow", "Estimates"))
        self.action_Corrections.setText(_translate("MainWindow", "Corrections"))
        self.actionShow_Stations.setText(_translate("MainWindow", "Show Stations"))
        self.action_from_CG5_observation_file.setText(_translate("MainWindow", "from CG5 observation file"))
        self.action_from_BEV_observation_file.setText(_translate("MainWindow", "from BEV observation file"))
from pyqtgraph import GraphicsLayoutWidget
