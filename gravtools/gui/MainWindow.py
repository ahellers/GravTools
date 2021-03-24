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
        MainWindow.resize(1337, 729)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
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
        self.splitter_obs.setToolTipDuration(1)
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
        self.groupBox_obs_data_manipulation = QtWidgets.QGroupBox(self.layoutWidget)
        self.groupBox_obs_data_manipulation.setEnabled(False)
        self.groupBox_obs_data_manipulation.setObjectName("groupBox_obs_data_manipulation")
        self.verticalLayout_10 = QtWidgets.QVBoxLayout(self.groupBox_obs_data_manipulation)
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.pushButton_obs_apply_autoselect_current_data = QtWidgets.QPushButton(self.groupBox_obs_data_manipulation)
        self.pushButton_obs_apply_autoselect_current_data.setObjectName("pushButton_obs_apply_autoselect_current_data")
        self.verticalLayout_10.addWidget(self.pushButton_obs_apply_autoselect_current_data)
        self.pushButton_obs_comp_setup_data = QtWidgets.QPushButton(self.groupBox_obs_data_manipulation)
        self.pushButton_obs_comp_setup_data.setObjectName("pushButton_obs_comp_setup_data")
        self.verticalLayout_10.addWidget(self.pushButton_obs_comp_setup_data)
        self.pushButton_obs_run_estimation = QtWidgets.QPushButton(self.groupBox_obs_data_manipulation)
        self.pushButton_obs_run_estimation.setObjectName("pushButton_obs_run_estimation")
        self.verticalLayout_10.addWidget(self.pushButton_obs_run_estimation)
        self.verticalLayout_5.addWidget(self.groupBox_obs_data_manipulation)
        self.groupBox_obs_view_options = QtWidgets.QGroupBox(self.layoutWidget)
        self.groupBox_obs_view_options.setEnabled(False)
        self.groupBox_obs_view_options.setObjectName("groupBox_obs_view_options")
        self.verticalLayout_9 = QtWidgets.QVBoxLayout(self.groupBox_obs_view_options)
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.checkBox_obs_show_all_columns = QtWidgets.QCheckBox(self.groupBox_obs_view_options)
        self.checkBox_obs_show_all_columns.setChecked(True)
        self.checkBox_obs_show_all_columns.setObjectName("checkBox_obs_show_all_columns")
        self.verticalLayout_9.addWidget(self.checkBox_obs_show_all_columns)
        self.checkBox_obs_plot_reduced_observations = QtWidgets.QCheckBox(self.groupBox_obs_view_options)
        self.checkBox_obs_plot_reduced_observations.setChecked(True)
        self.checkBox_obs_plot_reduced_observations.setObjectName("checkBox_obs_plot_reduced_observations")
        self.verticalLayout_9.addWidget(self.checkBox_obs_plot_reduced_observations)
        self.checkBox_obs_plot_setup_data = QtWidgets.QCheckBox(self.groupBox_obs_view_options)
        self.checkBox_obs_plot_setup_data.setChecked(True)
        self.checkBox_obs_plot_setup_data.setObjectName("checkBox_obs_plot_setup_data")
        self.verticalLayout_9.addWidget(self.checkBox_obs_plot_setup_data)
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
        self.tab_observations_setups = QtWidgets.QWidget()
        self.tab_observations_setups.setObjectName("tab_observations_setups")
        self.verticalLayout_11 = QtWidgets.QVBoxLayout(self.tab_observations_setups)
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.tableView_observations_setups = QtWidgets.QTableView(self.tab_observations_setups)
        self.tableView_observations_setups.setObjectName("tableView_observations_setups")
        self.verticalLayout_11.addWidget(self.tableView_observations_setups)
        self.tabWidget_observations.addTab(self.tab_observations_setups, "")
        self.verticalLayout_7.addWidget(self.splitter_obs)
        self.tabWidget_Main.addTab(self.tab_main_observations, "")
        self.tab_main_results = QtWidgets.QWidget()
        self.tab_main_results.setToolTipDuration(2)
        self.tab_main_results.setObjectName("tab_main_results")
        self.verticalLayout_14 = QtWidgets.QVBoxLayout(self.tab_main_results)
        self.verticalLayout_14.setObjectName("verticalLayout_14")
        self.splitter_results = QtWidgets.QSplitter(self.tab_main_results)
        self.splitter_results.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_results.setObjectName("splitter_results")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.splitter_results)
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout_12 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_12.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_12.setObjectName("verticalLayout_12")
        self.groupBox_results_lsm_runs = QtWidgets.QGroupBox(self.verticalLayoutWidget)
        self.groupBox_results_lsm_runs.setObjectName("groupBox_results_lsm_runs")
        self.verticalLayout_16 = QtWidgets.QVBoxLayout(self.groupBox_results_lsm_runs)
        self.verticalLayout_16.setObjectName("verticalLayout_16")
        self.comboBox_results_lsm_run_selection = QtWidgets.QComboBox(self.groupBox_results_lsm_runs)
        self.comboBox_results_lsm_run_selection.setObjectName("comboBox_results_lsm_run_selection")
        self.verticalLayout_16.addWidget(self.comboBox_results_lsm_run_selection)
        self.pushButton_results_delete_lsm_run = QtWidgets.QPushButton(self.groupBox_results_lsm_runs)
        self.pushButton_results_delete_lsm_run.setObjectName("pushButton_results_delete_lsm_run")
        self.verticalLayout_16.addWidget(self.pushButton_results_delete_lsm_run)
        self.verticalLayout_12.addWidget(self.groupBox_results_lsm_runs)
        self.groupBox_results_data_selection = QtWidgets.QGroupBox(self.verticalLayoutWidget)
        self.groupBox_results_data_selection.setObjectName("groupBox_results_data_selection")
        self.verticalLayout_17 = QtWidgets.QVBoxLayout(self.groupBox_results_data_selection)
        self.verticalLayout_17.setObjectName("verticalLayout_17")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label_2 = QtWidgets.QLabel(self.groupBox_results_data_selection)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.comboBox_results_selection_station = QtWidgets.QComboBox(self.groupBox_results_data_selection)
        self.comboBox_results_selection_station.setObjectName("comboBox_results_selection_station")
        self.comboBox_results_selection_station.addItem("")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.comboBox_results_selection_station)
        self.label_5 = QtWidgets.QLabel(self.groupBox_results_data_selection)
        self.label_5.setObjectName("label_5")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_5)
        self.comboBox_results_selection_survey = QtWidgets.QComboBox(self.groupBox_results_data_selection)
        self.comboBox_results_selection_survey.setObjectName("comboBox_results_selection_survey")
        self.comboBox_results_selection_survey.addItem("")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.comboBox_results_selection_survey)
        self.verticalLayout_17.addLayout(self.formLayout)
        self.verticalLayout_12.addWidget(self.groupBox_results_data_selection)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_12.addItem(spacerItem1)
        self.tabWidget_results = QtWidgets.QTabWidget(self.splitter_results)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget_results.sizePolicy().hasHeightForWidth())
        self.tabWidget_results.setSizePolicy(sizePolicy)
        self.tabWidget_results.setObjectName("tabWidget_results")
        self.tab_results_info = QtWidgets.QWidget()
        self.tab_results_info.setObjectName("tab_results_info")
        self.verticalLayout_23 = QtWidgets.QVBoxLayout(self.tab_results_info)
        self.verticalLayout_23.setObjectName("verticalLayout_23")
        self.verticalLayout_22 = QtWidgets.QVBoxLayout()
        self.verticalLayout_22.setObjectName("verticalLayout_22")
        self.groupBox_results_info = QtWidgets.QGroupBox(self.tab_results_info)
        self.groupBox_results_info.setObjectName("groupBox_results_info")
        self.verticalLayout_21 = QtWidgets.QVBoxLayout(self.groupBox_results_info)
        self.verticalLayout_21.setObjectName("verticalLayout_21")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.label_3 = QtWidgets.QLabel(self.groupBox_results_info)
        self.label_3.setObjectName("label_3")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.label_results_adjustment_method = QtWidgets.QLabel(self.groupBox_results_info)
        self.label_results_adjustment_method.setMinimumSize(QtCore.QSize(100, 0))
        self.label_results_adjustment_method.setText("")
        self.label_results_adjustment_method.setObjectName("label_results_adjustment_method")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.label_results_adjustment_method)
        self.label_4 = QtWidgets.QLabel(self.groupBox_results_info)
        self.label_4.setObjectName("label_4")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_4)
        self.label_results_time_and_date = QtWidgets.QLabel(self.groupBox_results_info)
        self.label_results_time_and_date.setMinimumSize(QtCore.QSize(100, 0))
        self.label_results_time_and_date.setText("")
        self.label_results_time_and_date.setObjectName("label_results_time_and_date")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.label_results_time_and_date)
        self.label_10 = QtWidgets.QLabel(self.groupBox_results_info)
        self.label_10.setObjectName("label_10")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_10)
        self.label_results_comment = QtWidgets.QLabel(self.groupBox_results_info)
        self.label_results_comment.setMinimumSize(QtCore.QSize(100, 0))
        self.label_results_comment.setText("")
        self.label_results_comment.setObjectName("label_results_comment")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.label_results_comment)
        self.verticalLayout_21.addLayout(self.formLayout_2)
        self.verticalLayout_22.addWidget(self.groupBox_results_info)
        self.groupBox_results_statistics = QtWidgets.QGroupBox(self.tab_results_info)
        self.groupBox_results_statistics.setObjectName("groupBox_results_statistics")
        self.verticalLayout_24 = QtWidgets.QVBoxLayout(self.groupBox_results_statistics)
        self.verticalLayout_24.setObjectName("verticalLayout_24")
        self.formLayout_3 = QtWidgets.QFormLayout()
        self.formLayout_3.setObjectName("formLayout_3")
        self.label_9 = QtWidgets.QLabel(self.groupBox_results_statistics)
        self.label_9.setObjectName("label_9")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_9)
        self.label_results_sig0 = QtWidgets.QLabel(self.groupBox_results_statistics)
        self.label_results_sig0.setMinimumSize(QtCore.QSize(100, 0))
        self.label_results_sig0.setText("")
        self.label_results_sig0.setObjectName("label_results_sig0")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.label_results_sig0)
        self.label_6 = QtWidgets.QLabel(self.groupBox_results_statistics)
        self.label_6.setObjectName("label_6")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_6)
        self.label_results_goodness_of_fit_test_status = QtWidgets.QLabel(self.groupBox_results_statistics)
        self.label_results_goodness_of_fit_test_status.setText("")
        self.label_results_goodness_of_fit_test_status.setObjectName("label_results_goodness_of_fit_test_status")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.label_results_goodness_of_fit_test_status)
        self.label_7 = QtWidgets.QLabel(self.groupBox_results_statistics)
        self.label_7.setObjectName("label_7")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_7)
        self.label_results_numbr_of_outliers = QtWidgets.QLabel(self.groupBox_results_statistics)
        self.label_results_numbr_of_outliers.setText("")
        self.label_results_numbr_of_outliers.setObjectName("label_results_numbr_of_outliers")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.label_results_numbr_of_outliers)
        self.verticalLayout_24.addLayout(self.formLayout_3)
        self.verticalLayout_22.addWidget(self.groupBox_results_statistics)
        self.groupBox_results_lsm_run_log = QtWidgets.QGroupBox(self.tab_results_info)
        self.groupBox_results_lsm_run_log.setObjectName("groupBox_results_lsm_run_log")
        self.verticalLayout_20 = QtWidgets.QVBoxLayout(self.groupBox_results_lsm_run_log)
        self.verticalLayout_20.setObjectName("verticalLayout_20")
        self.plainTextEdit_results_log = QtWidgets.QPlainTextEdit(self.groupBox_results_lsm_run_log)
        self.plainTextEdit_results_log.setReadOnly(True)
        self.plainTextEdit_results_log.setObjectName("plainTextEdit_results_log")
        self.verticalLayout_20.addWidget(self.plainTextEdit_results_log)
        self.verticalLayout_22.addWidget(self.groupBox_results_lsm_run_log)
        self.verticalLayout_23.addLayout(self.verticalLayout_22)
        self.tabWidget_results.addTab(self.tab_results_info, "")
        self.tab_results_obs_table = QtWidgets.QWidget()
        self.tab_results_obs_table.setObjectName("tab_results_obs_table")
        self.verticalLayout_15 = QtWidgets.QVBoxLayout(self.tab_results_obs_table)
        self.verticalLayout_15.setObjectName("verticalLayout_15")
        self.tableView_results_observations = QtWidgets.QTableView(self.tab_results_obs_table)
        self.tableView_results_observations.setObjectName("tableView_results_observations")
        self.verticalLayout_15.addWidget(self.tableView_results_observations)
        self.tabWidget_results.addTab(self.tab_results_obs_table, "")
        self.tab_results_stat_table = QtWidgets.QWidget()
        self.tab_results_stat_table.setObjectName("tab_results_stat_table")
        self.verticalLayout_19 = QtWidgets.QVBoxLayout(self.tab_results_stat_table)
        self.verticalLayout_19.setObjectName("verticalLayout_19")
        self.tableView_results_stations = QtWidgets.QTableView(self.tab_results_stat_table)
        self.tableView_results_stations.setObjectName("tableView_results_stations")
        self.verticalLayout_19.addWidget(self.tableView_results_stations)
        self.tabWidget_results.addTab(self.tab_results_stat_table, "")
        self.tab_results_drift = QtWidgets.QWidget()
        self.tab_results_drift.setObjectName("tab_results_drift")
        self.verticalLayout_18 = QtWidgets.QVBoxLayout(self.tab_results_drift)
        self.verticalLayout_18.setObjectName("verticalLayout_18")
        self.tableView_results_drift = QtWidgets.QTableView(self.tab_results_drift)
        self.tableView_results_drift.setObjectName("tableView_results_drift")
        self.verticalLayout_18.addWidget(self.tableView_results_drift)
        self.tabWidget_results.addTab(self.tab_results_drift, "")
        self.tab_results_obs_plots = QtWidgets.QWidget()
        self.tab_results_obs_plots.setObjectName("tab_results_obs_plots")
        self.verticalLayout_13 = QtWidgets.QVBoxLayout(self.tab_results_obs_plots)
        self.verticalLayout_13.setObjectName("verticalLayout_13")
        self.graphicsLayoutWidget_results_observations_plots = GraphicsLayoutWidget(self.tab_results_obs_plots)
        self.graphicsLayoutWidget_results_observations_plots.setObjectName("graphicsLayoutWidget_results_observations_plots")
        self.verticalLayout_13.addWidget(self.graphicsLayoutWidget_results_observations_plots)
        self.groupBox_plot_settings = QtWidgets.QGroupBox(self.tab_results_obs_plots)
        self.groupBox_plot_settings.setObjectName("groupBox_plot_settings")
        self.verticalLayout_25 = QtWidgets.QVBoxLayout(self.groupBox_plot_settings)
        self.verticalLayout_25.setObjectName("verticalLayout_25")
        self.formLayout_4 = QtWidgets.QFormLayout()
        self.formLayout_4.setObjectName("formLayout_4")
        self.label_8 = QtWidgets.QLabel(self.groupBox_plot_settings)
        self.label_8.setObjectName("label_8")
        self.formLayout_4.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_8)
        self.comboBox_results_obs_plot_select_data_column = QtWidgets.QComboBox(self.groupBox_plot_settings)
        self.comboBox_results_obs_plot_select_data_column.setObjectName("comboBox_results_obs_plot_select_data_column")
        self.formLayout_4.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.comboBox_results_obs_plot_select_data_column)
        self.verticalLayout_25.addLayout(self.formLayout_4)
        self.verticalLayout_13.addWidget(self.groupBox_plot_settings)
        self.tabWidget_results.addTab(self.tab_results_obs_plots, "")
        self.verticalLayout_14.addWidget(self.splitter_results)
        self.tabWidget_Main.addTab(self.tab_main_results, "")
        self.horizontalLayout.addWidget(self.tabWidget_Main)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1337, 22))
        self.menubar.setObjectName("menubar")
        self.menu_File = QtWidgets.QMenu(self.menubar)
        self.menu_File.setObjectName("menu_File")
        self.menuAdd_Survey = QtWidgets.QMenu(self.menu_File)
        self.menuAdd_Survey.setEnabled(False)
        self.menuAdd_Survey.setObjectName("menuAdd_Survey")
        self.menu_Observations = QtWidgets.QMenu(self.menubar)
        self.menu_Observations.setEnabled(False)
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
        self.action_Edit_Observations.setEnabled(True)
        self.action_Edit_Observations.setObjectName("action_Edit_Observations")
        self.actionResiduals = QtWidgets.QAction(MainWindow)
        self.actionResiduals.setObjectName("actionResiduals")
        self.actionEstimates = QtWidgets.QAction(MainWindow)
        self.actionEstimates.setObjectName("actionEstimates")
        self.action_Corrections = QtWidgets.QAction(MainWindow)
        self.action_Corrections.setEnabled(True)
        self.action_Corrections.setObjectName("action_Corrections")
        self.actionShow_Stations = QtWidgets.QAction(MainWindow)
        self.actionShow_Stations.setObjectName("actionShow_Stations")
        self.action_from_CG5_observation_file = QtWidgets.QAction(MainWindow)
        self.action_from_CG5_observation_file.setObjectName("action_from_CG5_observation_file")
        self.action_from_BEV_observation_file = QtWidgets.QAction(MainWindow)
        self.action_from_BEV_observation_file.setObjectName("action_from_BEV_observation_file")
        self.action_Autoselection_settings = QtWidgets.QAction(MainWindow)
        self.action_Autoselection_settings.setObjectName("action_Autoselection_settings")
        self.action_Estimation_settings = QtWidgets.QAction(MainWindow)
        self.action_Estimation_settings.setObjectName("action_Estimation_settings")
        self.actionEstimate_long_term_drift = QtWidgets.QAction(MainWindow)
        self.actionEstimate_long_term_drift.setEnabled(True)
        self.actionEstimate_long_term_drift.setObjectName("actionEstimate_long_term_drift")
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
        self.menu_Observations.addAction(self.action_Autoselection_settings)
        self.menu_Observations.addAction(self.actionEstimate_long_term_drift)
        self.menuEstimation_settings.addAction(self.action_Estimation_settings)
        self.menuStations.addAction(self.actionShow_Stations)
        self.menubar.addAction(self.menu_File.menuAction())
        self.menubar.addAction(self.menu_Observations.menuAction())
        self.menubar.addAction(self.menuEstimation_settings.menuAction())
        self.menubar.addAction(self.menuStations.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget_Main.setCurrentIndex(2)
        self.tab_Widget_Stations.setCurrentIndex(1)
        self.tabWidget_observations.setCurrentIndex(1)
        self.tabWidget_results.setCurrentIndex(4)
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
        self.groupBox_obs_data_manipulation.setTitle(_translate("MainWindow", "Data manipulation"))
        self.pushButton_obs_apply_autoselect_current_data.setToolTip(_translate("MainWindow", "Activate/deactivate observations of the current selected survey or setup according to the auto-selection options."))
        self.pushButton_obs_apply_autoselect_current_data.setText(_translate("MainWindow", "Apply Autoselection"))
        self.pushButton_obs_comp_setup_data.setToolTip(_translate("MainWindow", "Compute variance-weighted mean of active observations for all setups in campaign."))
        self.pushButton_obs_comp_setup_data.setText(_translate("MainWindow", "Compute setup data for campaign"))
        self.pushButton_obs_run_estimation.setToolTip(_translate("MainWindow", "Run parameter estimation according to the estimation settings."))
        self.pushButton_obs_run_estimation.setText(_translate("MainWindow", "Run estimation"))
        self.groupBox_obs_view_options.setTitle(_translate("MainWindow", "View Options"))
        self.checkBox_obs_show_all_columns.setText(_translate("MainWindow", "Show all columns"))
        self.checkBox_obs_plot_reduced_observations.setToolTip(_translate("MainWindow", "Plot reduced observatio data, if availble."))
        self.checkBox_obs_plot_reduced_observations.setText(_translate("MainWindow", "Plot reduced observations"))
        self.checkBox_obs_plot_setup_data.setToolTip(_translate("MainWindow", "Plot variance-weighted mean of observation data for all setups."))
        self.checkBox_obs_plot_setup_data.setText(_translate("MainWindow", "Plot setup data"))
        self.tabWidget_observations.setTabText(self.tabWidget_observations.indexOf(self.tab_observations_table), _translate("MainWindow", "Table"))
        self.tabWidget_observations.setTabText(self.tabWidget_observations.indexOf(self.tab_observations_plots), _translate("MainWindow", "Plots"))
        self.tabWidget_observations.setTabText(self.tabWidget_observations.indexOf(self.tab_observations_setups), _translate("MainWindow", "Setups"))
        self.tabWidget_Main.setTabText(self.tabWidget_Main.indexOf(self.tab_main_observations), _translate("MainWindow", "Observations"))
        self.groupBox_results_lsm_runs.setTitle(_translate("MainWindow", "LSM runs"))
        self.pushButton_results_delete_lsm_run.setText(_translate("MainWindow", "Delete selected lsm run"))
        self.groupBox_results_data_selection.setTitle(_translate("MainWindow", "Data selection"))
        self.label_2.setText(_translate("MainWindow", "Station"))
        self.comboBox_results_selection_station.setItemText(0, _translate("MainWindow", "All stations"))
        self.label_5.setText(_translate("MainWindow", "Survey"))
        self.comboBox_results_selection_survey.setItemText(0, _translate("MainWindow", "All surveys"))
        self.groupBox_results_info.setTitle(_translate("MainWindow", "Info"))
        self.label_3.setText(_translate("MainWindow", "Adjustment method:"))
        self.label_4.setText(_translate("MainWindow", "Time and Date:"))
        self.label_10.setText(_translate("MainWindow", "Comment: "))
        self.groupBox_results_statistics.setTitle(_translate("MainWindow", "Statistics and tests"))
        self.label_9.setText(_translate("MainWindow", "sig0 = "))
        self.label_6.setText(_translate("MainWindow", "Goodness-of-fit test:"))
        self.label_7.setText(_translate("MainWindow", "Number of outliers:"))
        self.groupBox_results_lsm_run_log.setTitle(_translate("MainWindow", "LSM run log"))
        self.tabWidget_results.setTabText(self.tabWidget_results.indexOf(self.tab_results_info), _translate("MainWindow", "Info"))
        self.tabWidget_results.setTabText(self.tabWidget_results.indexOf(self.tab_results_obs_table), _translate("MainWindow", "Observations table"))
        self.tabWidget_results.setTabText(self.tabWidget_results.indexOf(self.tab_results_stat_table), _translate("MainWindow", " Stations table"))
        self.tabWidget_results.setTabText(self.tabWidget_results.indexOf(self.tab_results_drift), _translate("MainWindow", "Drift table"))
        self.groupBox_plot_settings.setTitle(_translate("MainWindow", "Plot settings"))
        self.label_8.setText(_translate("MainWindow", "Select data column:"))
        self.tabWidget_results.setTabText(self.tabWidget_results.indexOf(self.tab_results_obs_plots), _translate("MainWindow", "Observations Plots"))
        self.tabWidget_Main.setTabText(self.tabWidget_Main.indexOf(self.tab_main_results), _translate("MainWindow", "Results"))
        self.tabWidget_Main.setTabToolTip(self.tabWidget_Main.indexOf(self.tab_main_results), _translate("MainWindow", "Setup data: Corrected and reduces observation data accumulated for each instrument setup."))
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
        self.action_Autoselection_settings.setText(_translate("MainWindow", "Autoselection settings"))
        self.action_Estimation_settings.setText(_translate("MainWindow", "Settings"))
        self.actionEstimate_long_term_drift.setText(_translate("MainWindow", "Estimate long-term drift"))
from pyqtgraph import GraphicsLayoutWidget
