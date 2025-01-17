# Makefile for project GravTools

# Package test (editable install in current virtual environment)
test_pack:
	pip install -e .

# Uninstall test package
test_pack_uninstall:
	pip uninstall gravtools

# Build package with setuptools (new version):
build:
	python -m build

# Upload package to pypi.org
pypi_push:
	twine upload --repository grav-toolbox --verbose dist/*

# Convert *.ui files from Qt Designer to Python files:
py_gui:
	pyuic5 -o gravtools/gui/MainWindow.py gravtools/gui/MainWindow.ui
	pyuic5 -o gravtools/gui/dialog_new_campaign.py gravtools/gui/dialog_new_campaign.ui
	pyuic5 -o gravtools/gui/dialog_load_stations.py gravtools/gui/dialog_load_stations.ui
	pyuic5 -o gravtools/gui/dialog_corrections.py gravtools/gui/dialog_corrections.ui
	pyuic5 -o gravtools/gui/dialog_autoselection_settings.py gravtools/gui/dialog_autoselection_settings.ui
	pyuic5 -o gravtools/gui/dialog_estimation_settings.py gravtools/gui/dialog_estimation_settings.ui
	pyuic5 -o gravtools/gui/dialog_export_results.py gravtools/gui/dialog_export_results.ui
	pyuic5 -o gravtools/gui/dialog_options.py gravtools/gui/dialog_options.ui
	pyuic5 -o gravtools/gui/dialog_setup_data.py gravtools/gui/dialog_setup_data.ui
	pyuic5 -o gravtools/gui/dialog_about.py gravtools/gui/dialog_about.ui
	pyuic5 -o gravtools/gui/dialog_gis_export_settings.py gravtools/gui/dialog_gis_export_settings.ui
	pyuic5 -o gravtools/gui/dialog_correction_time_series.py gravtools/gui/dialog_correction_time_series.ui
	pyuic5 -o gravtools/gui/dialog_load_tsf_file.py gravtools/gui/dialog_load_tsf_file.ui
	pyuic5 -o gravtools/gui/dialog_load_cg6_obs_files.py gravtools/gui/dialog_load_cg6_obs_files.ui
	pyuic5 -o gravtools/gui/dialog_calc_drift.py gravtools/gui/dialog_calc_drift.ui
