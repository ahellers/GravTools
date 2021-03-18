# Makefile for project gravtools

# Test run
test:
	$(info Test-run for this makefile!)
	$(info Yeah!!)

# Project initaialization
init:
	pip install -r requirements.txt

# test packaging (install in current virtual environment)
test_pack:
	python3 setup.py develop

# schwaus test (with actual dataset)
schwaus_test:
	schwaus --obs-file data/BEV/n20200701_1 --out-dir /tmp
	cat /tmp/n20200701_1_prot.txt
	cat /tmp/n20200701_1.nsb
	pdfinfo /tmp/n20200701_1_drift.pdf

# Convert *.ui files from Qt Designer to Python files:
py_gui:
	pyuic5 -o gravtools/gui/MainWindow.py gravtools/gui/MainWindow.ui
	pyuic5 -o gravtools/gui/dialog_new_campaign.py gravtools/gui/dialog_new_campaign.ui
	pyuic5 -o gravtools/gui/dialog_load_stations.py gravtools/gui/dialog_load_stations.ui
	pyuic5 -o gravtools/gui/dialog_corrections.py gravtools/gui/dialog_corrections.ui
	pyuic5 -o gravtools/gui/dialog_autoselection_settings.py gravtools/gui/dialog_autoselection_settings.ui
	pyuic5 -o gravtools/gui/dialog_estimation_settings.py gravtools/gui/dialog_estimation_settings.ui