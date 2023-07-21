# GravTools
GravTools is an open source software toolbox for processing relative gravity surveys developed at the Austrian
Federal Office of Metrology and Surveying (BEV).  

 - The source code is hosted on github.com: https://github.com/ahellers/GravTools
 - The python package is hosted at: https://pypi.org/project/grav-toolbox/

# Installation

## Install python package with pip:
`pip install grav-toolbox`

### Optional dependency for GIS data export
GravTools allows to export station- und observation-results to shapefiles for import and analysis in GIS programs.
To enable these features the optional package "geopandas" needs to be installed by executing:
`pip install grav-toolbox[gis]`
 
## Installation issues on Centos 8:
* **1. Carry out installation steps describes above** 
* **2. Install QT5, if required**
  * `sudo yum install qt5-qtbase-devel.x86_64`
  
## Installation issues on Ubuntu 20.04:
After just installing PyQt5 with pip3 the following error occurred when trying to actually run a PyQt GUI: qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This issue was resolved by installing the QT dev tools (Designer, etc.): 
`sudo apt-get install qttools5-dev-tools`

## Test installation
In a virtual environment.

* Editable installation with pip: `pip install -e .`
  * Or: `make test_pack`
* Uninstall: `pip uninstall grav-toolbox`
  * Or: `test_pack_uninstall`

# Packaging
With setuptools and the build package.

## Create a source and binary distribution (sdist and wheel)
* Set up setup.cfg and pyproject.toml in the project root directory
* With make and the predefined makefile: `make build`
* Without make: `python3 -m build`

The created package (wheel and source distribution) is located in the "dist" directory.

## Push package to pypi.org
Define user credentials in /home/.pypirc

* `twine upload --verbose dist/*` or
* `make pypi_push`


# Create a stand-alone Windows executable
For creating a Windows executable (stand-alone application, without python installation) the 
package "auto-py-to-exe" can be used (see: https://pypi.org/project/auto-py-to-exe/). This is a 
simple graphical user interface (based on the chrome browser) that allows to define all settings 
for the package "PyInstaller" (see: https://www.pyinstaller.org/). 

Creating a Windows executable with auto-py-to-exe was tested on Windows 10. 
Follow these steps to create an executable on a Windows machine:
* **1. Install Python3**
  * With installer from https://www.python.org/downloads/windows/
  * Add python3.x to the Windows search path and install pip!
* **2. Install virtualenv and create a virtual environment** in CMD
  * a. `pip install virtualenv`
  * b. Change to project directory (`cd ..`)
  * c. Create virtual environment: `virtualenv env`
  * d. Activate it: `env\Scripts\activate.bat`
    * Deactivate with `deactivate`
* **3. Install gravtools**: `pip install gravtools[gis]`
* **4. Install auto-py-to-exe**: `pip install auto-py-to-exe`
* **5. Try to run the GravTools**
  * Type `gt` in the command line interface (virtual environment must be active)
* **6. Create exe with auto-py-to-exe**
  * Run the CMD Window as administrator!
  * a. Start auto-py-to-exe in CMD: `auto-py-to-exe`
  * b. Select the script location (select: gravtools/scripts/run_gui.py)
  * c. Select "One File" and "Console based" (in addition to teh GUI a console will appear)
  * d. Start conversion py pressing the big blue button on the GUI bottom
  * e. The exe file will be save at a new "output" directory. Move the file to: "Windows_executables"

# Create HTML documentation with sphinx:
Sphinx is used to create an API documentation based on docstrings. Run make in the gravtools/doc directory: 
* `>>>make html_doc`

# Guidelines and conventions

## Code style:
* Respect the PEP conventions on python coding!
  * PEP 8 -- Style Guide for Python Code: https://www.python.org/dev/peps/pep-0008/
* The maximum line length is 120 characters
* Use **type hints**: https://www.python.org/dev/peps/pep-0484/
* Use docstrings according to the numpy standard: https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
  * They are useful to generate the documentation automatically
  * Example: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
* Comment code, if necessary!
* Use English language for the code, docstrings and comments

## Documentation and docstring style
* The API reference is created with sphinx (https://www.sphinx-doc.org/).
* Docstrings have to follow the numpy standard, see: https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
  * Examples: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
  * Contrary to the numpy standard license and copyright information is provided at least in the package and module docstrings.
* Package documentation via docstring in __ini__.py files
* Module documentation via docstring at first lines of py-file
* Documentation of classes, class methods and functions via docstrings
  
## Command line interface and executable scripts
* The command line interface is realized via entry points (console_scripts) in setuptools (python packaging tool)
  * Input arguments are handled with argparse
* Executable scripts are located in gravtools/scripts
  

# License and copyright

Copyright (C) 2021  Andreas Hellerschmied (<andreas.hellerschmied@bev.gv.at>)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.


