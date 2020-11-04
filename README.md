# gravtools
Collection of newly developed software for gravity surveys. Programmed with Python 3, utilizing the scipy environment (numpy, pandas, matplotlib, scikit-learn,...).

# Installation and setup
* **1.  Create virtual environment (venv)**
  * `python3 -m venv env`
* **2. Activate virtual environment**
  * `source env/bin/activate`
* **3. Clone git repository to local machine**
  * `git clone git@gitlab.com:Heller182/grav.git`
  * `cd grav`
* **4. Install required python packages using pip**
  * `pip install -r requirements.txt`


# Test installation with setuptools 
With command line interface.

* **1. Configure setup.py**
  * Define entry points (*console_scripts*)
* **2. Activate virtual environment**
  * e.g. `source env/bin/activate`
* **3. Run setup.py**
  * `python3 setup.py develop`

## Using make
`python3 setup.py develop`


# Create HTML documentation with sphinx:
Run make in the gravtools/doc directory: 
* `>>>make html_doc`

# Components

## drift
Replacement for the currently used fortran (77) program "drift2011". In the first place the goal ist to deliver the same results as the old Fortan program.
Drift fits a polynomial of degree 1 to n (commonly 1 to 3) to gravity observations (tide corrected readings of the relative gravimeter CG5 (and others)) using multiple linear regression (L2 norm optimization).


# Guidelines and conventions

## Code style:
* Respect the PEP conventions on python coding!
* The maximum line length is 120 characters
* Use **type hints**: https://www.python.org/dev/peps/pep-0484/
* Use docstrings according to the numpy standard: https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
  * They are useful to generate the documentation automatically
* Comment code, if necessary!
* Use English language for the code, docstrings and comments
  * German is allowed for user interfaces (GUI, command line), although English is preferred

## Documentation and docstring style
* The API reference is created with sphinx (https://www.sphinx-doc.org/).
* Docstrings have to follow the numpy standard, see: https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
  * Examples: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
* Package documentation via docstring in __ini__.py files
* Module documentation via docstring at first lines of py-file
* Documentation of classes, class methods and functions via docstrings
  
## Command line interface and executable scripts
* The command line interface is realized via entry points (console_scripts) in setuptools (python packaging tool)
  * Input arguments are handled with argparse
  * The code is located in the command_line module (gravtools/command_line.py)
* Executable scripts are located in gravtools/scripts

## Dependancies
* Required python packages are listed in gravtools/requirements.txt
  * created with `>>>pip freeze > requirements.txt`
  
## Version control with GIT
* Gitlab repository: https://gitlab.com/Heller182/grav
* Branching model:
  * **master** branch: Current release version
  * **develop** branch: Current working version. 
    * All team members merge their feature branches into develop (merge request via gitlab)
    * Make sure that the develop branch contains a fully functional version of the code!
  * **feature** branches: Branches of develop for the implementation of new features and other changes.
    * Code changes only in feature branches!
    * Naming convention: feature_<description of change/feature>, e.g. feature_new_tide_model
* Use gitignore files to prevent any data files (except example files), IDE control files, compiled python code, etc. from being stored in the GIT repository
    
## Packaging and distribution
* With setuptools 
    
    
# Changelog


# To Do

