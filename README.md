# gravtools

Collection of newly developed software for gravity surveys. Programmed with Python 3, utilizing the scipy environment (numpy, pandas, matplotlib, scikit-learn,...).

Required packages can be found in requirements.txt.


# Installation and setup

**1.  Create virtual environment (venv)**

`python3 -m venv env`


**2. Activate virtual environment**

`source env/bin/activate`


**3. Clone git repository to local machine**

`git clone git@gitlab.com:Heller182/grav.git`

`cd grav`


**4. Install required python packages using pip**

`pip install -r requirements.txt`




# Components

## drift
Replacement for the currently used fortran (77) program "drift2011". In the first place the goal ist to deliver the same results as the old fortan program.

drift fits a polynomial of degree 1 to n (commonly 1 to 3) to gravity observations (tide corrected readings of the relative gravimeter CG5 (and others)) using multiple linear regression (L2 norm optimization).



# Changelog


# To Do


# Guidlines

The API reference is created with sphinx (https://www.sphinx-doc.org/).