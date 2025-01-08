# GravTools
GravTools is an open source software toolbox for the analysis of relative gravity surveys. 
GravTools is developed at Austria's Federal Office of Metrology and Surveying (BEV).  

 - The source code is hosted on github.com: https://github.com/ahellers/GravTools
 - The python package is published on pypi.org: https://pypi.org/project/grav-toolbox/


# Installation

The python 3 package is available on the Python package index website (pypi.org) and can be easily installed using pip.

## Install python package with pip:
`pip install grav-toolbox`

## Optional dependency for GIS data export
GravTools allows users to export adjustment results to shapefiles for data visualization in external GIS programs (e.g. QGIS).
To enable these features the optional package "geopandas" needs to be installed by executing:
`pip install grav-toolbox[gis]`

# Release notes
More details are provided along with version tags on github.

## 0.3.3 (2025-01-08)
  - New dialog for drift determination based on the reduced observations as shown in the observations plot. One drift 
polynomial is fitted per station. This feature is useful for calculating the drift correction coefficient for Scintrex 
CG5 and CG6 meters based on stationary observations. 

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


