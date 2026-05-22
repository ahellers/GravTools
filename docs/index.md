# GravTools

GravTools (`grav-toolbox` on PyPI) is a Python toolbox for processing relative gravity surveys,
developed at Austria's Federal Office of Metrology and Surveying (BEV).

## Features

- Load gravimeter observation files (Scintrex CG-5, CG-6 in multiple formats)
- Manage survey stations and gravimeter metadata
- Apply tidal, atmospheric pressure, and height corrections
- Run least-squares adjustment to estimate gravity values at stations
- Estimate drift polynomials per survey
- Export results to CSV and GIS formats (shapefile via geopandas)

## Installation

```bash
pip install grav-toolbox
```

With optional GIS support:

```bash
pip install "grav-toolbox[gis]"
```

## Source Code

[https://github.com/ahellers/GravTools](https://github.com/ahellers/GravTools)


## Python Package

[https://pypi.org/project/grav-toolbox/](https://pypi.org/project/grav-toolbox/)


## Citing GravTools

If you use datasets computed with GravTools in a publication, please cite the following reference:

> Hellerschmied, A., Zehetmaier, B. (2023): *GravTools: A software toolbox for the analysis of relative gravity surveys*, XXVIII General Assembly of the International Union of Geodesy and Geophysics (IUGG) (Berlin 2023). [https://doi.org/10.57757/IUGG23-0864](https://doi.org/10.57757/IUGG23-0864)


## License and copyright

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