"""Models for the calculation of atmospheric pressure corrections.

Copyright (C) 2023  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

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
"""
from typing import Tuple, Any

def normal_air_pressure_iso(height_m: float) -> float:
    """Returns the normal air pressure according to ISO 2533:1975.

    Parameters
    ----------
    height_m : float
        Physical height [m].

    Notes
    -----
    Reference:
    Wziontek et al. (2021): Status of the International Gravity Reference System and Frame, Journal of
    Geodesy, 95:7, https://doi.org/10.1007/s00190-020-01438-9
    Torge (1989): Gravimetry. De Gruyter, ISBN 3-11-010702-3, p. 323, Ch. 9.2.3

    Returns
    -------
    float
        Normal air pressure [hPa]

    """
    return 1013.25 * (1 - 0.0065 * height_m / 288.15) ** 5.2559


def pressure_correction_iso(height_m: float, p_obs_hpa: float, admittance: float = 0.3) -> tuple[float | Any, float]:
    """Returns the atmospheric pressure correction determined by the station height and ISO normal air pressure.

    Parameters
    ----------
    height_m : float
        Physical height [m].
    p_obs_hpa : float
        Observed air pressure [hPa]
    admittance : float, optional (default = 0.3)
        Admittance factor for the determination of pressure corrections based on the difference between measured and
        normal air pressure.

    Notes
    -----
    Reference:
    Wziontek et al. (2021): Status of the International Gravity Reference System and Frame, Journal of
    Geodesy, 95:7, https://doi.org/10.1007/s00190-020-01438-9
    Torge (1989): Gravimetry. De Gruyter, ISBN 3-11-010702-3, p. 323, Ch. 9.2.3

    Returns
    -------
    dg_mugal : float
        Atmopheric pressure correction [µGal].
    pn_hpa : float
        Normal air pressure according to ISO 2533:1975 [hPa].
    """
    pn_hpa = normal_air_pressure_iso(height_m)
    dp_hpa = p_obs_hpa - pn_hpa
    dg_mugal = admittance * dp_hpa
    return dg_mugal, pn_hpa


def pressure_correction_iso_pandas_series(height_m, p_obs_hpa, admittance: float = 0.3):
    """Returns the atmospheric pressure correction determined by the station height and ISO normal air pressure.

        Parameters
        ----------
        height_m : pandas Series of float values
            Physical height [m].
        p_obs_hpa : pandas Series of float values
            Observed air pressure [hPa]
        admittance : float, optional (default = 0.3)
            Admittance factor for the determination of pressure corrections based on the difference between measured and
            normal air pressure.

        Notes
        -----
        Reference:
        Wziontek et al. (2021): Status of the International Gravity Reference System and Frame, Journal of
        Geodesy, 95:7, https://doi.org/10.1007/s00190-020-01438-9
        Torge (1989): Gravimetry. De Gruyter, ISBN 3-11-010702-3, p. 323, Ch. 9.2.3

        Returns
        -------
        dg_mugal : float
            Atmospheric pressure correction [µGal].
        pn_hpa : float
            Normal air pressure according to ISO 2533:1975 [hPa]. `NaN` is returned for each item with missing observed
            pressure value in `p_obs_hpa`.
        """

    pn_hpa = height_m.apply(normal_air_pressure_iso)
    dp_hpa = p_obs_hpa - pn_hpa
    dg_mugal = admittance * dp_hpa
    return dg_mugal, pn_hpa

