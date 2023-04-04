"""Misc methods for GravTools models.

Copyright (C) 2021  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

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

from typing import Tuple
import numpy as np
import sys


def unique_ordered_list(seq: list) -> list:
    """Returns a unique list whilst keeping the order.

    Parameters
    ----------
    seq : list
        Input list.

    Returns
    -------
    list : Unique input list with same order of items.
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def is_singular_matrix(matrix) -> Tuple[bool, float]:
    """Check if the input matrix is singular and return the condition number.

    Parameters
    ----------
    matrix : np.array()
        Matrix (n,n) to be checked for singularity by checking the rank.

    Returns
    -------
    Tuple[bool, float] : `True`, if the matrix is singular and the matrix condition number.

    """
    cond = np.linalg.cond(matrix)
    if cond < 1 / sys.float_info.epsilon:  # not singular
        return False, cond
    else:  # is singular
        return True, cond


def print_matrix(matrix):
    """Print a matrix nicely.
    Parameters
    ----------
    matrix : np.array
        Input matrix.
    """
    for item in matrix:
        temp = (' '.join(f'{s:10.2f}' for s in item))
        print(temp)


def format_seconds_to_hhmmss(seconds):
    """Convert timespan in seconds to hh:mm:ss format.

    Parameters
    ----------
    seconds : float
        Timespan in seconds.

    Returns
    -------
    String: hh:mm:ss

    """
    hours = seconds // (60*60)
    seconds %= (60*60)
    minutes = seconds // 60
    seconds %= 60
    return f'{hours:02d}:{minutes:02d}:{seconds:02d}'
