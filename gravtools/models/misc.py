"""Misc methods for gravtools models."""
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

