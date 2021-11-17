"""Misc methods for gravtools models."""


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
