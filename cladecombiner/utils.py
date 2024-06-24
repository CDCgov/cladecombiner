from collections.abc import Iterable
from typing import Any


def table(x: Iterable) -> dict[Any, int]:
    """
    Like R's base::table(), counts occurrences of elements in container

    Parameters
    ---------
    x : Iterable
        The container to be tabulated.

    Returns
    -------
    dict
        A dictionary of form object : count, where object is some element in x
        and count is an int specifying the number of times that object is seen
        in x.
    """
    unique = set(x)
    res = {}
    for u in unique:
        res[u] = 0
    for item in x:
        res[item] += 1
    return res


def table_equal(x1: Iterable, x2: Iterable) -> bool:
    """
    Checks if x1 and x2 have the same items the same number of times

    That is, is table(x1) equivalent to table(x2)?

    Parameters
    ---------
    x1 : Iterable
        One container to be checked.
    x2 : Iterable
        The other container to be checked.

    Returns
    -------
    bool
        True if both x1 and x2 contain the same objects the same number of
        times each.
    """
    t1 = table(x1)
    t2 = table(x2)
    if t1.keys() != t2.keys():
        return False
    for k1, v1 in t1.items():
        if v1 != t2[k1]:
            return False
    return True


def table_index_map(x: Iterable) -> dict[Any, list[int]]:
    """
    Like table, but maps to the indices with those elements

    Parameters
    ---------
    x : Iterable
        The container to be tabulated.

    Returns
    -------
    dict
        A dictionary of form object : [indices], where object is some element
        in x and [indices] is a list of indices in x where that object occurs.
    """
    unique = set(x)
    res: dict[Any, list[int]] = {}
    for u in unique:
        res[u] = []
    idx = 0
    for item in x:
        res[item] += [idx]
        idx += 1
    return res
