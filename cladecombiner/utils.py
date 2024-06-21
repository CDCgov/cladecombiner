from collections.abc import Iterable, Sequence
from os import path
from typing import Any, Optional

from .nomenclature import Nomenclature
from .taxon import Taxon
from .taxonomy_scheme import TaxonomyScheme


def read_taxa(
    fp: str,
    is_tip: bool | Sequence[bool],
    nomenclature: Optional[Nomenclature],
    taxonomy_scheme: Optional[TaxonomyScheme],
) -> Sequence[Taxon]:
    """
    Reads in taxa as a list of Taxon objects.

    Parameters
    ---------
    fp : str
        The file path to be read from
    is_tip : bool | Sequence[bool]
        Either one bool specifying whether all these are tip taxa or not,
        or one bool per taxon in the file specifying for each.
    nomenclature : Optional[Nomenclature]
        If specified, taxon names are checked for validity according to this
        nomenclature scheme, and an error is raised if an invalid taxon is
        found.
    taxonomy_scheme : Optional[TaxonomyScheme]
        If specified, taxon names are checked for validity according to this
        taxonomy scheme, and an error is raised if an invalid taxon is found.

    Returns
    -------
    Sequence[Taxon]
        Container of the taxa as Taxon objects.
    """

    ext = path.splitext(fp)[1]
    taxa = []
    if ext == ".txt":
        f = open(fp)
        lines = f.readlines()
        f.close()
        taxa = []
        if not isinstance(is_tip, Sequence):
            is_tip = [is_tip for _ in range(len(lines))]

        for line, i in zip(lines, range(len(lines))):
            if nomenclature:
                if not nomenclature.is_valid_name(line[:-1]):
                    raise RuntimeError(
                        "The name "
                        + line[:-1]
                        + " is not valid under the provided nomenclature ("
                        + str(nomenclature)
                        + ")"
                    )
            taxon = Taxon(line[:-1], is_tip[i])
            if taxonomy_scheme:
                if not taxonomy_scheme.is_valid_taxon(taxon):
                    raise RuntimeError(
                        "The name "
                        + str(taxon)
                        + " is not valid under the provided taxonomy scheme ("
                        + str(taxonomy_scheme)
                        + ")"
                    )
            taxa.append(taxon)
    return taxa


def printable_taxon_list(taxa: Sequence[Taxon], sep: str = "\n") -> str:
    """
    Prettier printing of lists of taxa.

    Parameters
    ---------
    taxa : sequence[Taxon]
        The Taxon objects to be printed.
    sep : str
        The separator for printing the list

    Returns
    -------
    str
        A string which may be fed to print().
    """
    print_str = ""
    for taxon in taxa:
        print_str += str(taxon) + sep
    return print_str


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
