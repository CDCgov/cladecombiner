from collections.abc import Iterable, Sequence
from functools import cmp_to_key
from os import path
from typing import Optional

from .nomenclature import Nomenclature
from .taxon import Taxon
from .taxonomy_scheme import TaxonomyScheme, TreelikeTaxonomyScheme


def read_taxa(
    fp: str,
    is_tip: bool | Sequence[bool] = True,
    nomenclature: Optional[Nomenclature] = None,
    taxonomy_scheme: Optional[TaxonomyScheme] = None,
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
    assert nomenclature is None or isinstance(nomenclature, Nomenclature)

    assert taxonomy_scheme is None or isinstance(
        taxonomy_scheme, TaxonomyScheme
    )

    ext = path.splitext(fp)[1]
    taxa = []
    if ext == ".txt":
        f = open(fp)
        lines = f.readlines()
        f.close()
        taxa = []
        if not isinstance(is_tip, Sequence):
            is_tip = [is_tip for _ in range(len(lines))]

        for i in range(len(lines)):
            taxon = Taxon(lines[i][:-1], is_tip[i])
            taxa.append(taxon)

    if nomenclature:
        nomenclature.validate([taxon.name for taxon in taxa])

    if taxonomy_scheme:
        taxonomy_scheme.validate([taxon for taxon in taxa])
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


def sort_taxa(
    taxa: Iterable[Taxon], taxonomy_scheme: TreelikeTaxonomyScheme
) -> list[Taxon]:
    """
    Sorts taxa into a phylogenetic preorder according to a taxonomy scheme, such that if
    taxon X contains Y, Y comes before X.

    For example, the Pango lineages [KP.1, JN.1, BA.2, BA.3] will be sorted such that (1) KP.1
    comes before both JN.1 and BA.2 and (2) JN.1 appears before BA.2. The ordering of these
    with respect to BA.3 is arbitrary, as BA.3 is sister to BA.2 and its descendants.

    Parameters
    ---------
    taxa : Iterable[Taxon]
        The Taxon objects to be sorted.
    taxonomy_scheme : TreelikeTaxonomyScheme
        The taxonomy scheme by which to sort the taxa.

    Returns
    -------
    list[Taxon]
        The sorted taxa.
    """
    taxonomy_scheme.validate(taxa)
    return sorted(
        taxa,
        key=cmp_to_key(
            lambda x, y: 1 if taxonomy_scheme.contains(x, y) else -1
        ),
    )
