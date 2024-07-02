from collections.abc import Sequence
from os import path
from typing import Optional

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
