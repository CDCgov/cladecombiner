from collections.abc import Sequence
from os import path
from typing import Optional

from .nomenclature import Nomenclature
from .taxon import Taxon
from .taxonomy_scheme import TaxonomyScheme


def read_taxa(
    fp: str,
    nomenclature: Optional[Nomenclature],
    taxonomy_scheme: Optional[TaxonomyScheme],
) -> Sequence[Taxon]:
    """Reads in taxa as a list of Taxon objects."""
    ext = path.splitext(fp)[1]
    taxa = []
    if ext == ".txt":
        f = open(fp)
        lines = f.readlines()
        f.close()
        taxa = []
        for line in lines:
            if nomenclature:
                if not nomenclature.is_valid_name(line[:-1]):
                    raise RuntimeError(
                        "The name "
                        + line[:-1]
                        + " is not valid under the provided nomenclature ("
                        + str(nomenclature)
                        + ")"
                    )
            taxon = Taxon(line[:-1])
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
    """Prettier printing of lists of taxa."""
    print_str = ""
    for taxon in taxa:
        print_str += str(taxon) + sep
    return print_str
