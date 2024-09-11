from functools import total_ordering
from typing import Any

from .taxon import Taxon
from .taxonomy_scheme import TreelikeTaxonomyScheme


@total_ordering
class OrderedTaxon(Taxon):
    """
    A Taxon equipped with a TreelikeTaxonomyScheme which can be sorted.

    For the purposes of comparison, "A < B" means "A is contained within B (but is not B)."
    """

    def __init__(
        self,
        name: str,
        is_tip: bool,
        taxonomy_scheme: TreelikeTaxonomyScheme,
        data: Any = None,
    ):
        super().__init__(name, is_tip, data)
        self.taxonomy_scheme = taxonomy_scheme

    @classmethod
    def from_taxon(cls, taxon: Taxon, taxonomy_scheme: TreelikeTaxonomyScheme):
        return cls(taxon.name, taxon.tip, taxonomy_scheme, taxon.data)

    def to_taxon(self):
        return Taxon(self.name, self.tip, self.data)

    def __lt__(self, other):
        assert (
            hasattr(other, "taxonomy_scheme")
            and self.taxonomy_scheme is other.taxonomy_scheme
        )
        if self == other:
            return False
        if self.taxonomy_scheme.contains(other, self):
            return True
        return False
