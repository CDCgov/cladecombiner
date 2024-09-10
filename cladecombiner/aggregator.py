from abc import ABC, abstractmethod
from collections.abc import Iterable, MutableSequence
from copy import copy
from dataclasses import dataclass

from .ordered_taxon import OrderedTaxon
from .taxon import Taxon
from .taxonomy_scheme import PhylogeneticTaxonomyScheme


@dataclass
class TaxonMapping:
    """
    Class for tracking proposed step in aggregation process
    """

    map: dict[Taxon, Taxon]
    "from : to dictionary of the mapping"
    done: dict[Taxon, bool]
    "to : add_to_stack dictionary tracking if the aggregated taxon should be added back to the stack"


class Aggregator(ABC):
    """ """

    @property
    @abstractmethod
    def input_taxa(self) -> Iterable[Taxon]:
        """
        The taxa which are to be aggregated.
        """
        pass

    @property
    @abstractmethod
    def messy_map(self) -> dict[Taxon, Taxon]:
        """
        A record of the entire history of aggregation of all taxa.
        """
        pass

    @property
    @abstractmethod
    def stack(self) -> MutableSequence[Taxon]:
        """
        Taxa that still need to be aggregated.
        """
        pass

    def aggregate(self):
        """
        Perform aggregation. When done, results may be collected with self.map().

        In the interim, self.messy_map is populated with the results of aggregation.
        """
        while self.stack:
            taxa_in = self.next_taxa()
            step = self.next_agg_step(taxa_in)
            self.check_agg_step(step)
            self.resolve_agg_step(step)

    def check_agg_step(self, mapping: TaxonMapping):
        """
        Check that this proposed merge move is acceptable
        """
        assert all(key in self.stack for key in mapping.map)

    @abstractmethod
    def next_taxa(self) -> Iterable[Taxon]:
        """
        Find the next taxon or taxa to merge.
        """
        pass

    @abstractmethod
    def next_agg_step(self, taxa: Iterable[Taxon]) -> TaxonMapping:
        """
        Find the next available step in the aggregation process for these taxa.
        """
        pass

    def map(self):
        """
        Create clean {Taxon from : Taxon to} aggregation dictionary.
        """
        if len(self.stack) > 0:
            self.aggregate()
        taxon_map = {}
        for taxon_from in self.input_taxa:
            taxon_to = taxon_from
            while taxon_to in self.messy_map.keys():
                if taxon_to == self.messy_map[taxon_to]:
                    # Self-mappings are possible
                    break
                taxon_to = self.messy_map[taxon_to]
            taxon_map[taxon_from.name] = taxon_to.name
        return taxon_map

    def resolve_agg_step(self, mapping: TaxonMapping):
        for k, v in mapping.map.items():
            self.messy_map[k] = v
            self.stack.remove(k)

        for k, v in mapping.done.items():
            if not v:
                self.stack.append(k)


class BoilerplateAggregator(Aggregator):
    """
    A mostly-instantiated Aggregator. Descendants need only implement next_agg_step() and next_taxa()
    """

    def __init__(
        self,
        in_taxa: Iterable[Taxon],
    ):
        self._input_taxa = in_taxa
        self._last_target = Taxon("", False)
        self._stack = list(copy(self._input_taxa))
        self._messy_map = {}

    @property
    def input_taxa(self):
        return self._input_taxa

    @property
    def messy_map(self) -> dict[Taxon, Taxon]:
        return self._messy_map

    @property
    def stack(self) -> MutableSequence[Taxon]:
        return self._stack


class FixedAggregator(BoilerplateAggregator):
    def __init__(self, names_to_aggregate: Iterable[str], map: dict[str, str]):
        super().__init__(
            in_taxa=[Taxon(nm, False) for nm in names_to_aggregate],
        )
        self.fixed_map = map

    def next_agg_step(self, taxa: Iterable[Taxon]) -> TaxonMapping:
        map = {}
        done = {}
        for taxon_from in taxa:
            taxon_to = Taxon(self.fixed_map[taxon_from.name], False)
            map[taxon_from] = taxon_to
            done[taxon_to] = True

        return TaxonMapping(map, done)

    def next_taxa(self) -> Iterable[Taxon]:
        return [self.stack[0]]


class BasicPhylogeneticAggregator(BoilerplateAggregator):
    def __init__(
        self,
        unaggregated: Iterable[Taxon],
        targets: Iterable[Taxon],
        taxonomy_scheme: PhylogeneticTaxonomyScheme,
        sort_clades: bool = True,
    ):
        super().__init__(in_taxa=[taxon for taxon in unaggregated])
        self.cached_target = None
        self.cached_taxa = None
        self.targets = [taxon for taxon in targets]
        if sort_clades:
            ordered_taxa = [
                OrderedTaxon("", False, taxonomy_scheme).from_taxon(
                    taxon, taxonomy_scheme
                )
                for taxon in targets
            ]
            ordered_taxa.sort()
            self.targets = [otaxon.to_taxon() for otaxon in ordered_taxa]
        # So .pop() executes these in the expected order
        self.targets.reverse()
        self.taxonomy_scheme = taxonomy_scheme

    def next_agg_step(self, taxa: Iterable[Taxon]) -> TaxonMapping:
        assert isinstance(self.cached_target, Taxon)
        assert taxa is self.cached_taxa
        map = {}
        for taxon in taxa:
            map[taxon] = self.cached_target
        done = {self.cached_target: True}
        return TaxonMapping(map, done)

    def next_taxa(self) -> Iterable[Taxon]:
        if len(self.targets) > 0:
            self.cached_target = self.targets.pop()
            children = self.taxonomy_scheme.descendants(
                self.cached_target, True
            )
            self.cached_taxa = [
                taxon for taxon in self.stack if taxon in children
            ]
        else:
            self.cached_target = Taxon("Other", False)
            self.cached_taxa = self.stack
        return self.cached_taxa
