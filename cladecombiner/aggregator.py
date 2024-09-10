from abc import ABC, abstractmethod
from collections.abc import Iterable, MutableSequence
from copy import copy
from dataclasses import dataclass

from .taxon import Taxon


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
            taxon_map[taxon_from] = taxon_to
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
            print("++++++ Mapping " + str(taxon_from) + " to " + str(taxon_to))
            map[taxon_from] = taxon_to
            done[taxon_to] = True

        return TaxonMapping(map, done)

    def next_taxa(self) -> Iterable[Taxon]:
        return [self.stack[0]]


# class FixedAggregator(Aggregator):

#     def __init__(self, taxa_in: Iterable[str], map: dict[str, str]):
#         self.fixed_map = map
#         self._input_taxa = [Taxon(taxon, False) for taxon in taxa_in]
#         self._last_target = Taxon("", False)
#         self._stack = self._input_taxa
#         self._map = {}
#         self._messy_map = {}

#     @property
#     def input_taxa(self):
#         return self._input_taxa

#     @property
#     def last_target(self):
#         return self._last_target

#     @property
#     def stack(self):
#         return self._stack

#     @property
#     def map(self):
#         return self._map

#     @property
#     def messy_map(self):
#         return self._messy_map
