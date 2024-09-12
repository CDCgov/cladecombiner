from abc import ABC, abstractmethod
from collections.abc import Iterable, MutableSequence
from copy import copy
from dataclasses import dataclass
from warnings import warn

from .ordered_taxon import OrderedTaxon
from .taxon import Taxon
from .taxonomy_scheme import PhylogeneticTaxonomyScheme
from .utils import table


@dataclass
class TaxonMapping:
    """
    Class for tracking proposed step in aggregation process.
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

    def _valid_agg_step(self, mapping: TaxonMapping):
        """
        Check that this proposed merge move is acceptable
        """
        unknown = [k for k in mapping.map.keys() if k not in self.stack]
        if len(unknown) > 0:
            raise RuntimeError(
                f"Tried to map {unknown} which is not in current stack {self.stack}."
            )

    def _valid_mapping(self):
        """
        Checks that all input taxa have been mapped exactly once.
        """
        in_str = [taxon.name for taxon in self.input_taxa]
        map = self.map()
        if set(map.keys()) != set(in_str):
            raise RuntimeError(
                "Mismatch between mapped taxa and input taxa. Input taxa are: "
                + str(in_str)
                + " but mapped taxa are "
                + str(map.keys())
            )
        tab = table(map)
        if not all(v == 1 for v in tab.values()):
            raise RuntimeError(
                "Found following taxa mapped more than once: "
                + str([k for k, v in tab.items() if v > 1])
            )

    def aggregate(self):
        """
        Perform aggregation. When done, results may be collected with self.map().

        In the interim, self.messy_map is populated with the results of aggregation.
        """
        while self.stack:
            taxa_in = self.next_taxa()
            step = self.next_agg_step(taxa_in)
            self._valid_agg_step(step)
            self.resolve_agg_step(step)

        self._valid_mapping()

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

    def map(self) -> dict[str, str]:
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
    A mostly-instantiated Aggregator. Descendants need only implement next_agg_step() and next_taxa().
    """

    def __init__(
        self,
        in_taxa: Iterable[Taxon],
    ):
        """
        BoilerplateAggregator constructor.

        This class is partially-abstract and should not be used directly.

        Parameters
        ----------
        in_taxa : Iterable[Taxon]
            The taxa we wish to aggregate.
        """
        self._input_taxa = in_taxa
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
    """
    Aggregation via a user-provided dictionary.
    """

    def __init__(
        self,
        unaggregated: Iterable[str],
        map: dict[str, str],
        create_tip_taxa: bool = False,
    ):
        """
        FixedAggregator constructor.

        Parameters
        ----------
        unaggregated : Iterable[str]
            The taxa we wish to aggregate, as strings.

        map : dict[str, str]
            Dictionary mapping the input taxa to their final units.

        create_tip_taxa : bool
            Should the taxa we create be labeled as tips or not?

        Returns
        -------
        Taxon
            The taxon's parent, or None if this is the root.
        """
        super().__init__(
            in_taxa=[Taxon(nm, False) for nm in unaggregated],
        )
        self.fixed_map = map
        self.make_tip = create_tip_taxa

    def next_agg_step(self, taxa: Iterable[Taxon]) -> TaxonMapping:
        map = {}
        done = {}
        for taxon_from in taxa:
            taxon_to = Taxon(self.fixed_map[taxon_from.name], self.make_tip)
            map[taxon_from] = taxon_to
            done[taxon_to] = True

        return TaxonMapping(map, done)

    def next_taxa(self) -> Iterable[Taxon]:
        return [self.stack[0]]


class BasicPhylogeneticAggregator(BoilerplateAggregator):
    """
    An aggregator which maps a set of input taxa to a fixed set of aggregation targets using a tree.
    """

    def __init__(
        self,
        unaggregated: Iterable[Taxon],
        targets: Iterable[Taxon],
        taxonomy_scheme: PhylogeneticTaxonomyScheme,
        sort_clades: bool = True,
        unmapped_are_other: bool = True,
        warn: bool = True,
    ):
        """
        BasicPhylogeneticAggregator constructor.

        Parameters
        ----------
        unaggregated : Iterable[Taxon]
            The taxa we wish to aggregate.

        targets : Iterable[Taxon]
            The taxa into which we wish to map the unaggregated taxa.

        taxonomy_scheme : PhylogeneticTaxonomyScheme
            The tree which we use to do the mapping.

        sort_clades : bool
            If False, mapping is done using the taxa as ordered in `targets`.
            If True, `targets` are taxonomically sorted so that so that larger
            `targets` do not override smaller ones. For example, if BA.2 and
            BA.2.86 are both aggregation targets, sort_clades = True would handle
            BA.2.86 first, such that JN.1 would map to BA.2.86, while BG.1 would
            map to BA.2. If BA.2 is processed first, both will map to it.

        unmapped_are_other : bool
            If True, any taxon found in the input list which does not fall into
            a target will be mapped to Taxon("other"). If False, any such taxa
            will be mapped to themselves.
        """
        super().__init__(in_taxa=[taxon for taxon in unaggregated])
        self._cached_target = None
        self._cached_taxa = None
        self.targets = [taxon for taxon in targets]
        self._input_targets = copy(self.targets)
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
        self.agg_other = unmapped_are_other
        self.warn = warn

    def next_agg_step(self, taxa: Iterable[Taxon]) -> TaxonMapping:
        assert isinstance(self._cached_target, Taxon)
        assert taxa is self._cached_taxa
        map = {}
        for taxon in taxa:
            map[taxon] = self._cached_target
        done = {self._cached_target: True}
        return TaxonMapping(map, done)

    def next_taxa(self) -> Iterable[Taxon]:
        if len(self.targets) > 0:
            self._cached_target = self.targets.pop()
            children = self.taxonomy_scheme.descendants(
                self._cached_target, True
            )
            self._cached_taxa = [
                taxon for taxon in self.stack if taxon in children
            ]
        elif self.agg_other:
            self._cached_target = Taxon("other", False)
            self._cached_taxa = self.stack
        else:
            self._cached_target = self.stack[-1]
            self._cached_taxa = [self._cached_target]
        return self._cached_taxa

    def _valid_mapping(self):
        super()._valid_mapping()
        if self.warn:
            used_targets = set(self.map().values())
            unused_targets = [
                target.name
                for target in self._input_targets
                if target.name not in used_targets
            ]
            if len(unused_targets) > 0:
                warn(
                    f"The aggregation does not make use of the following input targets: {unused_targets}."
                )
