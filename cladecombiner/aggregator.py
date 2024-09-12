from abc import ABC, abstractmethod
from collections.abc import Iterable
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
    """
    A partially-abstract base class for aggregation.

    Descendants need only implement next_agg_step() to create a usable aggregation object.
    """

    def __init__(
        self,
        in_taxa: Iterable[Taxon],
    ):
        """
        Constructor for aggregation infrastructure shared by all subclasses.

        Parameters
        ----------
        in_taxa : Iterable[Taxon]
            The taxa we wish to aggregate.
        """
        self.input_taxa = list(in_taxa)
        "The taxa we wish to aggregate."
        self.stack = copy(self.input_taxa)
        "Taxa that still need to be aggregated."
        self.messy_map = {}
        "A record of the entire history of aggregation of all taxa."

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
            step = self.next_agg_step()
            self._valid_agg_step(step)
            self.resolve_agg_step(step)

        self._valid_mapping()

    @abstractmethod
    def next_agg_step(self) -> TaxonMapping:
        """
        Find the next available step in the aggregation process.
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


class FixedAggregator(Aggregator):
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

    def next_agg_step(self) -> TaxonMapping:
        taxon_from = self.stack[0]
        taxon_to = Taxon(self.fixed_map[taxon_from.name], self.make_tip)
        return TaxonMapping({taxon_from: taxon_to}, {taxon_to: True})


class BasicPhylogeneticAggregator(Aggregator):
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

    def _valid_agg_step(self, mapping: TaxonMapping):
        super()._valid_agg_step(mapping)
        unknown_taxa = [
            taxon
            for taxon in mapping.done.keys()
            if not self.taxon_is_valid(taxon)
        ]
        if len(unknown_taxa) > 0:
            raise RuntimeError(
                f"Aggregation step resulted in the following taxa unknown to the provided taxonomy scheme: {unknown_taxa}"
            )

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

    def next_agg_step(self) -> TaxonMapping:
        if len(self.targets) > 0:
            target = self.targets.pop()
            children = self.taxonomy_scheme.descendants(target, True)
            taxa = [taxon for taxon in self.stack if taxon in children]
        elif self.agg_other:
            target = Taxon("other", False)
            taxa = self.stack
        else:
            target = self.stack[-1]
            taxa = [target]

        map = {}
        for taxon in taxa:
            map[taxon] = target
        done = {target: True}
        return TaxonMapping(map, done)

    def taxon_is_valid(self, taxon: Taxon) -> bool:
        if self.taxonomy_scheme.is_valid_taxon(taxon):
            return True
        elif self.agg_other and taxon == Taxon("other", False):
            return True
        return False
