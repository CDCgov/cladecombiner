from abc import ABC, abstractmethod
from collections.abc import Iterable
from copy import copy
from dataclasses import dataclass
from warnings import warn

from .ordered_taxon import OrderedTaxon
from .taxon import Taxon
from .taxonomy_scheme import PhylogeneticTaxonomyScheme
from .utils import table


class TaxonMapping:
    def __init__(self, input_taxa: Iterable[Taxon], map: dict = None):
        self.input_taxa = input_taxa

        if self.map is None:
            self.map = {taxon: None for taxon in self.input_taxa}
        else:
            self.map = map

    def validate(self):
        """Validate that the mapping is completed"""
        assert isinstance(self.map, dict)
        for key in self.map.keys():
            assert isinstance(key, Taxon)
        for value in self.map.values():
            assert isinstance(value, Taxon)


class Aggregator(ABC):
    """
    A partially-abstract base class for aggregation.

    Descendants need only implement next_agg_step() to create a usable aggregation object.
    """

    def __init__(
        self,
        input_taxa: Iterable[Taxon],
    ):
        """
        Constructor for aggregation infrastructure shared by all subclasses.

        Parameters
        ----------
        input_taxa : Iterable[Taxon]
            The taxa we wish to aggregate.
        """
        self.input_taxa = list(input_taxa)
        "The taxa we wish to aggregate."

    def aggregate(self) -> TaxonMapping:
        pass


class TrivialAggregator(Aggregator):
    """Map every taxon to itself; i.e., do nothing"""

    def __init__(self, input_taxa: Iterable[Taxon]):
        super().__init__(input_taxa)

    def aggregate(self):
        tm = TaxonMapping(
            input_taxa=self.input_taxa, map={x: x for x in self.input_taxa}
        )
        tm.validate()
        return tm


class FixedAggregator(Aggregator):
    """
    Aggregation via a user-provided dictionary.
    """

    def __init__(self, input_taxa: Iterable[Taxon], map: dict[Taxon, Taxon]):
        super().__init__(input_taxa)
        self.map = map

    def aggregate(self) -> TaxonMapping:
        tm = TaxonMapping(self.input_taxa, self.map)
        tm.validate()
        return tm


class BasicPhylogeneticAggregator(Aggregator):
    """
    An aggregator which maps a set of input taxa to a fixed set of aggregation targets using a tree.
    """

    def __init__(
        self,
        input_taxa: Iterable[Taxon],
        target_taxa: Iterable[Taxon],
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
        super().__init__(input_taxa)

        if sort_clades:
            self.target_taxa = [
                otaxon.to_taxon
                for otaxon in sorted(
                    [
                        OrderedTaxon("", False, taxonomy_scheme).from_taxon(
                            taxon, taxonomy_scheme
                        )
                        for taxon in target_taxa
                    ]
                )
            ]
        else:
            self.target_taxa = list(target_taxa)

        # So .pop() executes these in the expected order
        self.targets.reverse()
        self.taxonomy_scheme = taxonomy_scheme
        self.agg_other = unmapped_are_other
        self.warn = warn

    def aggregate(self):
        # initialize the aggregator state: initially, no taxa are "done"
        # (i.e., terminally mapped)
        done_taxa = []
        stack = [taxon for taxon in self.input_taxa if taxon not in done_taxa]
        # I'm not surehow to initialize the mapping
        step = None
        while self.stack:
            # confusingly, the "step" is actually a map. I think it would be
            # easier if .next_step() took in the current map and gave an updated
            # one
            step = self.next_step(step)
            # then we want to check that that updated mapping looks good
            self.validate_step(step)
            # not sure we need this, in that paradigm
            self.resolve_step(step)

        # the step is actually a TaxonMapping, which should now be complete
        step.validate()

        return step

    @staticmethod
    def validate_step(mapping: TaxonMapping):
        # I'm not sure if, at every step here, the mapping is "complete" in the
        # sense that we could validly stop the aggregation at that point. So we
        # might need some different kind of validation here.
        mapping.validate()

        unknown_taxa = [
            taxon for taxon in mapping.done.keys() if not self.taxon_is_valid(taxon)
        ]
        if len(unknown_taxa) > 0:
            raise RuntimeError(
                f"Aggregation step resulted in the following taxa unknown to the provided taxonomy scheme: {unknown_taxa}"
            )

        # I'm not sure if these belong together!
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
