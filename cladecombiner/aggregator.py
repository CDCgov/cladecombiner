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
    def aggregate(self, input_taxa: Iterable[Taxon]) -> TaxonMapping:
        pass


class TrivialAggregator(Aggregator):
    """Map every taxon to itself; i.e., do nothing"""

    def aggregate(self, input_taxa: Iterable[Taxon]) -> TaxonMapping:
        tm = TaxonMapping(
            input_taxa=self.input_taxa, map={x: x for x in self.input_taxa}
        )
        tm.validate()
        return tm


class FixedAggregator(Aggregator):
    """Aggregation via a user-provided dictionary"""

    def __init__(self, map: dict[Taxon, Taxon]):
        self.map = map

    def aggregate(self, input_taxa: Iterable[Taxon]) -> TaxonMapping:
        tm = TaxonMapping(self.input_taxa, self.map)
        tm.validate()
        return tm


class BasicPhylogeneticAggregator(Aggregator):
    """
    An aggregator which maps a set of input taxa to a fixed set of aggregation targets using a tree.
    """

    def __init__(
        self,
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
        self.target_taxa.reverse()
        self.taxonomy_scheme = taxonomy_scheme
        self.unmapped_are_other = unmapped_are_other
        self.warn = warn

    def aggregate(self, input_taxa: Iterable[Taxon]) -> TaxonMapping:
        # initialize the state; not sure what this should be
        mapping = None
        # stack are the incomplete taxa
        stack = None

        # iterate
        while stack:
            # update the mapping
            mapping = self.step(mapping)
            # validate the mapping
            self.validate_step(mapping)

        # validate the final mapping
        self.validate_complete()

        # the mapping should be complete
        mapping.validate()

        return mapping

    @staticmethod
    def validate_step(mapping: TaxonMapping):
        # check for unknown output taxa
        raise NotImplementedError

    @staticmethod
    def validate_complete(mapping: TaxonMapping):
        # check that all the input taxa were used?
        raise NotImplementedError

    def step(self, incomplete_targets, mapping: TaxonMapping) -> TaxonMapping:
        if len(self.target_taxa) > 0:
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


class SerialAggregator(Aggregator):
    def __init__(self, aggregators: Iterable[Aggregator]):
        self.aggregators = aggregators

    def aggregate(self, input_taxa: Iterable[Taxon]) -> TaxonMapping:
        # initialize the input taxa
        input_taxa = list(input_taxa)
        composite_mapping = {taxon: None for taxon in input_taxa}

        for aggregator in self.aggregators:
            # get the next mapping
            mapping = aggregator.aggregate(input_taxa)
            # input for the next step are targets from this step
            input_taxa = list(mapping.map.values())
            # update the composite mapping: say the composite map has A -> B
            # and this next mapping has B -> C, so we update the composite map
            # to be A -> C
            composite_mapping = TaxonMapping(
                {key: mapping.map[value] for key, value in composite_mapping.items()}
            )

        composite_mapping.validate()

        return composite_mapping
