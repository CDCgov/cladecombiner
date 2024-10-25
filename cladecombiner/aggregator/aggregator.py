from abc import ABC, abstractmethod
from collections.abc import Iterable
from warnings import warn

from ..taxon import Taxon
from ..taxon_utils import sort_taxa
from ..taxonomy_scheme import PhylogeneticTaxonomyScheme
from .utils import Aggregation


class Aggregator(ABC):
    """
    Aggregators return Aggregations, maps of input_taxon : aggregated_taxon
    """

    @abstractmethod
    def aggregate(self, input_taxa: Iterable[Taxon]) -> Aggregation:
        raise NotImplementedError()


class ArbitraryAggregator(Aggregator):
    """
    Aggregation via a user-provided dictionary.
    """

    def __init__(
        self,
        map: dict[Taxon, Taxon],
    ):
        """
        FixedAggregator constructor.

        Parameters
        ----------
        map : dict[Taxon, Taxon]
            Dictionary mapping the input taxa to their aggregated taxa.
        """
        self.map = map

    def aggregate(self, input_taxa: Iterable[Taxon]) -> Aggregation:
        return Aggregation(
            input_taxa, {taxon: self.map[taxon] for taxon in input_taxa}
        )


class AutomaticAggregator(Aggregator):
    """ "
    Aggregation based on objective criteria.
    """


class BasicPhylogeneticAggregator(Aggregator):
    """
    An aggregator which maps a set of input taxa to a fixed set of aggregation targets using a tree.
    """

    def __init__(
        self,
        targets: Iterable[Taxon],
        taxonomy_scheme: PhylogeneticTaxonomyScheme,
        sort_clades: bool = True,
        off_target: str = "other",
        warn: bool = True,
    ):
        """
        BasicPhylogeneticAggregator constructor.

        Parameters
        ----------
        targets : Iterable[Taxon]
            The taxa into which we wish to aggregate the input taxa.

        taxonomy_scheme : PhylogeneticTaxonomyScheme
            The tree which we use to do the mapping.

        sort_clades : bool
            If False, mapping is done using the taxa as ordered in `targets`.
            If True, `targets` are taxonomically sorted so that so that larger
            `targets` do not override smaller ones. For example, if BA.2 and
            BA.2.86 are both aggregation targets, sort_clades = True would handle
            BA.2.86 first, such that JN.1 would map to BA.2.86, while BG.1 would
            map to BA.2. If BA.2 is processed first, both will map to it.

        off_target : str
            Specifies what to do with taxa which do not belong to any target.
            Options are "other" for aggregating all such taxa into Taxon("other"),
            and "self" for aggregating all such taxa into themselves.
        """
        self.taxonomy_scheme = taxonomy_scheme
        self.taxonomy_scheme.validate(targets)
        self.targets = [taxon for taxon in targets]
        off_target_options = ["self", "other"]
        if off_target not in off_target_options:
            raise RuntimeError(
                f"Unrecognized value for `off_target`, options are:{off_target}"
            )
        self.off_target = off_target
        self.warn = warn
        if sort_clades:
            self.targets = sort_taxa(self.targets, self.taxonomy_scheme)

    def _check_missing(self, agg_map: dict[Taxon, Taxon]):
        if self.warn:
            used_targets = set(agg_map.values())
            unused_targets = [
                target for target in self.targets if target not in used_targets
            ]
            if len(unused_targets) > 0:
                warn(
                    f"The aggregation does not make use of the following input targets: {unused_targets}."
                )

    def aggregate(self, input_taxa: Iterable[Taxon]) -> Aggregation:
        self.taxonomy_scheme.validate(input_taxa)
        agg_map: dict[Taxon, Taxon] = {}
        stack = set(input_taxa)
        for target in self.targets:
            children = self.taxonomy_scheme.descendants(target, True)
            sub_map = {taxon: target for taxon in stack if taxon in children}
            agg_map = agg_map | sub_map
            stack.difference_update(set(agg_map.keys()))

        if len(stack) > 0:
            if self.off_target == "other":
                cleanup = HomogenousAggregator(
                    Taxon("other", False)
                ).aggregate(stack)
            else:
                cleanup = SelfAggregator().aggregate(stack)
            agg_map = agg_map | cleanup

        self._check_missing(agg_map)

        return Aggregation(input_taxa, agg_map)


class HomogenousAggregator(Aggregator):
    """
    Aggregation of every taxon to some catch-all taxon.
    """

    def __init__(self, taxon: Taxon):
        self.agg_taxon = taxon

    def aggregate(self, input_taxa: Iterable[Taxon]) -> Aggregation:
        return Aggregation(
            input_taxa, {taxon: self.agg_taxon for taxon in input_taxa}
        )


class SelfAggregator(Aggregator):
    """
    Aggregation of every taxon to itself
    """

    def __init__(self):
        pass

    def aggregate(self, input_taxa: Iterable[Taxon]) -> Aggregation:
        return Aggregation(input_taxa, {taxon: taxon for taxon in input_taxa})


class SerialAggregator(Aggregator):
    """
    A number of aggregators chained in serial.
    """

    def __init__(self, aggregators: Iterable[Aggregator]):
        self.aggregators = aggregators

    def aggregate(self, input_taxa: Iterable[Taxon]) -> Aggregation:
        taxa = list(input_taxa)
        comp_agg = SelfAggregator().aggregate(input_taxa)

        for aggregator in self.aggregators:
            agg = aggregator.aggregate(taxa)
            taxa = set(agg.values())
            comp_agg = {taxon: agg[comp_agg[taxon]] for taxon in input_taxa}

        return Aggregation(input_taxa, comp_agg)
