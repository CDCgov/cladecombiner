import datetime
from abc import ABC, abstractmethod
from collections.abc import Iterable
from typing import Collection
from warnings import warn

import dendropy

from .nomenclature import HistoryAwareNomenclature, NomenclatureVersioner
from .taxon import Taxon
from .taxon_utils import sort_taxa
from .taxonomy_scheme import PhylogeneticTaxonomyScheme


class Aggregation(dict[Taxon, Taxon]):
    """
    An object for aggregations, basically just a dictionary.
    """

    def _validate(
        self, input_taxa: Iterable[Taxon], taxon_map: dict[Taxon, Taxon]
    ):
        """
        Checks that all input taxa are in the mapping.
        """
        if set(taxon_map.keys()) != set(input_taxa):
            raise RuntimeError(
                "Mismatch between aggregated taxa and input taxa. Input taxa are: "
                + str(input_taxa)
                + " but aggregated taxa are "
                + str(taxon_map.keys())
            )

    def __init__(
        self, input_taxa: Iterable[Taxon], taxon_map: dict[Taxon, Taxon]
    ):
        self._validate(input_taxa, taxon_map)
        super().__init__(taxon_map)

    def to_str(self):
        """
        Get str : str map of taxa names
        """
        return {k.name: v.name for k, v in self.items()}


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


class AsOfAggregator(Aggregator):
    def __init__(
        self,
        taxonomy_scheme: PhylogeneticTaxonomyScheme,
        versioning_provider: HistoryAwareNomenclature,
        as_of: datetime.date,
    ):
        self.as_of = as_of
        self.taxonomy_scheme = taxonomy_scheme
        self.versioner = versioning_provider.get_versioner(as_of)
        self.targets = None
        """
        AsOfAggregator constructor.

        Note that when using an AsOfAggregator, the resulting aggregated taxa
        will only be tips if the taxon is a tip now and was a tip then. That
        is, if there is an ancestral/internal corresponding taxon available,
        a current tip will be mapped to that taxon.

        Parameters
        ----------
        taxonomy_scheme : PhylogeneticTaxonomyScheme
            The tree which we use to do the mapping.

        versioning_provider : HistoryAwareNomenclature
            A Nomenclature with a .get_versioner(as_of) method that can be used
            to determine whether a taxon was recognized as-of the `as_of` date.

        as_of : datetime.date
            The as-of date. The time in the past which defines the set of
            recognized taxa into which we wish to aggregate the input taxa.
        """

    @staticmethod
    def _get_versioned_taxa(
        node: dendropy.Node,
        versioner: NomenclatureVersioner,
        taxa: list[tuple[str, bool, bool]],
    ):
        """
        Recursively add (taxon name, is currently internal, was internal as-of,) tuple to our list
        """

        children_as_of = [
            child
            for child in node.child_node_iter()
            if child.label != node.label and versioner(child.label)
        ]

        was_internal = False
        if len(children_as_of) > 0:
            for child in children_as_of:
                AsOfAggregator._get_versioned_taxa(child, versioner, taxa)
            was_internal = True

        taxa.append(
            (
                node.label,
                not node.is_leaf(),
                was_internal,
            )
        )

    @staticmethod
    def get_versioned_taxa(
        tree: dendropy.Tree, versioner: NomenclatureVersioner
    ) -> Iterable[tuple[str, bool, bool]]:
        """
        Get a list of (taxon name, is internal, was internal,) tuples for all taxa in the
        current tree which were recognized on the as-of date.

        Parameters
        ----------
        tree : dendropy.Tree
            The tree which we use to do the mapping. Should come from a
            PhylogeneticTaxonomyScheme.

        versioner : NomenclatureVersioner
            Used to determine if a name was recognized on the as-of date.

        Returns
        ----------
        Iterable[tuple[str, bool, bool]]
            For each taxon that was recognized on the as-of date, its name,
            whether it has an internal node in the tree now, and whether it had
            an internal node in the tree then. Per cladecombiner style, a taxon
            name can belong to both an ancestral taxon and a tip. For such taxa,
            only the ancestor will be recorded in this list.
        """
        taxa = []
        root = tree.seed_node
        assert root is not None
        assert versioner(root.label)
        AsOfAggregator._get_versioned_taxa(root, versioner, taxa)
        assert (
            len(taxa) > 0
        ), "Found no ancestors of input taxa for given as-of date."
        return taxa

    def aggregate(self, input_taxa: Iterable[Taxon]) -> Aggregation:
        if self.targets is None:
            self.targets = self.get_targets()
        agg = BasicPhylogeneticAggregator(
            self.targets,
            self.taxonomy_scheme,
            sort_clades=True,
            off_target="self",
            warn=False,
        ).aggregate(input_taxa)
        assert all(self.versioner(taxon.name) for taxon in agg.values())
        return agg

    def get_targets(self) -> Collection[Taxon]:
        taxa_as_of = AsOfAggregator.get_versioned_taxa(
            self.taxonomy_scheme.tree, self.versioner
        )
        return [
            Taxon(name, not is_internal) for name, is_internal, _ in taxa_as_of
        ]


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
