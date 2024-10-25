from abc import ABC, abstractmethod
from collections.abc import Iterable

from ..taxon import Taxon
from ..taxonomy_scheme import PhylogeneticTaxonomyScheme
from ..utils import validate_kwargs


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


class EnumerateSteps(ABC):
    @abstractmethod
    def __call__(self, agg: Aggregation, **kwargs) -> Iterable[Aggregation]:
        raise NotImplementedError()


class Continue(ABC):
    @abstractmethod
    def __call__(self, agg: Aggregation, **kwargs) -> bool:
        raise NotImplementedError()


class EnumerateAgglomerative(EnumerateSteps):
    exp_kwargs = [
        ("pts", PhylogeneticTaxonomyScheme, True),
        ("max_polyphyly_cuts", int, False),
    ]

    def __call__(self, agg: Aggregation, **kwargs):
        validate_kwargs(EnumerateAgglomerative.exp_kwargs, **kwargs)

        pts: PhylogeneticTaxonomyScheme = kwargs["pts"]
        max_polyphyly = kwargs.get("max_polyphyly_cuts", 0)

        # @TODO check logic for handling of root
        active_taxa = set(agg.values())
        if pts.root() in [pts.taxon_to_node[taxon] for taxon in active_taxa]:
            return []

        active_parents = set([pts.parents(taxon) for taxon in active_taxa])
        # We don't want to "aggregate" t1 to any of the "+" nodes
        # /---+---+---+--- t1
        # @
        # \--- t2
        # This is in-progress PhylogeneticTaxonomyScheme.MRBA()

        # If t1, t2, and t3 are monophyletic, and we want a monophyletic grouping,
        # we can only aggregate to node "@" not node "*"
        #    /--- t1
        #    |
        # ---@
        #    |   /--- t2
        #    \---*
        #        \--- t3
        # This should be offloaded to PhylogeneticTaxonomyScheme
        remove = set()
        for parent in active_parents:
            assert isinstance(parent, Taxon)
            if pts.parents(parent) in active_parents:
                remove.add(parent)
        active_parents = active_parents.difference(remove)

        next_steps = []
        for parent in active_parents:
            assert isinstance(parent, Taxon)  # pylance paranoia
            children = list(pts.children(parent))
            monophyletic_next = {}
            for utax, atax in agg.items():
                # @TODO is this right, or should all children be expected to be in the map?
                if atax in children:
                    monophyletic_next[utax] = parent
                else:
                    monophyletic_next[utax] = atax
            next_steps.append(Aggregation(agg.keys(), monophyletic_next))

            if max_polyphyly == 0:
                continue

            if max_polyphyly > 1:
                raise NotImplementedError()

            # for node in pts.taxon_to_node[parent].child_node_iter():
            #     paraphyletic_next = deepcopy(agg)
