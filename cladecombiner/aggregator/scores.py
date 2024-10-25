import operator
from abc import ABC, abstractmethod
from typing import Any, Iterable

import dendropy
from param import Callable

from cladecombiner.taxon import Taxon
from cladecombiner.utils import validate_kwargs

from .utils import Aggregation


class AggregationScorer(ABC):
    """
    An AggregationScorer scores an aggregation. Higher scores are better aggregations.
    """

    def __init__(self, **kwargs):
        raise NotImplementedError()

    @abstractmethod
    def score(self, agg: Aggregation, **kwargs) -> float:
        raise NotImplementedError()

    @abstractmethod
    def process_for_components(
        self, agg: Aggregation, **kwargs
    ) -> dict["str", Any]:
        raise NotImplementedError()


class DataScorer(ABC):
    """
    A DataScorer scores the part of an aggregation which has to do with the data.

    For example, it might enforce that all taxa have taxon.data >= threshold,
    favor larger values of taxon.data, or both.
    """

    def __init__(self, **kwargs):
        raise NotImplementedError()

    @abstractmethod
    def score(self, **kwargs) -> float:
        raise NotImplementedError()


class TaxonScorer(ABC):
    """
    A TaxonScorer scores the part of an aggregation which has to do
    """

    def __init__(self, **_):
        raise NotImplementedError()

    @abstractmethod
    def score(self, **kwargs) -> float:
        raise NotImplementedError()


class DecomposableAggregationScorer(AggregationScorer):
    """
    An AggregationScorer that can be written as
    """

    exp_kwargs = [
        ("taxon_scorer", TaxonScorer, True),
        ("data_scorer", DataScorer, True),
        ("fun", Callable, True),
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        validate_kwargs(DecomposableAggregationScorer.exp_kwargs, **kwargs)
        self.taxon_scorer = kwargs["taxon_scorer"]
        self.data_scorer = kwargs["data_scorer"]
        self.score_combiner = kwargs["fun"]

    def score(self, agg: Aggregation, **kwargs) -> float:
        kwargs = self.process_for_components(agg, **kwargs)
        return self.score_combiner(
            self.taxon_scorer(**kwargs), self.data_scorer(**kwargs)
        )


class AdditiveAggregationScorer(DecomposableAggregationScorer):
    """
    Convenience shortcut class for when the combined score is the sum of the taxon and data scores.
    """

    def __init__(self, **kwargs):
        kwargs["fun"] = operator.add
        super().__init__(**kwargs)


class PolyphyleticStarScorer(AggregationScorer):
    exp_kwargs = [("polyphyletic_taxa", Iterable[Taxon], True)]

    def __init__(self, **kwargs):
        validate_kwargs(PolyphyleticStarScorer.exp_kwargs, **kwargs)
        self.polyphyletic_taxa = kwargs["polyphyletic_taxa"]
        super().__init__(**kwargs)

    def tree_from_kwargs(self, **kwargs) -> dendropy.Tree:
        pass

    # def process_for_components(self, agg: Aggregation, **kwargs) -> dict[]:
    #     tree = self.tree_from_kwargs(**kwargs)

    #     # Get forest from agg and tree
    #     # Return the forest


class PolyphyleticStarDataScorer(DataScorer):
    """
    A DataScorer based on a tree of taxa, where some taxa are removed from the tree and placed into some non-tree-based polyphyletic taxa.

    All polyphyletic taxa are treated as star trees for the purposes of combining data
    """

    exp_kwargs = [
        ("polyphyletic_taxa", Iterable[Taxon], True),
        ("data_combiner", Callable, True),
        ("data_scorer", Callable, True),
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        validate_kwargs(PolyphyleticStarDataScorer.exp_kwargs, **kwargs)
        self.polyphyletic_taxa = kwargs["polyphyletic_taxa"]
        "Names of the taxa that are polyphyletic"
        self.data_combiner = kwargs["data_combiner"]
        "A function that gets a node's data from its children's data"
        self.data_scorer = kwargs["data_scorer"]
        "A function that scores a node's data"
