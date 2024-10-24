import math
from copy import deepcopy
from typing import Iterable

import dendropy

from ..taxon import Taxon
from ..taxonomy_scheme import PhylogeneticTaxonomyScheme
from .aggregator import Aggregator


def score_clade_count(
    node: dendropy.Node, pts: PhylogeneticTaxonomyScheme, const: float = 1.0
) -> float:
    """
    Score of a subtree is -(const * # clades contained in subtree)
    """
    # Super inefficient
    return -const * (1 + len(pts.descendants(pts.node_to_taxon[node], False)))


def score_data_entropy(
    node: dendropy.Node,
    pts: PhylogeneticTaxonomyScheme,
    const: float = 1.0,
    pow: float = 10.0,
) -> float:
    """ """
    x = pts.node_to_taxon[node].data
    return math.pow(-const * x * math.log(x), pow)


def score_data_pow(
    node: dendropy.Node,
    pts: PhylogeneticTaxonomyScheme,
    const: float = 1.0,
    pow: float = 10.0,
) -> float:
    """ """
    return const * math.pow(pts.node_to_taxon[node].data, pow)


def score_data_log(
    node: dendropy.Node,
    pts: PhylogeneticTaxonomyScheme,
    const: float = 1.0,
    pow: float = 5.0,
) -> float:
    """ """
    return const * math.pow(math.log(pts.node_to_taxon[node].data), pow)


class MonoAutoAgg(Aggregator):
    def __init__(self, taxonomy_scheme: PhylogeneticTaxonomyScheme):
        self.intact_pts = taxonomy_scheme
        self.mutable_pts = deepcopy(self.intact_pts)
        self.combine_data = sum  # some callable on a list
        self.combine_scores = sum
        self.score_subtree = score_clade_count
        self.subtree_score_const = 0
        self.data_score_const = 0
        self.score_data = score_data_entropy

    def aggregate(self, taxa: Iterable[Taxon]):
        self.mutable_pts = deepcopy(self.intact_pts)
        assert set(self.mutable_pts.taxon_to_node.keys()).issuperset(set(taxa))
        for node in self.mutable_pts.tree.leaf_node_iter():
            if self.mutable_pts.node_to_taxon[node] not in taxa:
                self.mutable_pts.prune_subtree(
                    self.mutable_pts.node_to_taxon[node]
                )
        root = self.mutable_pts.tree.seed_node
        assert isinstance(root, dendropy.Node)
        self.subtree_score_const = abs(
            1.0 / self.score_subtree(root, self.mutable_pts)
        )

        self.prepare_data(taxa)
        # self.data_score_const = abs(1.0 / self.score_data(root, self.mutable_pts))
        mean_data_score = sum(
            [
                self.score_data(node, self.mutable_pts)
                for node in self.mutable_pts.node_to_taxon.keys()
            ]
        ) / len(self.mutable_pts.taxon_to_node.keys())
        print(f">>>>>>>> mean data score = {mean_data_score} <<<<<<<<")
        self.data_score_const = abs(1.0 / (mean_data_score))
        return self.ttr_greedy(self.score_nodes())

    def combine_node_data(self, node: dendropy.Node) -> None:
        child_data = []
        for child in node.child_node_iter():
            if self.mutable_pts.node_to_taxon[child].data is None:
                self.combine_node_data(child)
            child_data.append(self.mutable_pts.node_to_taxon[child].data)
        self.mutable_pts.node_to_taxon[node].data = self.combine_data(
            child_data
        )

    def prepare_data(self, taxa: Iterable[Taxon]):
        tree_taxa = list(self.mutable_pts.taxon_to_node.keys())
        for taxon in taxa:
            if taxon.data is None:
                raise RuntimeError(f"Input taxon {taxon} has no data.")
            tree_taxon = tree_taxa[tree_taxa.index(taxon)]
            # TODO Should probably check more than just not none?
            tree_taxon.data = taxon.data
        assert all(
            self.mutable_pts.node_to_taxon[leaf].data is not None
            for leaf in self.mutable_pts.tree.leaf_node_iter()
        )
        for node in self.mutable_pts.tree.preorder_internal_node_iter():
            self.mutable_pts.node_to_taxon[node].data = None
        assert isinstance(self.mutable_pts.tree.seed_node, dendropy.Node)
        self.combine_node_data(self.mutable_pts.tree.seed_node)

    def score_nodes(self) -> dict[Taxon, float]:
        scores: dict[Taxon, float] = {}
        for node in self.mutable_pts.tree.postorder_node_iter():
            taxon = self.mutable_pts.node_to_taxon[node]
            data_score = self.score_data(
                node, self.mutable_pts, self.data_score_const
            )
            subtree_score = self.score_subtree(
                node, self.mutable_pts, self.subtree_score_const
            )
            scores[taxon] = data_score + subtree_score
            print(
                f"++ Scoring taxon {taxon}: {scores[taxon]} = {data_score} + {subtree_score}"
            )
        return scores

    def ttr_greedy(self, scores: dict[Taxon, float]) -> Iterable[Taxon]:
        use_flags: dict[dendropy.Node, bool] = {}
        for node in self.mutable_pts.tree.leaf_node_iter():
            use_flags[node] = True
        for node in self.mutable_pts.tree.postorder_internal_node_iter():
            if all(use_flags[child] for child in node.child_node_iter()):
                child_score_tot = 0.0
                for child in node.child_node_iter():
                    child_score_tot += scores[
                        self.mutable_pts.node_to_taxon[child]
                    ]
                if not (
                    scores[self.mutable_pts.node_to_taxon[node]]
                    < child_score_tot
                ):  # a hack because math.inf is neither > or < itself
                    use_flags[node] = True
                    for child in node.child_node_iter():
                        use_flags[child] = False
                else:
                    use_flags[node] = False
            else:
                use_flags[node] = False
        return [
            taxon
            for taxon in self.mutable_pts.taxon_to_node.keys()
            if use_flags[self.mutable_pts.taxon_to_node[taxon]]
        ]
