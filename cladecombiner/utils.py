from collections.abc import Iterable, Sequence
from os import path
from typing import Any, Optional

import dendropy

from .nomenclature import Nomenclature
from .taxon import Taxon
from .taxonomy_scheme import TaxonomyScheme


def fully_labeled_subtrees_same(
    node1: dendropy.Node, node2: dendropy.Node
) -> bool:
    """
    Are two subtrees with every node labeled topologically equivalent?
    """
    children1 = node1.child_nodes()
    children2 = node2.child_nodes()

    child_labels1 = [child.label for child in children1].sort()
    child_labels2 = [child.label for child in children2].sort()

    if not child_labels1 == child_labels2:
        return False

    for child1 in children1:
        for child2 in children2:
            if child1.label == child2.label:
                fully_labeled_subtrees_same(child1, child2)

    return True


def fully_labeled_trees_same(
    tree1: dendropy.Tree, tree2: dendropy.Tree
) -> bool:
    """
    Are two trees with every node labeled topologically equivalent?
    """
    if isinstance(tree1.seed_node, dendropy.Node) and isinstance(
        tree2.seed_node, dendropy.Node
    ):
        if tree1.seed_node.label != tree2.seed_node.label:
            return False
        else:
            return fully_labeled_subtrees_same(
                tree1.seed_node, tree2.seed_node
            )
    else:
        # Should never hit, required for type checking
        raise RuntimeError("Malformed tree, seed_node must be a dendropy.Node")


def read_taxa(
    fp: str,
    is_tip: bool | Sequence[bool],
    nomenclature: Optional[Nomenclature],
    taxonomy_scheme: Optional[TaxonomyScheme],
) -> Sequence[Taxon]:
    """
    Reads in taxa as a list of Taxon objects.
    """
    ext = path.splitext(fp)[1]
    taxa = []
    if ext == ".txt":
        f = open(fp)
        lines = f.readlines()
        f.close()
        taxa = []
        if not isinstance(is_tip, Sequence):
            is_tip = [is_tip for _ in range(len(lines))]

        for line, i in zip(lines, range(len(lines))):
            if nomenclature:
                if not nomenclature.is_valid_name(line[:-1]):
                    raise RuntimeError(
                        "The name "
                        + line[:-1]
                        + " is not valid under the provided nomenclature ("
                        + str(nomenclature)
                        + ")"
                    )
            taxon = Taxon(line[:-1], is_tip[i])
            if taxonomy_scheme:
                if not taxonomy_scheme.is_valid_taxon(taxon):
                    raise RuntimeError(
                        "The name "
                        + str(taxon)
                        + " is not valid under the provided taxonomy scheme ("
                        + str(taxonomy_scheme)
                        + ")"
                    )
            taxa.append(taxon)
    return taxa


def printable_taxon_list(taxa: Sequence[Taxon], sep: str = "\n") -> str:
    """Prettier printing of lists of taxa."""
    print_str = ""
    for taxon in taxa:
        print_str += str(taxon) + sep
    return print_str


def table(x: Iterable) -> dict:
    """Like R's base::table(), counts occurrences of elements in container"""
    unique = set(x)
    res = {}
    for u in unique:
        res[u] = 0
    for item in x:
        res[item] += 1
    return res


def table_equal(x1: Iterable, x2: Iterable) -> bool:
    """Checks if x1 and x2 have the same items the same number of times"""
    t1 = table(x1)
    t2 = table(x2)
    if t1.keys() != t2.keys():
        return False
    for k1, v1 in t1.items():
        if v1 != t2[k1]:
            return False
    return True


def table_index_map(x: Iterable) -> dict:
    """Like table, but maps to the indices with those elements"""
    unique = set(x)
    res: dict[Any, list[int]] = {}
    for u in unique:
        res[u] = []
    idx = 0
    for item in x:
        res[item] += [idx]
        idx += 1
    return res
