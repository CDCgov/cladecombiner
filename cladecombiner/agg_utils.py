from typing import Iterable

import dendropy

from .nomenclature import NomenclatureVersioner


def _get_clades_as_of(
    node: dendropy.Node, versioner: NomenclatureVersioner, clades: list[str]
):
    if node.is_leaf() and versioner(node.label):
        clades.append(node.label)
    elif versioner(node.label):
        for child in node.child_node_iter():
            _get_clades_as_of(child, versioner, clades)


def get_clades_as_of(
    tree: dendropy.Tree, versioner: NomenclatureVersioner
) -> Iterable[str]:
    clades = []
    root = tree.seed_node
    assert root is not None
    _get_clades_as_of(root, versioner, clades)
    return clades
