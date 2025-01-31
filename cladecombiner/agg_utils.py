from typing import Iterable

import dendropy

from .nomenclature import NomenclatureVersioner


def _get_versioned_tips(
    node: dendropy.Node,
    versioner: NomenclatureVersioner,
    clades: list[tuple[str, bool]],
):
    children_as_of = [
        child
        for child in node.child_node_iter()
        if child.label != node.label and versioner(child.label)
    ]

    if len(children_as_of) > 0:
        for child in children_as_of:
            _get_versioned_tips(child, versioner, clades)
    else:
        clades.append(
            (
                node.label,
                node.is_leaf(),
            )
        )


def get_versioned_tip_taxa(
    tree: dendropy.Tree, versioner: NomenclatureVersioner
) -> Iterable[tuple[str, bool]]:
    tips = []
    root = tree.seed_node
    assert root is not None
    assert versioner(root.label)
    _get_versioned_tips(root, versioner, tips)
    return tips
