import json

import pytest

from cladecombiner.tree_utils import (
    add_paraphyletic_tips,
    prune_nonancestral,
    tree_from_edge_dict,
)

r"""
The json file defines the following tree
  /-----PATHG
  |
MONTY      /-----LIFE-----OF-----BRIAN
  |        |
  \-----PYTHONS    /-----CARNIVAL
           |       |
           \-----FLYING
                   |
                   \-----CIRCUS
"""


@pytest.fixture
def tree():
    with open("tests/toy_arbitrary_tree.json") as f:
        tree = tree_from_edge_dict(json.load(f))
    return tree


def test_tree_correct(tree):
    expected_tips = {"PATHG", "BRIAN", "CARNIVAL", "CIRCUS"}
    tips = {leaf.label for leaf in tree.leaf_node_iter()}
    assert expected_tips == tips

    expected_internal = {"MONTY", "PYTHONS", "FLYING", "LIFE", "OF"}
    internal = {leaf.label for leaf in tree.leaf_node_iter()}
    assert not expected_internal == internal

    expected_children = {
        "MONTY": {"PYTHONS", "PATHG"},
        "PYTHONS": {"FLYING", "LIFE"},
        "FLYING": {"CARNIVAL", "CIRCUS"},
        "LIFE": {"OF"},
        "OF": {"BRIAN"},
    }

    for node in tree.postorder_internal_node_iter():
        observed_children = set(
            child.label for child in node.child_node_iter()
        )
        assert expected_children[node.label] == observed_children


def test_prune(tree):
    tree = add_paraphyletic_tips(tree, ["MONTY", "OF"])
    target_tips = ["MONTY", "OF", "CIRCUS"]
    tree = prune_nonancestral(tree, target_tips)
    print(
        tree.as_ascii_plot(plot_metric="level", show_internal_node_labels=True)
    )
    expected_tips = set(target_tips)
    tips = {leaf.label for leaf in tree.leaf_node_iter()}
    assert expected_tips == tips

    expected_internal = {"MONTY", "PYTHONS", "FLYING", "LIFE"}
    internal = {leaf.label for leaf in tree.leaf_node_iter()}
    assert not expected_internal == internal

    expected_children = {
        "MONTY": {"PYTHONS", "MONTY"},
        "PYTHONS": {"FLYING", "LIFE"},
        "FLYING": {"CIRCUS"},
        "LIFE": {"OF"},
    }

    for node in tree.postorder_internal_node_iter():
        observed_children = set(
            child.label for child in node.child_node_iter()
        )
        assert expected_children[node.label] == observed_children
