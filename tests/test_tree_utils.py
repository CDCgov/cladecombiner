from cladecombiner.tree_utils import (
    add_paraphyletic_tips,
    fully_labeled_trees_same,
    prune_nonancestral,
    tree_from_edge_table_string,
)


def test_tree_correct(arbitrary_taxonomy_tree):
    expected_tips = {"PATHG", "BRIAN", "CARNIVAL", "CIRCUS"}
    tips = {leaf.label for leaf in arbitrary_taxonomy_tree.leaf_node_iter()}
    assert expected_tips == tips

    expected_internal = {"MONTY", "PYTHONS", "FLYING", "LIFE", "OF"}
    internal = {
        leaf.label for leaf in arbitrary_taxonomy_tree.leaf_node_iter()
    }
    assert not expected_internal == internal

    expected_children = {
        "MONTY": {"PYTHONS", "PATHG"},
        "PYTHONS": {"FLYING", "LIFE"},
        "FLYING": {"CARNIVAL", "CIRCUS"},
        "LIFE": {"OF"},
        "OF": {"BRIAN"},
    }

    for node in arbitrary_taxonomy_tree.postorder_internal_node_iter():
        observed_children = set(
            child.label for child in node.child_node_iter()
        )
        assert expected_children[node.label] == observed_children


def test_prune(arbitrary_taxonomy_tree):
    tree = add_paraphyletic_tips(arbitrary_taxonomy_tree, ["MONTY", "OF"])
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


def test_treefile_parsing(arbitrary_taxonomy_tree):
    edge_table = """PYTHONS\tMONTY
FLYING\tPYTHONS
CIRCUS\tFLYING
LIFE\tPYTHONS
OF\tLIFE
BRIAN\tOF
CARNIVAL\tFLYING
PATHG\tMONTY"""
    parsed_tree = tree_from_edge_table_string(edge_table, "\t", 1, 0)
    assert fully_labeled_trees_same(parsed_tree, arbitrary_taxonomy_tree)
