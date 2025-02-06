import json

from cladecombiner import Taxon
from cladecombiner.utils import arbitrary_taxonomy_helper


def test_nomenclature():
    # We've stipulated the root should be 19A, but are intentionally breaking that in the tests
    with open("tests/toy_arbitrary_tree.json") as f:
        nomenclature, _ = arbitrary_taxonomy_helper(
            json.load(f), expected_root="MONTY", name="test"
        )

    assert nomenclature.is_valid_name("PYTHONS")
    assert not nomenclature.is_valid_name("24A")

    assert nomenclature.is_root("MONTY")

    # # Mostly just to check they don't error out
    # assert not nomenclature.is_ambiguous("MONTY")
    # assert not nomenclature.is_hybrid("MONTY")


def test_taxonomy():
    # Mostly just a check that we get the tree right, tests of PTS handle the downstream behavior
    with open("tests/toy_arbitrary_tree.json") as f:
        _, taxonomy_scheme = arbitrary_taxonomy_helper(
            json.load(f), expected_root="MONTY", name="test"
        )

    tip_str = {"PATHG", "BRIAN", "CARNIVAL", "CIRCUS"}
    internal_str = {"MONTY", "PYTHONS", "LIFE", "OF", "FLYING"}

    assert tip_str == set(
        node.label for node in taxonomy_scheme.tree.leaf_node_iter()
    )
    assert internal_str == set(
        node.label
        for node in taxonomy_scheme.tree.preorder_internal_node_iter()
    )
    assert taxonomy_scheme.is_valid_taxon(Taxon("PATHG", True))
    assert not taxonomy_scheme.is_valid_taxon(Taxon("PATHG", False))
    assert taxonomy_scheme.parents(Taxon("FLYING", False)) == Taxon(
        "PYTHONS", False
    )
