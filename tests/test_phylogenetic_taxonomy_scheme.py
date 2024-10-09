import copy

import dendropy
import pytest

from cladecombiner import Taxon
from cladecombiner.taxon_utils import sort_taxa
from cladecombiner.tree_utils import fully_labeled_trees_same
from cladecombiner.utils import table_equal


def test_ancestors(pango_phylo_taxonomy):
    expected = [
        Taxon("XA.1.2", is_tip=False),
        Taxon("XA.1", is_tip=False),
        Taxon("XA", is_tip=False),
        Taxon("B.2", is_tip=False),
        Taxon("A.2.1.2", is_tip=False),
        Taxon("A.2.1", is_tip=False),
        Taxon("A.2", is_tip=False),
        Taxon("A", is_tip=False),
        Taxon("", is_tip=False),
    ]
    observed = pango_phylo_taxonomy.ancestors(Taxon("XA.1.2.3", is_tip=True))
    for o, e in zip(observed, expected):
        assert o == e


def test_children(pango_phylo_taxonomy):
    xa_children_exp = [
        Taxon("XA", is_tip=True),
        Taxon("XA.1", is_tip=False),
        Taxon("XA.12", is_tip=False),
    ]
    xa_children_obs = pango_phylo_taxonomy.children(Taxon("XA", is_tip=False))
    assert table_equal(xa_children_obs, xa_children_exp)

    assert len(pango_phylo_taxonomy.children(Taxon("A", is_tip=True))) == 0


def test_contains(pango_phylo_taxonomy):
    assert not pango_phylo_taxonomy.contains(
        focal=Taxon("A.1.1.2", is_tip=True),
        target=Taxon("A.1.1", is_tip=False),
    )

    assert pango_phylo_taxonomy.contains(
        focal=Taxon("A", is_tip=False), target=Taxon("A.1.1", is_tip=False)
    )

    assert pango_phylo_taxonomy.contains(
        focal=Taxon("A", is_tip=False), target=Taxon("XA", is_tip=False)
    )

    assert pango_phylo_taxonomy.contains(
        focal=Taxon("A.2.1", is_tip=False),
        target=Taxon("XA.12.3", is_tip=True),
    )


def test_descendants(pango_phylo_taxonomy):
    xa_desc_exp = [
        Taxon("XA", is_tip=True),
        Taxon("XA.1", is_tip=False),
        Taxon("XA.12", is_tip=False),
        Taxon("XA.12.3", is_tip=True),
        Taxon("XA.1.2", is_tip=False),
        Taxon("XA.1.2.3", is_tip=True),
    ]
    xa_desc_obs = pango_phylo_taxonomy.descendants(
        Taxon("XA", is_tip=False), tip_only=False
    )
    assert table_equal(xa_desc_obs, xa_desc_exp)

    xa_desc_exp = [
        Taxon("XA", is_tip=True),
        Taxon("XA.12.3", is_tip=True),
        Taxon("XA.1.2.3", is_tip=True),
    ]
    xa_desc_obs = pango_phylo_taxonomy.descendants(
        Taxon("XA", is_tip=False), tip_only=True
    )
    assert table_equal(xa_desc_obs, xa_desc_exp)


def test_mrca(pango_phylo_taxonomy):
    assert pango_phylo_taxonomy.mrca(
        [Taxon(r"A.1.1.2", is_tip=True), Taxon(r"A.1.1.3", is_tip=True)]
    ) == Taxon(r"A.1.1", is_tip=False)

    assert pango_phylo_taxonomy.mrca(
        [Taxon(r"A.1.1.2", is_tip=True), Taxon(r"A.2", is_tip=False)]
    ) == Taxon(r"A", is_tip=False)

    assert pango_phylo_taxonomy.mrca(
        [Taxon(r"XA.1.2", is_tip=False), Taxon(r"XA.12", is_tip=False)]
    ) == Taxon(r"XA", is_tip=False)

    assert pango_phylo_taxonomy.mrca(
        [Taxon(r"XA.1.2", is_tip=False), Taxon(r"B.2.4.8", is_tip=True)]
    ) == Taxon(r"B.2", is_tip=False)

    assert pango_phylo_taxonomy.mrca(
        [Taxon(r"A.1.1.2", is_tip=True), Taxon(r"B.2.4.8", is_tip=True)]
    ) == Taxon(r"A", is_tip=False)


def test_parents(pango_phylo_taxonomy):
    with pytest.raises(Exception):
        pango_phylo_taxonomy.parents(Taxon("A.1.1.2", is_tip=False))

    with pytest.raises(Exception):
        pango_phylo_taxonomy.parents(Taxon("B.2", is_tip=True))

    assert pango_phylo_taxonomy.parents(
        Taxon("A.1.1.2", is_tip=True)
    ) == Taxon("A.1.1", is_tip=False)

    assert pango_phylo_taxonomy.parents(Taxon("A.2", is_tip=False)) == Taxon(
        "A", is_tip=False
    )

    assert pango_phylo_taxonomy.parents(Taxon("XA.12", is_tip=False)) == Taxon(
        "XA", is_tip=False
    )

    assert pango_phylo_taxonomy.parents(Taxon("XA", is_tip=False)) == Taxon(
        "B.2", is_tip=False
    )


def test_prune(pango_phylo_taxonomy):
    copy_pango = copy.deepcopy(pango_phylo_taxonomy)
    copy_pango.prune_subtree(Taxon("A.2.1", False))

    a = dendropy.Node(label="A")
    a_1 = dendropy.Node(label="A.1")
    a_1_1 = dendropy.Node(label="A.1.1")
    a_1_1_2 = dendropy.Node(label="A.1.1.2")
    a_1_1_3 = dendropy.Node(label="A.1.1.3")

    a_2 = dendropy.Node(label="A.2")

    a_2_2 = dendropy.Node(label="A.2.2")
    a_2_2_3 = dendropy.Node(label="A.2.2.3")

    a.add_child(a_1)
    a.add_child(a_2)

    a_1.add_child(a_1_1)
    a_1_1.add_child(a_1_1_2)
    a_1_1.add_child(a_1_1_3)

    a_2.add_child(a_2_2)
    a_2_2.add_child(a_2_2_3)

    a_other = dendropy.Node(label="A")
    a.add_child(a_other)

    expected_tree = dendropy.Tree()
    # Always true, makes pylance happy about node access
    assert isinstance(expected_tree.seed_node, dendropy.Node)
    expected_tree.seed_node.label = ""
    expected_tree.seed_node.add_child(a)

    # print(
    #     copy_pango.tree.as_ascii_plot(
    #         plot_metric="level", show_internal_node_labels=True
    #     )
    # )

    assert fully_labeled_trees_same(copy_pango.tree, expected_tree)


def test_root(pango_phylo_taxonomy):
    assert pango_phylo_taxonomy.is_root(Taxon("", is_tip=False))


def test_valid(pango_phylo_taxonomy):
    expect_valid = [
        Taxon("XA.1.2", is_tip=False),
        Taxon("XA.1", is_tip=False),
        Taxon("XA", is_tip=False),
        Taxon("B.2", is_tip=False),
        Taxon("A.2.1.2", is_tip=False),
        Taxon("A.2.1", is_tip=False),
        Taxon("A.2", is_tip=False),
        Taxon("A", is_tip=False),
        Taxon("", is_tip=False),
    ]

    for taxon in expect_valid:
        assert pango_phylo_taxonomy.is_valid_taxon(taxon)

    expect_invalid = [
        Taxon("XA.1.2", is_tip=True),
        Taxon("XA.12.3", is_tip=False),
        Taxon("A.2", is_tip=True),
        Taxon("B.2.42", is_tip=True),
        Taxon("A.47", is_tip=False),
    ]
    for taxon in expect_invalid:
        assert not pango_phylo_taxonomy.is_valid_taxon(taxon)


def test_sort(pango_phylo_taxonomy):
    taxa = [
        Taxon("XA", is_tip=False),  # 0
        Taxon("A.1", is_tip=False),  # 1
        Taxon("A.2.2", is_tip=False),  # 2
        Taxon("A.2", is_tip=False),  # 3
        Taxon("A", is_tip=True),  # 4
        Taxon("A", is_tip=False),  # 5
        Taxon("", is_tip=False),  # 6
    ]

    staxa = list(sort_taxa(taxa, pango_phylo_taxonomy))

    # General correctness
    assert staxa.index(taxa[6]) == 6

    assert staxa.index(taxa[5]) > staxa.index(taxa[4])
    assert staxa.index(taxa[5]) > staxa.index(taxa[3])
    assert staxa.index(taxa[5]) > staxa.index(taxa[2])
    assert staxa.index(taxa[5]) > staxa.index(taxa[1])
    assert staxa.index(taxa[5]) > staxa.index(taxa[0])

    assert staxa.index(taxa[3]) > staxa.index(taxa[2])

    # Specific order we may get used to
