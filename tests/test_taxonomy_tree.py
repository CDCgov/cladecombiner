import dendropy

from cladecombiner.tree_utils import fully_labeled_trees_same


def test_tree_making(tax_tree):
    a = dendropy.Node(label="A")
    a_1 = dendropy.Node(label="A.1")
    a_1_1 = dendropy.Node(label="A.1.1")
    a_1_1_2 = dendropy.Node(label="A.1.1.2")
    a_1_1_3 = dendropy.Node(label="A.1.1.3")

    a_2 = dendropy.Node(label="A.2")
    a_2_1 = dendropy.Node(label="A.2.1")
    a_2_1_2 = dendropy.Node(label="A.2.1.2")

    b_1 = dendropy.Node(label="B.1")
    b_1_1 = dendropy.Node(label="B.1.1")
    b_1_1_3 = dendropy.Node(label="B.1.1.3")

    b_2 = dendropy.Node(label="B.2")
    b_2_4 = dendropy.Node(label="B.2.4")
    b_2_4_8 = dendropy.Node(label="B.2.4.8")

    a_2_2 = dendropy.Node(label="A.2.2")
    a_2_2_3 = dendropy.Node(label="A.2.2.3")

    xa = dendropy.Node(label="XA")
    xa_12 = dendropy.Node(label="XA.12")
    xa_12_3 = dendropy.Node(label="XA.12.3")
    xa_1 = dendropy.Node(label="XA.1")
    xa_1_2 = dendropy.Node(label="XA.1.2")
    xa_1_2_3 = dendropy.Node(label="XA.1.2.3")

    expected_tree = dendropy.Tree()
    # Always true, makes pylance happy about node access
    assert isinstance(expected_tree.seed_node, dendropy.Node)
    expected_tree.seed_node.label = ""
    expected_tree.seed_node.add_child(a)

    a.add_child(a_1)
    a.add_child(a_2)

    a_1.add_child(a_1_1)
    a_1_1.add_child(a_1_1_2)
    a_1_1.add_child(a_1_1_3)

    a_2.add_child(a_2_1)
    a_2_1.add_child(a_2_1_2)

    a_2_1_2.add_child(b_1)
    b_1.add_child(b_1_1)
    b_1_1.add_child(b_1_1_3)

    a_2_1_2.add_child(b_2)
    b_2.add_child(b_2_4)
    b_2_4.add_child(b_2_4_8)

    b_2.add_child(xa)
    xa.add_child(xa_12)
    xa_12.add_child(xa_12_3)

    xa.add_child(xa_1)
    xa_1.add_child(xa_1_2)
    xa_1_2.add_child(xa_1_2_3)

    a_2.add_child(a_2_2)
    a_2_2.add_child(a_2_2_3)

    # print(
    #     expected_tree.as_ascii_plot(
    #         plot_metric="level", show_internal_node_labels=True
    #     )
    # )
    assert fully_labeled_trees_same(expected_tree, tax_tree[0])

    # Add tips for ancestral taxa which we've observed
    a_other = dendropy.Node(label="A")
    xa_other = dendropy.Node(label="XA")
    a.add_child(a_other)
    xa.add_child(xa_other)

    # print(
    #     expected_tree.as_ascii_plot(
    #         plot_metric="level", show_internal_node_labels=True
    #     )
    # )
    assert fully_labeled_trees_same(expected_tree, tax_tree[1])
