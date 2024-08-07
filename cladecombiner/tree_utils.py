import copy
from collections.abc import Sequence

import dendropy


def add_paraphyletic_tips(
    phy: dendropy.Tree, tips: Sequence[str]
) -> dendropy.Tree:
    """
    Disambiguates ancestral versus tip taxa by adding tips explicitly.

    Assumes all nodes have labels.

    In nomenclatures for evolving pathogens, naming a new taxon will make
    a previously-named taxon paraphyletic. There can then be ambiguity with
    respect to whether that previous taxon name is being used to refer to
    the monophyletic group comprising this taxon and all its descendants,
    or the non-monophyletic group of the previous taxon except its newly
    named descendant.

    This function adds a tip to the phylogeny to represent the
    non-monophyletic group which has been split by subsequently-named
    taxa.

    For example, the SARS-CoV-2 Pango taxon JN.1 could mean the higher
    taxon JN.1 (which includes many more specifically-named taxa, such as
    JN.1.11.1 (KP) and JN.1.30.1 (KU)), or JN.1 as something we can
    observe as a label for sampled sequences. The latter of these means a
    non-more-specifically-named JN.1 lineage, some part of the tree of JN.1
    which has not been named more specifically. This also occurs with
    NextStrain clades, for example the SARS-CoV-2 clade 23I was made
    paraphyletic with respect to 24A, which was in turn made paraphyletic
    by 24B. So 23I can mean an ancestral taxon, comprising all lineages in
    any of these clades, or a non-more-specifically named part of the 23I
    tree, which we could see in a sample at the same time as we see 24A.

    Parameters
    ---------
    phy : dendropy.Tree with a label for all nodes
        The tree to which we will add the tips.
    tips : Sequence[str]
        The names of taxa that should exist as both ancestral and tip taxa.

    Returns
    -------
    dendropy.Tree
        The tree with all added tips.
    """
    tree = copy.deepcopy(phy)
    to_add = []
    for node in tree.preorder_node_iter():
        if node.is_internal():
            if node.label in tips:
                tip = dendropy.Node(label=node.label)
                to_add.append(
                    (
                        node,
                        tip,
                    )
                )
    for nt in to_add:
        nt[0].add_child(nt[1])

    return tree


def fully_labeled_subtrees_same(
    node1: dendropy.Node, node2: dendropy.Node
) -> bool:
    """
    Are two subtrees with every node labeled topologically equivalent?

    Used by fully_labeled_trees_same().

    Recursive function, calls itself until either a difference is seen or all
    tips in the subtree in both tree 1 and tree 2 are seen.

    Parameters
    ---------
    node1 : dendropy.Node
        Node defining the subtree in tree 1.
    node2 : dendropy.Node
        Node defining the subtree in tree 2.

    Returns
    -------
    bool
        True if the subtrees are the same.
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

    Standard topological identity means that two trees portray the same
    evolutionary relationships between the tips. This function assumes that
    every internal node is labeled and checks the relationships between all
    nodes.

    Calls fully_labeled_subtrees_same() to recursively evaluate subtrees.

    Parameters
    ---------
    tree1 : dendropy.Tree
        One tree to compare.
    tree2 : dendropy.Tree
        The other tree to compare.

    Returns
    -------
    bool
        True if the trees are the same.
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
