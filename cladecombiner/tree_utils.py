import csv
from collections.abc import Sequence
from typing import Iterable

import dendropy


def _add_edict_children(
    node: dendropy.Node, name: str, parent_child: dict
) -> None:
    """
    Recursively add to tree from a parent-child edge dictionary.

    Parameters
    ---------
    node : dendropy.Node
        Tip of subtree being added to.
    name : str
        name of this node
    parent_child : dict
        Dict of all children of all nodes that will eventually be in the tree.
    """
    if name in parent_child:
        for child_name in parent_child[name]:
            child = node.new_child()
            child.label = child_name
            _add_edict_children(child, child_name, parent_child)


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
    tree = phy.clone(2)
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


def edge_dict_to_tree(
    child_parent: dict[str, str],
) -> dendropy.Tree:
    """
    Turns a {child : parent} dictionary of taxon names into a dendropy Tree.

    Parameters
    ---------
    child_parent : Optional[dict[str, list[str]]]
        Dict giving parent taxon for all non-root taxa.

    Returns
    -------
    PhylogeneticTaxonomyScheme
        A PhylogeneticTaxonomyScheme using the tree specified by the given dict.
    """
    children = set(child_parent.keys())
    parents = set(child_parent.values())

    roots = parents.difference(children)
    assert (
        len(roots) == 1
    ), f"There should be one root, not {len(roots)}. Found {roots}."
    str_root = roots.pop()

    parent_child = {}
    for k, v in child_parent.items():
        if v in parent_child:
            parent_child[v].append(k)
        else:
            parent_child[v] = [k]

    str_taxa = list(children)
    str_taxa.append(str_root)
    taxon_namespace = dendropy.TaxonNamespace(str_taxa)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)

    assert tree.seed_node is not None
    tree.seed_node.label = str_root
    _add_edict_children(tree.seed_node, str_root, parent_child)

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


def tree_from_edge_table_string(
    edge_table: str,
    delimiter: str,
    parent_col: str | int,
    child_col: str | int,
):
    """ """
    if isinstance(parent_col, str) and isinstance(child_col, str):
        use_names = True
    elif isinstance(parent_col, int) and isinstance(child_col, int):
        use_names = False
    else:
        raise TypeError(
            "Must specify either indices or names for extracting columns."
        )

    if use_names:
        reader = csv.DictReader(edge_table.split("\n"), delimiter=delimiter)
    else:
        reader = csv.reader(edge_table.split("\n"), delimiter=delimiter)

    child_parent = {row[child_col]: row[parent_col] for row in reader}  # type: ignore #pylance can't track that we've already sanitized this

    return edge_dict_to_tree(child_parent)


def _find_unobserved_subtrees(
    node: dendropy.Node, target_tips: set[str], subtrees: list
):
    """
    Recursively find nodes defining subtrees with no tips in `target_tips`
    """
    tip_descendants = set(n.label for n in node.leaf_nodes())
    if not tip_descendants.isdisjoint(target_tips):
        for child in node.child_node_iter():
            # Don't remove paraphyletic tips we've already added
            if child.label != node.label:
                _find_unobserved_subtrees(child, target_tips, subtrees)
    else:
        subtrees.append(node)


def prune_nonancestral(phy: dendropy.Tree, tips: Iterable[str]):
    """
    Prune a tree to only contain portions ancestral to provided tip taxa.

    Tree is assumed to be node labeled as a PhylogeneticTaxonomyScheme.tree.

    Parameters
    ---------
    phy : dendropy.Tree with a label for all nodes
        The tree which we are to prune.
    tips : Iterable[str]
        The names of the only tips to be found in the desired pruned tree.

    Returns
    -------
    dendropy.Tree
        The tree pruned to only portions ancestral to the tips.
    """
    tree = phy.clone(2)
    assert all(node.label is not None for node in tree.postorder_node_iter())

    target_tips = set(tips)
    tree_tips = set([leaf.label for leaf in tree.leaf_node_iter()])
    assert target_tips.issubset(
        tree_tips
    ), f"Tree is missing target tips {target_tips.difference(tree_tips)}"

    unobserved = []
    assert isinstance(tree.seed_node, dendropy.Node)

    _find_unobserved_subtrees(tree.seed_node, target_tips, unobserved)
    for node in unobserved:
        tree.prune_subtree(node, suppress_unifurcations=False)
