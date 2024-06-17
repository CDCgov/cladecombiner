from abc import ABC, abstractmethod
from collections.abc import Collection, Iterable, Sequence

import dendropy

from .taxon import Taxon


class TaxonomyScheme(ABC):
    """
    Abstract class for most general casting of Taxonomy

    Allows hybridization-induced multiple ancestry.
    """

    @abstractmethod
    def ancestors(self, taxon: Taxon) -> Collection[Taxon]:
        """All taxa which are between this taxon and the root (including the root)."""
        anc = set()
        queue = list(self.parents(taxon))
        while queue:
            tax = queue.pop(0)
            anc.add(tax)
            queue = [*queue, list(self.parents(tax))]
        return anc

    @abstractmethod
    def children(self, taxon: Taxon) -> Collection[Taxon]:
        """
        All taxa which are direct children of this taxon.

        Returns empty container if this taxon has no children.
        """
        pass

    @abstractmethod
    def descendants(self, taxon: Taxon) -> Collection[Taxon]:
        """All taxa which are contained by this taxon."""
        desc = set()
        queue = list(self.children(taxon))
        while queue:
            tax = queue.pop(0)
            desc.add(tax)
            queue = [*queue, list(self.children(tax))]
        return desc

    @abstractmethod
    def is_root(self, taxon: Taxon) -> bool:
        """Is this the largest taxon that contains all others?"""
        pass

    @abstractmethod
    def is_valid_taxon(self, taxon: Taxon) -> bool:
        """Does the scheme recognize this Taxon?"""
        pass

    @abstractmethod
    def parents(self, taxon: Taxon) -> Collection[Taxon]:
        """
        All parent taxa of taxon, e.g. ancestors exactly one level above this taxon.

        Hybridization allows a taxon to have multiple parent taxa.

        Returns empty container if this taxon is the root (which has no parents).
        """
        pass


class TreelikeTaxonomyScheme(TaxonomyScheme):
    """
    Abstract class for hybrid-free Taxonomy
    """

    def ancestors(self, taxon: Taxon) -> Sequence[Taxon]:
        """Postorder sequence of taxa between this taxon and the root (including the root)."""
        taxon = self.parents(taxon)
        anc = [taxon]
        while not self.is_root(taxon):
            taxon = self.parents(taxon)
            anc.append(taxon)
        return anc

    @abstractmethod
    def contains(self, focal: Taxon, target: Taxon) -> bool:
        """Does the focal taxon contain the target taxon?"""
        pass

    @abstractmethod
    def mrca(self, taxa: Iterable[Taxon]) -> Taxon:
        """Find the MRCA of a set of taxa"""
        pass

    @abstractmethod
    def parents(self, taxon: Taxon) -> Taxon:
        """A taxon has only one parent if the scheme is treelike"""
        pass


class PhylogeneticTaxonomyScheme(TreelikeTaxonomyScheme):
    """
    A TaxonomyScheme powered by a phylogeny.
    """

    def __init__(self, tree: dendropy.Tree):
        """
        PhylogeneticTaxonomyScheme constructor

        Required arguments:
        tree: the phylogeny to be used, internal nodes must be labeled
        """

        for node in tree.preorder_node_iter():
            if node.label is None:
                raise RuntimeError(
                    "TaxonomyTree constructor requires all nodes have labels."
                )
        self.tree = tree
        "The tree describing the relationships between taxa"
        self.node_to_taxon: dict[dendropy.Node, Taxon] = {}
        "The taxon represented by each node, for ease of access"
        self.taxon_to_node: dict[Taxon, dendropy.Node] = {}
        "The node representing each taxon, for ease of access"

        self.map_from_tree()

    #######################
    # Setup and utilities #
    #######################
    def map_from_tree(self) -> None:
        """
        Make Node<->Taxon maps

        By using these maps, we can avoid searching the tree repeatedly.
        """
        self.node_to_taxon = {}
        self.taxon_to_node = {}
        for node in self.tree.preorder_node_iter():
            self.node_to_taxon[node] = Taxon(node.label)
            self.taxon_to_node[self.node_to_taxon[node]] = node

    def node_path_to_root(self, taxon: Taxon) -> Sequence[dendropy.Node]:
        """
        Get all nodes between given taxon and the root (inclusive of the root)
        """
        path = []
        node = self.taxon_to_node[taxon]
        while node is not self.tree.seed_node:
            path.append(node)
            node = node.parent_node
        return path

    # @TODO: this class should have more of its own node management, which outsources to dendropy and cleans maps up
    def prune_subtree(self, taxon: Taxon) -> None:
        """
        Remove subtree corresponding to this taxon and clean up maps
        """
        # print(self.tree.as_ascii_plot(plot_metric="level"))
        node = self.taxon_to_node[taxon]
        node.parent_node.remove_child(node)
        self.map_from_tree()
        # print(self.tree.as_ascii_plot(plot_metric="level"))

    def root(self) -> dendropy.Node:
        if not isinstance(self.tree.seed_node, dendropy.Node):
            raise RuntimeError("Malformed tree has no seed_node")
        else:
            return self.tree.seed_node

    ########################
    # Superclass overrides #
    ########################
    def ancestors(self, taxon: Taxon) -> Sequence[Taxon]:
        nodes = self.node_path_to_root(taxon)
        return [node.label for node in nodes]

    def contains(self, focal: Taxon, target: Taxon) -> bool:
        if focal not in self.taxon_to_node:
            return False
        node_x = self.taxon_to_node[focal]

        if target not in self.taxon_to_node:
            return False
        node = self.taxon_to_node[target]

        while node is not self.tree.seed_node:
            if node is node_x:
                return True
            node = node.parent_node

        return False

    def descendants(self, taxon: Taxon, tip_only: bool) -> Collection[Taxon]:
        desc = []
        queue = [self.taxon_to_node[taxon]]
        while queue:
            node = queue.pop(0)
            has_kids = False
            for child in node.child_node_iter():
                has_kids = True
                queue.append(child)
            if (not tip_only) or (not has_kids):
                desc.append(self.node_to_taxon[node])
        if not tip_only:
            desc = desc[1:]
        return desc

    def is_root(self, taxon: Taxon) -> bool:
        return taxon == self.node_to_taxon[self.root()]

    def is_valid_taxon(self, taxon: Taxon) -> bool:
        return taxon in self.taxon_to_node.keys()

    def mrca(self, taxa: Sequence[Taxon]) -> Taxon:
        paths = [self.node_path_to_root(taxon) for taxon in taxa]
        # paths are tip->root, we need root->tip
        for path in paths:
            list(path).reverse()
        max_idx = max([len(path) for path in paths])
        idx = 0
        while idx < max_idx:
            tax_at_lvl = set([path[idx] for path in paths])
            if len(tax_at_lvl) == 1:
                idx += 1
            else:
                idx -= 1
                break
        if idx < 0:
            raise RuntimeError("Provided taxa do not have MRCA in the tree.")
        return self.node_to_taxon[paths[0][idx]]

    def parents(self, taxon: Taxon):
        node_x = self.taxon_to_node[taxon]
        if node_x is self.tree.seed_node:
            raise RuntimeError("The root has no containing group.")
        return self.node_to_taxon[node_x.parent_node]
