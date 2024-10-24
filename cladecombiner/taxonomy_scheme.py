from abc import ABC, abstractmethod
from collections.abc import Collection, Iterable, Sequence

import dendropy

from .taxon import Taxon


class TaxonomyScheme(ABC):
    """
    Abstract class for most general casting of Taxonomy

    Allows hybridization-induced multiple ancestry.
    """

    def ancestors(self, taxon: Taxon) -> Collection[Taxon]:
        """
        All taxa which are between this taxon and the root (including the root).

        Parameters
        ----------
        taxon : Taxon
            The taxon whose ancestors we want.

        Returns
        -------
        Collection[Taxon]
            All unique taxa between this taxon and the root.
            Empty container if this taxon is the root.
        """
        if self.is_root(taxon):
            return set()

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

        Parameters
        ----------
        taxon : Taxon
            The taxon whose children we want.

        Returns
        -------
        Collection[Taxon]
            A collection of the taxa that are direct descendants of this taxon.
            Returns empty container if this taxon has no children (i.e., if
            this taxon is a tip taxon).
        """
        raise NotImplementedError()

    @abstractmethod
    def descendants(self, taxon: Taxon, tip_only: bool) -> Collection[Taxon]:
        """
        All taxa which are contained by this taxon.

        Parameters
        ----------
        taxon : Taxon
            The taxon whose descendants we want.
        tip_only : bool
            Do we want only tip descendants of this taxon?

        Returns
        -------
        Collection[Taxon]
            If tip_only == True, all tips that are descended from this taxon.
            Otherwise, a collection of the taxa that descend from this taxon.
            That is, its children, and its childrens' children, and so forth.
            Returns empty container if this taxon is a tip.
        """
        desc = set()
        queue = list(self.children(taxon))
        while queue:
            tax = queue.pop(0)
            desc.add(tax)
            queue = [*queue, list(self.children(tax))]
        return desc

    @abstractmethod
    def is_root(self, taxon: Taxon) -> bool:
        """
        Is this the largest taxon that contains all others?

        Parameters
        ----------
        taxon: Taxon
            The taxon to be checked.

        Returns
        -------
        bool
            True if this taxon is the root.
        """
        raise NotImplementedError()

    @abstractmethod
    def is_valid_taxon(self, taxon: Taxon) -> bool:
        """
        Does the scheme recognize this Taxon?

        Parameters
        ----------
        taxon: Taxon
            The taxon to be checked.

        Returns
        -------
        bool
            True if this taxon is valid.
        """
        raise NotImplementedError()

    @abstractmethod
    def parents(self, taxon: Taxon) -> Collection[Taxon]:
        """
        All parent taxa of taxon, e.g. ancestors exactly one level above this taxon.

        Hybridization allows a taxon to have multiple parent taxa.

        Parameters
        ----------
        taxon : Taxon
            The taxon whose parents we want.

        Returns
        -------
        Collection[Taxon]
            A collection of the taxa that are direct parents of this taxon.
            Returns empty container if this taxon is the root.
        """
        raise NotImplementedError()


class TreelikeTaxonomyScheme(TaxonomyScheme):
    """
    Abstract class for hybrid-free Taxonomy.

    Common taxonomic notions that are either ill-defined or require
    generalization in the face of hybridization are defined here, such as the
    MRCA of a set of taxa.
    """

    ########################
    # Superclass overrides #
    ########################

    def ancestors(self, taxon: Taxon) -> Sequence[Taxon]:
        """
        Postorder sequence of taxa between this taxon and the root (including
        the root).

        Parameters
        ----------
        taxon : Taxon
            The taxon whose ancestors we want.

        Returns
        -------
        Sequence[Taxon]
            All unique taxa between this taxon and the root, in that order,
            and including the root. Returns empty container if this taxon is
            the root.
        """
        anc = []
        parent = self.parents(taxon)
        while parent is not None:
            anc.append(parent)
            parent = self.parents(parent)
        return anc

    @abstractmethod
    def parents(self, taxon: Taxon) -> Taxon | None:
        """
        A taxon has only one parent if the scheme is treelike.

        Parameters
        ----------
        taxon : Taxon
            The taxon whose parents we want.

        Returns
        -------
        Taxon
            The taxon's parent, or None if this is the root.
        """
        raise NotImplementedError()

    #################
    # Class methods #
    #################

    @abstractmethod
    def contains(self, focal: Taxon, target: Taxon) -> bool:
        """
        Does the focal taxon contain the target taxon?

        That is, is target a descendant of focal?

        Parameters
        ----------
        focal : Taxon
            This taxon may or may not contain the target taxon.
        target : Taxon
            The taxon which may or may not be contained by the focal taxon.

        Returns
        -------
        bool
            True if focal contains target.
        """
        raise NotImplementedError()

    @abstractmethod
    def mrca(self, taxa: Iterable[Taxon]) -> Taxon:
        """
        Find the MRCA of a set of taxa

        The MRCA is the most recent common ancestor of a set of taxa. There
        are potentially many common ancestors of a particular group of taxa,
        but this is the one which contains the fewest other taxa possible.

        Parameters
        ----------
        taxa : Iterable[Taxon]
            The taxa for which we want the MRCA.

        Returns
        -------
        Taxon
            The MRCA.
        """
        raise NotImplementedError()


class PhylogeneticTaxonomyScheme(TreelikeTaxonomyScheme):
    """
    A TaxonomyScheme powered by a phylogeny.

    Errors are provoked when a PhylogeneticTaxonomyScheme is queried about taxa
    that are not in the phylogeny.

    Internally, a dendropy.Tree object is used to represent the taxonomic
    relationships.
    """

    def __init__(self, tree: dendropy.Tree):
        """
        PhylogeneticTaxonomyScheme constructor

        Parameters
        ----------
        tree : dendropy.Tree
            The phylogeny to be used, internal nodes must be labeled.
        """

        for node in tree.preorder_node_iter():
            if node.label is None:
                raise ValueError(
                    "TaxonomyTree constructor requires all nodes have labels."
                )
        self.tree = tree
        "The tree describing the relationships between taxa"
        self.node_to_taxon: dict[dendropy.Node, Taxon] = {}
        "The taxon represented by each node, for ease of access"
        self.taxon_to_node: dict[Taxon, dendropy.Node] = {}
        "The node representing each taxon, for ease of access"

        self.map_from_tree()

    ########################
    # Superclass overrides #
    #                      #
    # These change only    #
    # how methods work,    #
    # not what they return #
    ########################

    def ancestors(self, taxon: Taxon) -> Sequence[Taxon]:
        nodes = self.node_path_to_root(taxon)
        return [self.node_to_taxon[node] for node in nodes[1:]]

    def children(self, taxon: Taxon) -> Collection[Taxon]:
        child_nodes = self.taxon_to_node[taxon].child_nodes()
        children = []
        if child_nodes:
            for node in child_nodes:
                children.append(self.node_to_taxon[node])
        return children

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

        if node is node_x:
            return True

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
        max_idx = max([len(path) for path in paths])
        idx = 0
        while idx < max_idx:
            tax_at_lvl = set([path[-(1 + idx)] for path in paths])
            if len(tax_at_lvl) == 1:
                idx += 1
            else:
                idx -= 1
                break
        if idx < 0:
            raise RuntimeError("Provided taxa do not have MRCA in the tree.")
        return self.node_to_taxon[paths[0][-(1 + idx)]]

    def parents(self, taxon: Taxon) -> Taxon | None:
        node_x = self.taxon_to_node[taxon]
        if node_x is self.tree.seed_node:
            return None
        return self.node_to_taxon[node_x.parent_node]

    #################
    # Class methods #
    #################

    def map_from_tree(self) -> None:
        """
        Make Node<->Taxon maps

        By using these maps, we can avoid searching the tree repeatedly.

        Returns
        -------
        None
            Modifies self.node_to_taxon and self.taxon_to_node in-place.
        """

        self.node_to_taxon = {}
        self.taxon_to_node = {}
        for node in self.tree.preorder_node_iter():
            assert isinstance(node.is_leaf(), bool)  # Pylance paranoia
            taxon = Taxon(node.label, is_tip=node.is_leaf())
            self.node_to_taxon[node] = taxon
            self.taxon_to_node[taxon] = node

    def node_path_to_root(self, taxon: Taxon) -> Sequence[dendropy.Node]:
        """
        Get all nodes between given taxon and the root (inclusive of the root and this node)

        Parameters
        ----------
        taxon : Taxon
            The taxon for which we want the path to the root.

        Returns
        -------
        Sequence[dendropy.Node]
            Path of nodes in self.tree from this taxon (inclusive) to the root
            (inclusive).
        """
        path = []
        node = self.taxon_to_node[taxon]
        while node is not self.tree.seed_node:
            path.append(node)
            node = node.parent_node
        path.append(self.root())
        return path

    def prune_subtree(self, taxon: Taxon) -> None:
        """
        Remove subtree corresponding to this taxon and clean up maps

        Parameters
        ----------
        taxon : Taxon
            The taxon which is the base of the subtree to be removed.

        Returns
        -------
        None
            Edits self.tree in-place.
        """
        node = self.taxon_to_node[taxon]
        node.parent_node.remove_child(node)
        self.map_from_tree()

    def root(self) -> dendropy.Node:
        """
        Typing-safe function to access root, always returns a dendropy.Node.

        Returns
        -------
        dendropy.Node
            The root node of the phylogeny underlying this taxonomy scheme.
        """
        if not isinstance(self.tree.seed_node, dendropy.Node):
            raise RuntimeError("Malformed tree has no seed_node")
        else:
            return self.tree.seed_node
