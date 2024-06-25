import json
import string
import urllib.request
import warnings
from abc import ABC, abstractmethod
from collections.abc import Collection, Container, MutableSequence, Sequence
from sys import maxsize as integer_inf
from typing import Any, Callable, Optional
import validators

import dendropy

from .taxon import Taxon
from .tree_utils import add_paraphyletic_tips


class Nomenclature(ABC):
    """
    Abstract class for most general casting of Nomenclature

    Nomenclature concerns rules for naming taxa, and what names may imply about those taxa.
    """

    @staticmethod
    @abstractmethod
    def is_ambiguous(name: str) -> bool:
        """
        Does this name indicate an ambiguous taxon?

        An ambiguous taxon is known only to a higher level than to which resolution is possible.

        Returns
        -------
        bool
            True if this name indicates an ambiguous taxon.
        """
        pass

    @staticmethod
    @abstractmethod
    def is_hybrid(name: str) -> bool:
        """
        Does this name indicate a hybrid?

        Parameters
        ----------
        name : string specifying name of the taxon

        Returns
        -------
        bool
            True if this name indicates a hybrid taxon.
        """
        pass

    @staticmethod
    @abstractmethod
    def is_root(name: str) -> bool:
        """
        Does this string specify the root taxon (to which all other taxa belong)?

        Parameters
        ----------
        name : string specifying name of the taxon

        Returns
        -------
        bool
            True if this name indicates the root taxon.
        """
        pass

    @staticmethod
    @abstractmethod
    def is_valid_name(name: str) -> bool:
        """
        Is this name valid in the nomenclature scheme?

        Parameters
        ----------
        name : string specifying name of the taxon

        Returns
        -------
        bool
            True if this is a valid name under the nomenclature.
        """
        pass


class AlgorithmicNomenclature(Nomenclature):
    """
    Abstract class Nomenclature schemes which encode a taxon's history in its
    name in some form.

    The primary exemplar is the Pango nomenclature, which descends from this
    class via the more-general PangoLikeNomenclature.

    This class assumes that the history of a set of taxa can be decoded (in
    some way), the result being for each taxon a Sequence of taxa linking
    the root to it. A method is provided for constructing from these histories
    a tree suitable for use in PhylogeneticTaxonomyScheme.
    """

    @abstractmethod
    def history(self, taxon: str, stop_at_hybrid: bool = False) -> Sequence[str]:
        """
        Get the sequence of names of ancestors from the root to this taxon.

        Parameters
        ----------
        taxa : str
            taxon name
        stop_at_hybrid : boolean
            If True, the history starts at the most recent
            hybridization event in its ancestry. If False, we extract a linear
            history by taking the ancestry through the first indicated parent
            every time.

        Returns
        -------
        Sequence[str]
            list of taxon names
        """
        pass

    @staticmethod
    def add_history_to_tree(root: dendropy.Node, history: Sequence[str]) -> None:
        """
        Add a root-to-taxon history to an existing tree.

        Parameters
        ----------
        root : dendropy.Node
            root node
        history: Sequence[str]
            list of taxon names, from root down

        Returns
        -------
        None
            Modifies tree in-place.
        """
        # if no history, nothing to do
        if len(history) == 0:
            return None

        # root should be what we think it is
        assert root.label == history[0]

        # go down the history, following existing nodes where they exist,
        # creating new ones when needed
        node = root

        for taxon in history[1:]:
            # look at the children of the working node. do any of these match the taxon name?
            child_names = [x.label for x in node.children]
            taxon_exists = taxon in child_names

            if taxon_exists:
                # find which node it is, and set that as the new working taxon
                i = child_names.index(taxon)
                node = node.children[i]
            else:
                # create a new node, and set that as the working node
                new_node = dendropy.Node(label=taxon)
                node.add_child(new_node)
                node = new_node

            # double check that we're at the right place in the tree
            assert node.label == taxon

    @classmethod
    def taxonomy_tree(
        cls,
        taxa: Sequence[Taxon],
        insert_tips: bool,
        warn: bool = True,
    ) -> dendropy.Tree:
        """
        Makes a taxonomy tree for a set of taxa.

        A taxonomy tree is the core object of a PhylogeneticTaxonomyScheme,
        being a phylogenetic representation of the relationships between
        all taxa. It takes the form of a dendropy.Tree object where every
        node has a label.

        Parameters
        ----------
        taxa : Sequence[Taxon]
            We will build the tree of these taxa.
        insert_tips : boolean
            If True, where a Taxon in the provided taxa is an internal node,
            a tip is added to represent any paraphyletic observations of this
            taxon using add_paraphyletic_tips().
        warn : bool
            Should we warn the user if any taxa are dropped in the process
            of making the tree?

        Returns
        -------
        dendropy.Tree object with all nodes labeled
            The taxonomy tree is given by the phylogeny and all nodes are
            labeled with the taxon they represent. This tree may have nodes
            with only one descendant.
        """
        unique_names = list(set([taxon.name for taxon in taxa if taxon.tip]))
        if warn and (len(unique_names) < len(taxa)):
            warnings.warn("Removed non-unique and/or non-tip taxa to build tree.")

        if len(taxa) == 0:
            raise RuntimeError("Need at least one taxon")

        histories = [cls.history(name) for name in unique_names]
        # check that all histories have the same root
        if len(set(history[0] for history in histories)) > 1:
            raise RuntimeError("Cannot start tree, not all histories have same root")

        all_names = [taxon for history in histories for taxon in history]

        namespace = dendropy.TaxonNamespace(list(all_names))
        phy = dendropy.Tree(taxon_namespace=namespace)
        root = phy.seed_node
        if not isinstance(root, dendropy.Node):
            # Should never hit, required for type checking
            raise RuntimeError(
                "Cannot start tree because seed_node is not a dendropy.Node"
            )

        # set root label
        root.label = histories[0][0]

        # build out the tree from the histories
        for history in histories:
            cls.add_history_to_tree(root, history)

        # validate and clean tree node names
        for node in phy.preorder_node_iter():
            node.label = cls.clean_name(node.label)
            assert cls.is_valid_name(node.label)

        # add paraphyletic tips, if requested
        if insert_tips:
            phy = add_paraphyletic_tips(phy, unique_names)

        # check for duplicated names
        tip_names = [node.label for node in phy.leaf_node_iter()]
        int_names = [node.label for node in phy.preorder_internal_node_iter()]

        duplicate_tip_names = [x for x in tip_names if tip_names.count(x) > 1]
        duplicate_int_names = [x for x in int_names if int_names.count(x) > 1]

        if len(duplicate_int_names) > 0:
            raise RuntimeError(f"Duplicate tip names: {duplicate_tip_names}")

        if len(duplicate_int_names) > 0:
            raise RuntimeError(f"Duplicate internal node names: {duplicate_int_names}")

        return phy

    @classmethod
    def clean_name(cls, name):
        # by default, the name cleaning is to do nothing
        return name


class PangoLikeNomenclature(AlgorithmicNomenclature):
    def __init__(
        self,
        name: str,
        alias_map_path: str,
        alias_map_hybrid: Collection[type],
        charsets: Sequence[set],
        cumulative_alias: bool,
        max_sublevels: int,
        root: str,
        sep: str,
        special: Container,
        alias_map_inv: dict = {},
    ):
        """
        A Pango-like nomenclature is an AlgorithmicNomenclature with more specific
        assumptions about the encoding of history.

        Specifically, we assume that the name encodes the history in a string such
        that the name is a series of (sub)levels denoted by a consistent set of
        characters (say, digits) separated by a consistent separator (say, r".").
        The first portion of the name is assumed to be an alias, which is a set of
        different characters (say, upper case letters) which serve as shorthand
        for a longer series of levels. The alias is allowed to be cumulative (such
        as in RSV nomenclature) or not (such as in Pango nomenclature).

        An external file storing the alias shortcuts is required.

        Args:
            name: nomenclature system name
            alias_map_path (str): URL or filepath. The indicated place should
              be a json, which maps from the alias to a longer form of the name
              (but not necessarily the longest possible form)
            alias_map_hybrid (Collection[type]): Container type(s) used in alias map when hybrid ancestry is indicated
            charsets (Sequence[set]): what's allowed in alias names [0] and sublevel names [1]
            cumulative_alias (bool): Does the alias accumulate (like RSV system) or not (like Pango)
            max_sublevels (int): maximum number of sublevels before aliasing must be done
            root (str): name of root taxon
            sep (str): what separates the levels of the name
            special (Container): what aliases are allowed to appear alone
        """
        self.name = name
        self.alias_map_path = alias_map_path
        self.alias_map_hybrid = alias_map_hybrid
        self.charsets: charsets
        self.cumulative_alias = cumulative_alias
        self.max_sublevels = max_sublevels
        self.root = root
        self.sep = sep
        self.special = special

        # load the raw alias map
        if validators.url(self.alias_map_path):
            with urllib.request.urlopen(self.alias_map_path) as response:
                self.alias_map_raw = json.loads(response.read().decode("utf8"))
        else:
            with open(self.alias_map_path) as f:
                self.alias_map_raw = dict(json.load(f))

        # clean the map
        self.alias_map = self.sanitize_alias_map(self.alias_map_raw)
        # invert map, to go from shorter to longer
        self.alias_map_inv = self.invert_alias_map(self.alias_map)

    ##############################
    # Superclass implementations #
    ##############################

    def history(self, name: str, stop_at_hybrid: bool) -> Sequence[str]:
        """
        Get a path of ancestry from the root to this taxon.

        This is different than a long-form name because it allows us to pass
        through hybridization (recombination) events. In the face of
        recombination, when stop_at_hybrid == False, we follow the ancestry of
        the 5'-most portion of the genome.

        Parameters
        ----------
        name : str
            A taxon's name.
        stop_at_hybrid : bool
            Raise error if we hit a hybrid?

        Returns
        -------
        int
            Position of match in taxa, or -1 if no match
        """
        if stop_at_hybrid:
            raise NotImplementedError()

        assert self.alias_map is not None

        if not self.is_valid_name(name):
            raise ValueError(f"{name} is not a valid name")
        history = []
        self.extend_history(name, history, stop_at_hybrid)
        history.reverse()
        return history

    def is_valid_name(
        self,
        name: str,
        min_sublevels: int = 1,
        max_sublevels: Optional[int] = None,
    ) -> bool:
        parts = self.partition_name(name)
        # Check aliasing portion of name
        n_alias = len(parts[0])
        if n_alias < 1:
            return False
        if n_alias > 1 and not self.cumulative_alias:
            return False
        for a in parts[0]:
            if not set(a) < self.charsets[0]:
                return False
        # Check sublevels
        n_lvl = len(parts[1])
        if n_lvl < min_sublevels:
            return False
        if max_sublevels is None:
            if n_lvl > self.max_sublevels:
                return False
        elif n_lvl > max_sublevels:
            return False
        for lvl in parts[1]:
            if not set(lvl) < self.charsets[1]:
                return False
        return True

    ########################
    # Superclass overrides #
    ########################

    def taxonomy_tree(
        self,
        taxa: Sequence[Taxon],
        insert_tips: bool = True,
        warn: bool = True,
    ) -> dendropy.Tree:
        return super().taxonomy_tree(
            taxa=taxa,
            insert_tips=insert_tips,
            name_cleanup_fun=self.coax_name,
            warn=warn,
        )

    #################
    # Class methods #
    #################

    @classmethod
    def clean_name(cls, name: str) -> str:
        """
        Coax a potentially too-short or too-long name to proper format.

        For example, we might coax the SARS-CoV-2 Pango name from
        B.1.1.529.2.86.1.1.11.1.3 (which encodes the entire history but is too)
        long under the scheme to be proper, to KP.3. Alternately, we might coax
        the too-short KP to JN.1.11.1

        Parameters
        ----------
        name : str
            The name of the taxon.

        Returns
        -------
        str
            The name, without too many or too few sublevels.
        """
        if cls.is_root(name):
            return name
        else:
            return cls.shorter_name(cls.longer_name(name))

    @classmethod
    def equals_ignore_alias(cls, x: str, y: str) -> bool:
        """
        Are two names the same, accounting for aliasing?

        For example, the Pango SARS-CoV-2 names JN.1.11.1.3 and KP.3 both
        encode the history of the same taxon, KP3.

        Parameters
        ----------
        x : str
            A taxon's name.
        y : str
            A putatively equivalent name for the taxon

        Returns
        -------
        bool
            Are the names the same ignoring aliasing?
        """
        return x == cls.longer_name(y) or y == cls.longer_name(x)

    def extend_history(
        self, name: str, history: MutableSequence[str], stop_at_hybrid: bool
    ) -> None:
        """
        Recursively extend a path of ancestry from this taxon to the root.

        Parameters
        ----------
        taxon : str
            A taxon's name.
        history : MutableSequence[str]
            The history we are in the process of building
        stop_at_hybrid : bool
            Should we consider hybridization to start a new tree or not? If
            not, we break hybridization by following the first listed parent.

        Returns
        -------
        None
            Adds history to the history argument and then returns or calls
            itself if not done.
        """
        name = self.longer_name(name)
        comp = self.partition_name(name)
        if not comp[0]:
            raise ValueError(f"Invalid name: {repr(name)}")
        # Digest sublevels
        if comp[1]:
            for i in range(1, len(comp[1]) + 1)[::-1]:
                history.append(self.unpartition_name([comp[0], comp[1][:i]]))
        # Handle alias
        alias = self.join_name(comp[0])
        if self.is_root(alias):
            history.append(self.root)
        else:
            if self.is_special(alias):
                history.append(alias)
            if not self.is_hybrid(alias):
                self.extend_history(self.alias_map[alias], history, stop_at_hybrid)
            elif not stop_at_hybrid:
                self.extend_history(self.alias_map[alias][0], history, stop_at_hybrid)

    @classmethod
    def invert_alias_map(cls, x: dict) -> dict:
        """
        Inverts the shorter->longer alias map

        The inverted alias map is incapable of handling hybridization.

        Returns
        -------
        dict
            inverted (longer->shorter) map
        """
        rev_map = {}
        for k, v in x.items():
            if not isinstance(v, list):
                v = [v]
            for vi in v:
                # Don't add empty root alias
                if (not cls.is_root(vi)) and (not cls.is_hybrid(k)):
                    if vi in rev_map:
                        raise RuntimeError(
                            "Alias list cannot be inverted. Trying to add "
                            f"inverse alias for {vi}, which is an alias of"
                            f"{k}, but reversed map already has it"
                        )
                    rev_map[vi] = k

        return rev_map

    def is_alias_map_hybrid(self, alias_value: Any) -> bool:
        """
        Is this lineage a hybrid according to the alias map?

        Checks whether a value (rather than a key) from an alias map indicates
        a taxon has hybrid ancestry by checking if it is a container.

        Parameters
        ----------
        alias_value : the value (as opposed to the key) for some taxon in
            self.alias_map, i.e., self.alias_map[<some key>]

        Returns
        -------
        bool
            True if the alias map indicates this is a hybrid.
        """
        return any([isinstance(alias_value, t) for t in self.alias_map_hybrid])

    @abstractmethod
    def is_special(self, name: str) -> bool:
        """
        Is this a recognized special-purpose ancestor?

        Special-purpose ancestors are allowed to be used with 0 sublevels.

        Under the Pango nomenclature, direct root descendants and recombinants
        are special-purpose ancestors. Thus for Pango SARS-CoV-2, a special
        lineage is A, B, or any recombinant such as XBB (but not a descendant,
        such as XBB.1).

        Parameters
        ----------
        name : str
            A taxon's name.

        Returns
        -------
        bool
            True if this taxon is a special taxon.
        """
        pass

    def is_valid_alias(self, alias: str) -> bool:
        """
        Does this string specify a valid shortcut/alias for a taxon's history?

        A valid alias should contain only characters allowed in the aliasing
        portion of the name, possibly with separators if the alias is
        cumulative.

        Parameters
        ----------
        alias : str
            String to be checked for validity as alias.

        Returns
        -------
        bool
            True if this is a valid alias.
        """
        if self.cumulative_alias:
            return all([set(a) < self.charsets[0] for a in self.split_name(alias)])
        else:
            return set(alias) < self.charsets[0]

    def join_name(self, comp: Sequence[str]) -> str:
        """
        Join list of component levels into name.

        The inverse of self.split_name(name), such that
        self.join_name(self.split_name(name)) == name.

        Parameters
        ----------
        comp : Sequence[str]
            Components of a taxon's name.

        Returns
        -------
        str
            The name a a single string.
        """
        return self.sep.join(comp)

    def longer_name(self, name: str) -> str:
        """
        Get non-aliased form of an aliased name.

        A long-form name stops at the most recent hybridization event in a
        taxon's ancestry if there is such an event, otherwise at the special
        root descendent taxa.

        For example, the Pango SARS-CoV-2 taxon JN.1.11 would become
        B.1.1.529.2.86.1.1.11.

        Parameters
        ----------
        name : str
            A taxon's name.

        Returns
        -------
        str
            The taxon's name in the longest form of history.
        """
        if not self.alias_map:
            raise RuntimeError(
                "Cannot construct long form of name without an alias list."
            )
        if self.is_root(name):
            return name
        alias_levels = list(self.partition_name(name))
        next_alias = self.alias_map[alias_levels[0][-1]]
        while not self.is_root(next_alias) and (
            not self.is_hybrid(alias_levels[0][-1])
        ):
            parts = self.partition_name(next_alias)
            alias_levels[0] = parts[0]
            alias_levels[1] = [*parts[1], *alias_levels[1]]
            next_alias = self.alias_map[parts[0][-1]]
        return self.unpartition_name(alias_levels)

    def next_shorter_alias(self, name: str, depth: int) -> str:
        """
        Get the next shortest name available to a taxon.

        This removes one "layer" of self.max_sublevels from a name. For
        example, the Pango SARS-CoV-2 lineage B.1.1.529.2.86.1.1.11 would
        become BA.2.86.1.1.11 because BA is an alias for B.1.1.529.

        Parameters
        ----------
        name : str
            A expanded taxon name to be contracted
        depth : int
            How many levels of aliasing deep is this name? Starting at 1 for
            longest (fully de-aliased) name and increasing as the name gets
            shorter.

        Returns
        -------
        str
            The taxon's name with one fewer levels of aliasing.
        """

        parts = self.partition_name(name)
        n = self.max_sublevels * depth
        alias = None

        for k, v in self.alias_map_inv.items():
            kl = self.partition_name(k)
            if kl[1] == parts[1] or (
                self.partition_name(self.longer_name(k))[1][:n]
                == self.partition_name(self.longer_name(name))[1][:n]
            ):
                alias = v
                break

        if not alias:
            raise RuntimeError(f"Cannot find shorter alias for {repr(name)}")

        return alias

    def num_sublevels(self, name: str) -> int:
        """
        How many sublevels does this name contain?

        For a Pango SARS-CoV-2 example, the names XBB, XBB.1, XBB.1.5, and
        XBB.1.5.39 contain 0, 1, 2, and 3 sublevels respectively.

        Parameters
        ----------
        name : str
            The taxon's name.

        Returns
        -------
        int
            The number of sublevels the name contains.
        """
        return len(self.partition_name(name)[1])

    def partition_name(self, name: str) -> Sequence[Sequence[str]]:
        """
        Splits name into alias and sublevels, each as a sequence of components

        This function assumes that the name is ordered alias, sublevels, and does not check correctness.

        Parameters
        ----------
        name : str
            The taxon's name to be partitioned.

        Returns
        -------
        Sequence[Sequence[str]]
            First element is Sequence of components in the aliasing portion of
            the taxon's name, second element is Sequence of sublevels.
        """
        comp = self.split_name(name)
        if not self.cumulative_alias:
            alias = [] if len(comp) == 0 else [comp[0]]
            sublevels = [] if len(comp) < 2 else comp[1:]
            return [alias, sublevels]

        if not (set(comp[-1]) < self.charsets[1]):
            return [comp, []]
        if not (set(comp[0]) < self.charsets[0]):
            return [[], comp]
        n = 1
        while set(comp[-n]) < self.charsets[1]:
            n += 1

        alias = [comp[0]] if n == len(comp) else comp[: (len(comp) - n + 1)]
        sublevels = [comp[-1]] if n == 1 else comp[(len(comp) - n + 1) :]
        return [alias, sublevels]

    @classmethod
    def sanitize_alias_map(cls, x: dict) -> dict:
        """
        Drop ambiguity markers and check all names are valid.

        For the purposes of determining ancestry, an unknown sublineage is
            effectively just its ancestor, and we treat it as such.

        Returns
        -------
        dict
            new alias map
        """
        # make a copy of the input map, to not edit in-place
        x = dict(x)

        for k, v in x.items():
            # check for validity, ambiguity, and valid root mapping
            if not cls.is_valid_alias(k):
                raise RuntimeError(
                    f"Found invalid taxon {repr(k)} as key in alias list"
                )
            if cls.is_ambiguous(k):
                raise RuntimeError(
                    f"Found ambiguous taxon {repr(k)} as key in alias list"
                )
            if (
                not cls.is_alias_map_hybrid(v)
                and cls.is_root(v)
                and not cls.is_special(k)
            ):
                raise RuntimeError(
                    f"Found alias for root in taxon not listed as special: {repr(k)}"
                )

            # ensure list
            if not cls.is_alias_map_hybrid(v) and not cls.is_root(v):
                v = [v]

            for i in range(len(v)):
                if cls.is_ambiguous(v[i]):
                    v[i] = v[i][:-1]
                if not cls.is_valid_name(v[i], max_sublevels=integer_inf):
                    raise RuntimeError(
                        f"Found invalid taxon {repr(v[i])} as value in alias list "
                        f"(with key {repr(k)})"
                    )

        return x

    def shorter_name(self, name: str) -> str:
        """
        Get shortest form of a maximally-long name using aliases

        For example, the SARS-CoV-2 Pango name B.1.1.529.2.86.1.1.11.1.3 will
        be made into KP.3, and B.1.1.529.2.86.1.1.11.1 will be made into
        JN.1.11.1. Both of these are the shortest-possible valid forms of the
        names, having neither too many nor too few sublevels.

        Parameters
        ----------
        name : str
            The taxon's name to be shortened.

        Returns
        -------
        str
            Shortest valid form of the name for this taxon.
        """
        if not self.alias_map:
            raise RuntimeError("Cannot get shorter name without an alias list.")
        comp = list(self.partition_name(name))
        lvl = 1
        while len(comp[1]) > self.max_sublevels:
            alias = self.next_shorter_alias(self.unpartition_name(comp), lvl)
            comp[0] = [alias]
            comp[1] = comp[1][3:]
            lvl += 1
        return self.unpartition_name(comp)

    def split_name(self, name: str) -> Sequence[str]:
        """
        Split name into component levels

        The inverse of self.join_name(name), such that
        self.split_name(self.join_name(components)) == components.

        Parameters
        ----------
        name : str
            The name a a single string.

        Returns
        -------
        Sequence[str]
            Components of a taxon's name.
        """
        return name.split(self.sep)

    def unpartition_name(self, components: Sequence[Sequence[str]]) -> str:
        """
        Undoes partition_name

        Parameters
        ----------
        Sequence[Sequence[str]]
            First element is Sequence of components in the aliasing portion of
            the taxon's name, second element is Sequence of sublevels.

        Returns
        -------
        str
            The taxon's name as a single string.
        """
        return self.join_name([*components[0], *components[1]])


class PangoNomenclature(PangoLikeNomenclature):
    """
    Pango nomenclature in the general sense, absent SARS-CoV-2- or mpox-specific features.

    See: https://doi.org/10.1038/s41564-020-0770-5
    """

    def __init__(self, alias_map_hybrid, max_sublevels, special):
        super().__init__(
            alias_map_hybrid=alias_map_hybrid,
            charsets=[set(string.ascii_uppercase), set(string.digits)],
            cumulative_alias=False,
            max_sublevels=max_sublevels,
            root="",
            sep=r".",
            special=special,
        )
        self.ambiguity = r"*"
        self.tip = r"$"

    ##############################
    # Superclass implementations #
    ##############################

    def is_ambiguous(self, name: str) -> bool:
        """
        Does this name specify an ambiguous taxon?

        Pango taxa are ambiguous if the name ends in *, such that JN.1* means
        some unknown or unspecified sublineage of JN.1.

        Parameters
        ----------
        name : str
            The name of the taxon

        Returns
        -------
        bool
            True if the name is ambiguous.
        """
        return not self.is_root(name) and str(name)[-1] == self.ambiguity

    def is_hybrid(self, name: str) -> bool:
        """
        Does this name specify a hybrid taxon?

        Hybrids are recombinants, and recombinant names start with X:
        https://virological.org/t/pango-lineage-nomenclature-provisional-rules-for-naming-recombinant-lineages/657

        Parameters
        ----------
        name : str
            The name of the taxon

        Returns
        -------
        bool
            True if the name is a hybrid.
        """
        return not self.is_root(name) and name[0] == "X"

    def is_special(self, name: str) -> bool:
        return name in self.special or (
            self.is_hybrid(name) and self.num_sublevels(name) == 0
        )

    ########################
    # Superclass overrides #
    ########################

    def is_valid_name(
        self,
        name: str,
        min_sublevels: int = 1,
        max_sublevels: int | None = None,
    ) -> bool:
        """
        Is this name valid in the Pango nomenclature?

        A valid name must have >1 and <= self.max_sublevels sublevels unless
        it is a special-purpose ancestor such as a recombinant or a directly-
        named root descendant, in which case they may have 0 sublevels.

        Parameters
        ----------
        name : string specifying name of the taxon

        Returns
        -------
        bool
            True if this is a valid name under the Pango nomenclature.
        """
        return (
            self.is_special(name)
            or self.is_hybrid(name)
            or super().is_valid_name(name, min_sublevels, max_sublevels)
        )


class PangoSc2Nomenclature(PangoNomenclature):
    """
    Pango nomenclature for SARS-CoV-2.

    A PangoNomenclature with a specific .name() method, a known url for the
    alias map, maximally 3 sublevels, and the special root descendants A and B.

    See: https://doi.org/10.1038/s41564-020-0770-5
    """

    def __init__(
        self,
        gh_alias_url="https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json",
    ):
        """
        PangoSc2Nomenclature constructor
        """

        super().__init__(alias_map_hybrid=[list], max_sublevels=3, special=["A", "B"])
        self.gh_alias_url: str = gh_alias_url

    ##############################
    # Superclass implementations #
    ##############################
    def name(self) -> str:
        return "PangoNomenclature(SARS-CoV-2)"

    ########################
    # Superclass overrides #
    ########################

    def setup_alias_map(
        self, fp: Optional[str] = None, url: Optional[str] = None
    ) -> None:
        """
        Aliasing schemes come from either local or remote json files.

        The json is assumed to be a map from the alias to a longer form of the name (but not necessarily the longest possible form)

        Parameters
        ----------
        fp : Optional[str]
            File path of local json file containing the alias map.
        url : Optional[str]
            URL of online json file containing the alias map.
            If None, pulls from the Pango SARS-CoV-2 GitHub repo.

        Must  supply one of fp or url.

        Returns
        -------
        None
            Reads the alias map and stores it in self.alias_map, then calls
            self.sanitize_map() and self.invert_map().
        """
        if url is None:
            url = self.gh_alias_url
        super().setup_alias_map(fp=fp, url=url)
