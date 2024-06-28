import json
import string
import urllib.request
import warnings
from abc import ABC, abstractmethod
from collections.abc import Collection, Container, MutableSequence, Sequence
from sys import maxsize as integer_inf
from typing import Any, Callable, Optional

import dendropy

from .taxon import Taxon
from .tree_utils import add_paraphyletic_tips
from .utils import table


class Nomenclature(ABC):
    """
    Abstract class for most general casting of Nomenclature

    Nomenclature concerns rules for naming taxa, and what names may imply about
    those taxa.
    """

    @abstractmethod
    def is_ambiguous(self, name: str) -> bool:
        """
        Does this name indicate an ambiguous taxon?

        Ambiguity means a taxon specified only to a higher level than to which
        resolution is possible.

        Returns
        -------
        bool
            True if this name indicates an ambiguous taxon.
        """
        pass

    @abstractmethod
    def is_hybrid(self, name: str) -> bool:
        """
        Does this name indicate a hybrid?

        Hybrid taxa have more than one parent taxon.

        Parameters
        ----------
        name : string specifying name of the taxon

        Returns
        -------
        bool
            True if this name indicates a hybrid taxon.
        """
        pass

    @abstractmethod
    def is_root(self, name: str) -> bool:
        """
        Does this string specify the root taxon?

        The root taxon includes all taxa in the nomenclature scheme.

        Parameters
        ----------
        name : string specifying name of the taxon

        Returns
        -------
        bool
            True if this name indicates the root taxon.
        """
        pass

    @abstractmethod
    def is_valid_name(self, name: str) -> bool:
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

    @abstractmethod
    def name(self) -> str:
        """
        Name of this nomenclature scheme

        Returns
        -------
        string
            The name of this taxonomy scheme.
        """
        pass

    def __str__(self):
        return self.name()


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
    def full_histories(
        self, taxa: Sequence[str], stop_at_hybrid: bool = False
    ) -> Sequence[Sequence[str]]:
        """
        For each taxon, get the sequence of names of ancestors from the root
        to it.

        Parameters
        ----------
        taxa : Sequence[str]
            Each string is the name of one taxon for which we want the full
            history.
        stop_at_hybrid : boolean
            If True, the history for a taxon starts at the most recent
            hybridization event in its ancestry. If False, we extract a linear
            history by taking the ancestry through the first indicated parent
            every time.

        Returns
        -------
        Sequence[Sequence[str]]
            For each input taxon, the history from the root to the taxon as a
            Sequence of names of taxa.
        """
        pass

    def subtree_from_histories(
        self, node: dendropy.Node, lvl: int, histories: Sequence[Sequence[str]]
    ) -> None:
        """
        Recursive building of taxonomic tree from taxon-specific histories.

        Parameters
        ----------
        node : dendropy.Node
            Node defining the subtree to operate on.
        lvl : int
            How many levels deep from the root are we?
        histories: Sequence[Sequence[str]]
            The histories of all taxa in this subtree for which we are
            attempting to construct the subtree.

        Returns
        -------
        None
            Modifies tree in-place recursively.
        """
        next_step = set([history[lvl] for history in histories])
        if len(next_step) == 1:
            child = dendropy.Node(label=next_step.pop())
            node.add_child(child)
            next_histories = [
                history for history in histories if len(history) > lvl + 1
            ]
            if next_histories:
                self.subtree_from_histories(child, lvl + 1, next_histories)
        elif len(next_step) > 1:
            for step in next_step:
                child = dendropy.Node(label=step)
                node.add_child(child)
                next_histories = [
                    history
                    for history in histories
                    if len(history) > lvl + 1 and history[lvl] == step
                ]
                if next_histories:
                    self.subtree_from_histories(child, lvl + 1, next_histories)

    def taxonomy_tree(
        self,
        taxa: Sequence[Taxon],
        insert_tips: bool,
        name_cleanup_fun: Optional[Callable[[str], str]] = None,
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
        name_cleanup_fun : Optional[Callable]
            A function applied to all node labels after the tree is
            constructed, to ensure validity of all names.
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
            warnings.warn(
                "Removed non-unique and/or non-tip taxa to build tree."
            )

        histories = self.full_histories(unique_names)

        all_names: set[str] = set()
        for history in histories:
            for taxon in history:
                all_names.add(taxon)

        namespace = dendropy.TaxonNamespace(list(all_names))
        phy = dendropy.Tree(taxon_namespace=namespace)
        node = phy.seed_node
        if not isinstance(node, dendropy.Node):
            # Should never hit, required for type checking
            raise RuntimeError(
                "Cannot start tree because seed_node is not a dendropy.Node"
            )

        # Support for forests, where we break trees at recombination, could be added
        first_step = set([history[0] for history in histories])
        if len(first_step) != 1:
            raise RuntimeError(
                "Cannot start tree, not all histories have same root"
            )
        node.label = first_step.pop()

        self.subtree_from_histories(node, 1, histories)

        if name_cleanup_fun is not None:
            for node in phy.preorder_node_iter():
                node.label = name_cleanup_fun(node.label)

        if insert_tips:
            phy = add_paraphyletic_tips(phy, unique_names)

        tip_names = [node.label for node in phy.leaf_node_iter()]
        int_names = [node.label for node in phy.preorder_internal_node_iter()]

        if len(set(tip_names)) != len(tip_names):
            tab = table(tip_names)
            mults = ", ".join(
                [
                    str(k) + " (x" + str(v) + ")"
                    for k, v in tab.items()
                    if v > 1
                ]
            )
            raise RuntimeError(
                "Malformed tree has multiples of tip taxa: " + mults
            )

        if len(set(int_names)) != len(int_names):
            tab = table(int_names)
            mults = ", ".join(
                [
                    str(k) + " (x" + str(v) + ")"
                    for k, v in tab.items()
                    if v > 1
                ]
            )
            raise RuntimeError(
                "Malformed tree has multiples of internal taxa: " + mults
            )

        return phy


class PangoLikeNomenclature(AlgorithmicNomenclature):
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
    """

    def __init__(
        self,
        alias_map_hybrid,
        charsets,
        cumulative_alias,
        max_sublevels,
        root,
        sep,
        special,
    ):
        """ """
        self.alias_map: dict = {}
        "Defines mapping to make longer names from shorter ones"
        self.alias_map_hybrid: Collection[type] = alias_map_hybrid
        "Container type(s) used in alias map when hybrid ancestry is indicated"
        self.alias_map_inv: dict = {}
        "Defines mapping to make shorter names from longer ones"
        self.charsets: Sequence[set] = charsets
        "Defines what's allowed in alias names [0] and sublevel names [1]"
        self.cumulative_alias: bool = cumulative_alias
        "Does the alias accumulate (like RSV system) or not (like Pango)"
        self.max_sublevels: int = max_sublevels
        "Defines maximum number of sublevels before aliasing must be done"
        self.root: str = root
        "Name for the root taxon"
        self.sep: str = sep
        "Defines what separates the levels of the name"
        self.special: Container = special
        "Defines what aliases are allowed to appear alone"

    ##############################
    # Superclass implementations #
    ##############################

    def full_histories(
        self, taxa: Sequence[str], stop_at_hybrid: bool = False
    ) -> Sequence[Sequence[str]]:
        if stop_at_hybrid:
            raise NotImplementedError(
                "Forests of histories are not currently implemented or supported."
            )
        return [self.get_history(taxon, stop_at_hybrid) for taxon in taxa]

    def is_root(self, name: str) -> bool:
        return name == self.root

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

    def coax_name(self, name: str) -> str:
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
        if self.is_root(name):
            return name
        return self.shorter_name(self.longer_name(name))

    def equals_ignore_alias(self, x: str, y: str) -> bool:
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
        return self.find_ignore_alias(x, [y]) == 0

    def find_ignore_alias(self, taxon: str, taxa: Sequence[str]) -> int:
        """
        Find taxon in taxa, comparing long-form names, -1 if no match

        Assumes there is only one match. If there are multiple matches, returns
        the index of the first.

        For a Pango SARS-CoV-2 example, if taxon="KP.3" and
        taxa=["JN.1.11.1.3", "JN.1"], we will get 0, because "JN.1.11.1.3" is
        an equivalent name to "KP.3", ignoring the shortening due to aliases.

        Parameters
        ----------
        taxon : str
            The taxon to be located.
        taxa : Sequence[str]
            The taxa in which we search for the taxon

        Returns
        -------
        int
            Position of match in taxa, or -1 if no match
        """
        match = -1
        target_name = self.longer_name(taxon)
        for i in range(len(taxa)):
            if target_name == self.longer_name(taxa[i]):
                match = i
                break
        return match

    def get_history(self, name: str, stop_at_hybrid: bool) -> Sequence[str]:
        """
        Get a path of ancestry from the root to this taxon.

        This is different than a long-form name because it allows us to pass
        through hybridization (recombination) events. In the face of
        recombination, when stop_at_hybrid == False, we follow the ancestry of
        the 5'-most portion of the genome.

        Parameters
        ----------
        taxon : str
            A taxon's name.
        taxa : Sequence[str]
            The taxa in which we search for the taxon

        Returns
        -------
        int
            Position of match in taxa, or -1 if no match
        """
        if not self.alias_map:
            raise RuntimeError(
                "Cannot obtain histories until setup_alias_map() has been called."
            )
        if not self.is_valid_name(name):
            raise ValueError(name + " is not a valid name in " + self.name())
        history = []
        self.extend_history(name, history, stop_at_hybrid)
        history.reverse()
        return history

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
            raise ValueError("Invalid name: " + name)
        # Digest sublevels
        if comp[1]:
            for i in range(1, len(comp[1]) + 1)[::-1]:
                history.append(self.unpartition_name([comp[0], comp[1][:i]]))
        # Handle alias
        alias = self.join(comp[0])
        if self.is_root(alias):
            history.append(self.root)
        else:
            if self.is_special(alias):
                history.append(alias)
            if not self.is_hybrid(alias):
                self.extend_history(
                    self.alias_map[alias], history, stop_at_hybrid
                )
            elif not stop_at_hybrid:
                self.extend_history(
                    self.alias_map[alias][0], history, stop_at_hybrid
                )

    def invert_map(self) -> None:
        """
        Inverts the shorter->longer self.alias_map

        The inverted alias map is incapable of handling hybridization.

        Returns
        -------
        None
            The inverted map is stored as self.alias_map_inv
        """
        rev_map = {}
        for k, v in self.alias_map.items():
            if not isinstance(v, list):
                v = [v]
            for vi in v:
                # Don't add empty root alias
                if (not self.is_root(vi)) and (not self.is_hybrid(k)):
                    if vi in rev_map:
                        raise RuntimeError(
                            "Alias list cannot be inverted. "
                            + "Trying to add inverse alias for "
                            + vi
                            + ", which is an alias of "
                            + k
                            + ", but reversed map already has it"
                        )
                    rev_map[vi] = k
        self.alias_map_inv = rev_map

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
        for t in self.alias_map_hybrid:
            if isinstance(alias_value, t):
                return True
        return False

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
            return all([set(a) < self.charsets[0] for a in self.split(alias)])
        else:
            return set(alias) < self.charsets[0]

    def join(self, comp: Sequence[str]) -> str:
        """
        Join list of component levels into name.

        The inverse of self.split(name), such that
        self.join(self.split(name)) == name.

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
        # print("Trying to find shorter form of " + name)
        for k, v in self.alias_map_inv.items():
            kl = self.partition_name(k)
            # print(">> " + str(k) + " : " + str(v))
            # print(">>>> " + str(kl[1]) + " == " + str(parts[1]))
            # print(">>>> " + str(self.partition_name(self.longer_name(k))[1][:n]) + " == " + str(parts[1][:n]))
            if (
                kl[1] == parts[1]
                or self.partition_name(self.longer_name(k))[1][:n]
                == self.partition_name(self.longer_name(name))[1][:n]
            ):
                alias = v
                break
        if not alias:
            raise RuntimeError("Cannot find shorter alias for " + name)
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
        comp = self.split(name)
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

    def sanitize_map(self) -> None:
        """
        Drop ambiguity markers and check all names are valid.

        For the purposes of determining ancestry, an unknown sublineage is
            effectively just its ancestor, and we treat it as such.

        Returns
        -------
        None
            Modifies self.alias_map in-place
        """
        if not self.alias_map:
            raise RuntimeError(
                "Missing self.alias_map when trying to sanitize."
            )

        for k, v in self.alias_map.items():
            if not self.is_valid_alias(k):
                raise RuntimeError(
                    "Found invalid taxon as key in alias list: " + k
                )
            if self.is_ambiguous(k):
                raise RuntimeError(
                    "Found ambiguous taxon as key in alias list: " + k
                )
            if not self.is_alias_map_hybrid(v):
                if self.is_root(v):
                    if not self.is_special(k):
                        raise RuntimeError(
                            'Found alias for root in taxon not listed as special: "'
                            + k
                            + '"'
                        )
                else:
                    v = [v]
            for i in range(len(v)):
                if self.is_ambiguous(v[i]):
                    v[i] = v[i][:-1]
                if not self.is_valid_name(v[i], max_sublevels=integer_inf):
                    raise RuntimeError(
                        'Found invalid taxon as value in alias list: "'
                        + v[i]
                        + '" (for key "'
                        + k
                        + '")'
                    )

    def setup_alias_map(self, fp: Optional[str], url: Optional[str]) -> None:
        """
        Aliasing schemes come from either local or remote json files

        The json is assumed to be a map from the alias to a longer form of the name (but not necessarily the longest possible form)

        Parameters
        ----------
        fp : Optional[str]
            File path of local json file containing the alias map.
        url : Optional[str]
            URL of online json file containing the alias map.

        Must  supply one of fp or url.

        Returns
        -------
        None
            Reads the alias map and stores it in self.alias_map, then calls
            self.sanitize_map() and self.invert_map().
        """
        # Should we be thinking about encoding and/or defensive measures?
        if fp:
            alias_file = open(fp)
            alias = json.load(alias_file)
            alias_file.close()
            self.alias_map = dict(alias)
        elif url:
            with urllib.request.urlopen(url) as response:
                self.alias_map = json.loads(response.read().decode("utf8"))
        else:
            raise RuntimeError(
                "Must supply either file path or URL to alias json"
            )

        self.sanitize_map()
        self.invert_map()

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
            raise RuntimeError(
                "Cannot get shorter name without an alias list."
            )
        comp = list(self.partition_name(name))
        lvl = 1
        while len(comp[1]) > self.max_sublevels:
            alias = self.next_shorter_alias(self.unpartition_name(comp), lvl)
            comp[0] = [alias]
            comp[1] = comp[1][3:]
            lvl += 1
        return self.unpartition_name(comp)

    def split(self, name: str) -> Sequence[str]:
        """
        Split name into component levels

        The inverse of self.join(name), such that
        self.split(self.join(components)) == components.

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
        return self.join([*components[0], *components[1]])


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
        if self.is_root(name):
            return False
        elif str(name)[-1] == self.ambiguity:
            return True
        else:
            return False

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
        if self.is_root(name):
            return False
        elif name[0] == "X":
            return True
        else:
            return False

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
        if self.is_special(name) or self.is_hybrid(name):
            return True
        return super().is_valid_name(name, min_sublevels, max_sublevels)


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

        super().__init__(
            alias_map_hybrid=[list], max_sublevels=3, special=["A", "B"]
        )
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
