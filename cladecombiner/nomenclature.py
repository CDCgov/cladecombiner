import json
import string
import urllib.request
from abc import ABC, abstractmethod
from collections.abc import Collection, Container, Sequence
from sys import maxsize as integer_inf
from typing import Any, Optional

import dendropy


class Nomenclature(ABC):
    """
    Abstract class for most general casting of Nomenclature
    """

    @abstractmethod
    def is_ambiguous(self, name: str) -> bool:
        """
        Does this name indicate an ambiguous taxon?

        Ambiguity means a taxon known only to a higher level than to which resolution is possible.
        """
        pass

    @abstractmethod
    def is_hybrid(self, name: str) -> bool:
        """Does this name indicate a hybrid?"""
        pass

    @abstractmethod
    def is_root(self, name: str) -> bool:
        """Does this string specify the root taxon (to which all other taxa belong)?"""
        pass

    @abstractmethod
    def is_valid_name(self, name: str) -> bool:
        """Is this name valid in the nomenclature scheme?"""
        pass

    @abstractmethod
    def name(self) -> str:
        """Name of this nomenclature scheme"""
        pass

    def __str__(self):
        return self.name()


class AlgorithmicNomenclature(Nomenclature):
    """
    Abstract class Nomenclature schemes which encode a taxon's history in its name in some form.
    """

    @abstractmethod
    def full_histories(
        self, taxa: Sequence[str], stop_at_hybrid: bool = False
    ) -> Sequence[Sequence[str]]:
        """
        For each taxon, get the sequence of names of ancestors from the root to it.

        Arguments:
        taxa: names of taxa for which we want the full histories
        stop_at_hybrid: If True, the history for a taxon starts at the most recent hybridization event in its ancestry.
        If False, we extract a linear history by taking the ancestry through the first indicated parent every time.
        """
        pass

    def add_tips_for_ancestral(
        self, phy: dendropy.Tree, tips: Sequence[str]
    ) -> None:
        raise NotImplementedError()

    def subtree_from_histories(
        self, node: dendropy.Node, lvl: int, histories: Sequence[Sequence[str]]
    ) -> None:
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
        self, taxa: Sequence[str], insert_tips: bool = True
    ) -> dendropy.Tree:
        """ """
        histories = self.full_histories(taxa)

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

        first_step = set([history[0] for history in histories])
        if len(first_step) != 1:
            raise RuntimeError(
                "Cannot start tree, not all histories have same root"
            )
        node.label = first_step.pop()

        self.subtree_from_histories(node, 1, histories)

        if insert_tips:
            self.add_tips_for_ancestral(phy, taxa)

        return phy


class PangoLikeNomenclature(AlgorithmicNomenclature):
    """
    A Pango-like nomenclature is an algorithmic nomenclature with more specific assumptions about the encoding of history.

    has two parts to any name.
    The first is an aliasing portion which is a shorter stand-in for a longer string describing its history.
    The second is a series of sublevels relative to the aliased starting point.
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

    #######
    # Implementations/specializations of super-class methods
    #######

    def full_histories(
        self, taxa: Sequence[str], stop_at_hybrid: bool = False
    ) -> Sequence[Sequence[str]]:
        long_names = [
            self.longer_name(taxon, stop_at_hybrid=stop_at_hybrid)
            for taxon in taxa
        ]
        histories = []
        for name in long_names:
            history = [""]
            levels = self.split(name)
            for i in range(len(levels)):
                history.append(self.join(levels[: (i + 1)]))
            histories.append(history)
        return histories

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

    def taxonomy_tree(self, taxa: Sequence[str]) -> dendropy.Tree:
        raise NotImplementedError()
        # phy = super().taxonomy_tree(taxa)
        # clean up names

    #######
    # Implementations of class methods
    #######

    def equals_ignore_alias(self, x: str, y: str) -> bool:
        """
        Are two names the same, accounting for aliasing?
        """
        return self.find_ignore_alias(x, [y]) == 1

    def find_ignore_alias(self, taxon: str, taxa: Sequence[str]) -> int:
        """
        Find taxon in taxa, comparing long-form names, -1 if no match
        """
        match = -1
        target_name = self.longer_name(taxon)
        for i in range(len(taxa)):
            if target_name == self.longer_name(taxa[i]):
                match = i
                break
        return match

    def invert_map(self) -> None:
        """Inverts the shorter->longer self.alias_map"""
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

    def is_alias_map_hybrid(self, parent_info: Any) -> bool:
        for t in self.alias_map_hybrid:
            if isinstance(parent_info, t):
                return True
        return False

    def is_special(self, name: str) -> bool:
        """
        Does this string specify a name that is allowed not to have further levels?
        """
        return name in self.special

    def is_valid_alias(self, alias: str) -> bool:
        """
        Does this string specify a valid key for aliasing?
        """
        if self.cumulative_alias:
            return all([set(a) < self.charsets[0] for a in self.split(alias)])
        else:
            return set(alias) < self.charsets[0]

    def join(self, comp: Sequence[str]) -> str:
        """
        Join list of component levels into name
        """
        return self.sep.join(comp)

    # def longer_name(self, name: str, stop_at_hybrid: bool = True) -> str:
    #     """
    #     Get non-aliased form of an aliased name
    #     """
    #     if not self.alias_map:
    #         raise RuntimeError(
    #             "Cannot construct long form of name without an alias list."
    #         )
    #     alias_levels = list(self.partition_name(name))
    #     next_alias = self.alias_map[alias_levels[0][-1]]
    #     while (not self.is_root(next_alias)) and (not self.is_hybrid(alias_levels[0][-1])):
    #         parts = self.partition_name(next_alias)
    #         alias_levels[0] = parts[0]
    #         alias_levels[1] = [*parts[1], *alias_levels[1]]
    #         next_alias = self.alias_map[parts[0][-1]]
    #     return self.join([*alias_levels[0], *alias_levels[1]])

    def longer_name(self, name: str, stop_at_hybrid: bool = True) -> str:
        """
        Get non-aliased form of an aliased name
        """
        if not self.alias_map:
            raise RuntimeError(
                "Cannot construct long form of name without an alias list."
            )
        alias_levels = list(self.partition_name(name))
        next_alias = self.alias_map[alias_levels[0][-1]]
        while not self.is_root(next_alias):
            if self.is_hybrid(alias_levels[0][-1]):
                if stop_at_hybrid:
                    break
                else:
                    next_alias = next_alias[0]
            parts = self.partition_name(next_alias)
            alias_levels[0] = parts[0]
            alias_levels[1] = [*parts[1], *alias_levels[1]]
            next_alias = self.alias_map[parts[0][-1]]
        return self.unpartition_name(alias_levels)

    def next_shorter_alias(self, name: str, depth: int) -> str:
        """
        Get the next shortest name available to a taxon

        Arguments:
            name: an expanded taxon name to be contracted
            depth: how many levels of aliasing deep is this name? Starting at 1 for longest (fully de-aliased) name and increasing as the name gets shorter
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

    def partition_name(self, name: str) -> Sequence[Sequence[str]]:
        """
        Splits name into alias and sublevels, each as a sequence of components

        This function assumes that the name is ordered alias, sublevels, and does not check correctness.
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
                    if k not in self.special:
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
        """
        # @TODO should we be thinking about encoding and/or defensive measures?
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
        """Get shortest form of a maximally-long name using aliases"""
        if not self.alias_map:
            raise RuntimeError(
                "Cannot get shorter name without an alias list."
            )

        alias_levels = list(self.partition_name(name))
        lvl = 1
        while len(alias_levels[1]) > self.max_sublevels:
            alias = self.next_shorter_alias(
                self.unpartition_name(alias_levels), lvl
            )
            alias_levels[0] = [alias]
            alias_levels[1] = alias_levels[1][3:]
            lvl += 1
        return self.unpartition_name(alias_levels)

    def split(self, name: str) -> Sequence[str]:
        """Split name into component levels"""
        return name.split(self.sep)

    def unpartition_name(self, components: Sequence[Sequence[str]]) -> str:
        return self.join([*components[0], *components[1]])


class PangoNomenclature(PangoLikeNomenclature):
    """
    Pango nomenclature in the general sense, absent SARS-CoV-2- or mpox-specific rules

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

    #####
    # Implementation of superclass methods
    #####

    def is_ambiguous(self, name: str) -> bool:
        if self.is_root(name):
            return False
        elif str(name)[-1] == self.ambiguity:
            return True
        else:
            return False

    def is_hybrid(self, name: str) -> bool:
        # Hybrids are recombinants, and recombinant names start with X
        # https://virological.org/t/pango-lineage-nomenclature-provisional-rules-for-naming-recombinant-lineages/657
        # Root has name "" with len(0)
        if self.is_root(name):
            return False
        elif name[0] == "X":
            return True
        else:
            return False

    ########################
    # Superclass overrides #
    ########################
    def is_special(self, name):
        # Under the Pango scheme, recombinants are special-purpose ancestors
        return name in self.special or self.is_hybrid(name)


class PangoSc2Nomenclature(PangoNomenclature):
    def __init__(
        self,
        gh_alias_url="https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json",
    ):
        super().__init__(
            alias_map_hybrid=[list], max_sublevels=3, special=["A", "B"]
        )
        self.gh_alias_url: str = gh_alias_url

    #####
    # Implementation of superclass methods
    #####
    def name(self) -> str:
        return "PangoNomenclature(SARS-CoV-2)"

    #####
    # Overrides of superclass methods
    #####
    def is_valid_name(
        self,
        name: str,
        min_sublevels: int = 1,
        max_sublevels: int | None = None,
    ) -> bool:
        if self.is_hybrid(name) or self.is_special(name):
            min_sublevels = 0
        return super().is_valid_name(name, min_sublevels, max_sublevels)

    def setup_alias_map(self, **kwargs) -> None:
        fp = None
        if "fp" in kwargs:
            fp = kwargs["fp"]
        super().setup_alias_map(fp=fp, url=self.gh_alias_url)
