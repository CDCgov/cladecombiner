# import csv
# import urllib.request
# from collections.abc import Callable

# import dendropy

# from .nomenclature import ArbitraryNomenclature
# from .taxonomy_scheme import PhylogeneticTaxonomyScheme


# def arbitrary_taxonomy_helper(
#     clade_hierarchy: dict | str,
#     expected_root: str,
#     name: str,
#     ambiguous_fun: Callable = lambda _: False,
#     hybrid_fun: Callable = lambda _: False,
# ) -> tuple[ArbitraryNomenclature, PhylogeneticTaxonomyScheme]:
#     """
#     Constructs a Nomenclature and TaxonomyScheme from a {child : parent} dict of names of taxa.

#     Parameters
#     ----------
#     clade_hierarchy : Container[str]
#         The names of all taxa recognized by this nomenclature system.
#     expected_root: str
#         The expected name of the taxon that includes all others, used for checking results.
#     name : str
#         The name of the resultant nomenclature system.
#     ambiguous_fun : Callable
#         As in cladecombiner.ArbitraryNomenclatureScheme. Default is sensible for schemes that cannot specify ambiguity.
#     hybrid_fun : Callable
#         As in cladecombiner.ArbitraryNomenclatureScheme. Default is sensible for schemes that cannot specify hybrids.

#     Returns
#     -------
#     ArbitraryNomenclature, PhylogeneticTaxonomyScheme
#         An ArbitraryNomenclature and a PhylogeneticTaxonomyScheme for use with Nextstrain clades.
#     """
#     if isinstance(clade_hierarchy, dict):
#         child_parent = clade_hierarchy
#     elif isinstance(clade_hierarchy, str):
#         with urllib.request.urlopen(clade_hierarchy) as response:
#             tsv_str = response.read().decode("utf8")
#             tsv_reader = csv.DictReader(tsv_str.split("\n"), delimiter="\t")
#         child_parent = {row["clade"]: row["parent"] for row in tsv_reader}
#     else:
#         raise TypeError(
#             "Argument `clade_hierarchy` should specify either the clade hierarchy as a dictionary or "
#         )

#     taxonomy = edge_dict_to_taxonomy(child_parent)
#     assert taxonomy.tree.seed_node is not None

#     nomenclature = ArbitraryNomenclature(
#         known_taxa={
#             node.label for node in taxonomy.tree.postorder_node_iter()
#         },
#         root=taxonomy.tree.seed_node.label,
#         name=name,
#         ambiguous_fun=ambiguous_fun,
#         hybrid_fun=hybrid_fun,
#     )

#     assert (
#         nomenclature.root == expected_root
#     ), f"Found unexpected root {nomenclature.root} instead of {expected_root}"

#     return nomenclature, taxonomy


# def get_nextstrain_sc2(
#     clade_hierarchy: dict
#     | str = "https://raw.githubusercontent.com/nextstrain/ncov/refs/heads/master/defaults/clade_hierarchy.tsv",
# ) -> tuple[ArbitraryNomenclature, PhylogeneticTaxonomyScheme]:
#     """
#     Get Nomenclature and TaxonomyScheme for use with Nextstrain clades for SARS-CoV-2.

#     Parameters
#     ---------
#     clade_hierarchy : dict[str, list[str]]
#         Either a dictionary specifying for every clade its parent clade, or the
#         path to clade_hierarchy.tsv which specifies such.

#     Returns
#     -------
#     ArbitraryNomenclature, PhylogeneticTaxonomyScheme
#         An ArbitraryNomenclature and a PhylogeneticTaxonomyScheme for use with Nextstrain clades.
#     """
#     return arbitrary_taxonomy_helper(
#         clade_hierarchy=clade_hierarchy,
#         expected_root="19A",
#         name="NextstrainClades(SARS-CoV-2)",
#     )
