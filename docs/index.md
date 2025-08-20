# cladecombiner

This package provides functionality for working with taxa, with a focus on phylogenetic aggregation to higher-order taxa.

While the functionality is intended to be general, current out-of-the-box support is limited to [Pango lineages](https://www.nature.com/articles/s41564-020-0770-5) and [Nextstrain clades](https://nextstrain.org/blog/2021-01-06-updated-SARS-CoV-2-clade-naming) for SARS-CoV-2.

## Taxon, clade, or lineage?

For the purposes of cladecombiner and its documentation, a taxon is a unit of taxonomy. Nextstrain clades and Pango lineages for SARS-CoV-2 are both examples of taxa.
We reserve the terms "clade" and "lineage" (unless referring to Nextstrain clades or Pango lineages) for the usual sense of monophyletic groupings defined by an ancestor and all its descendants.

As clades can contain clades and lineages can contain lineages, taxa can contain taxa.
Cladecombiner assumes that there is some base level of taxonomy which can be represented as a phylogeny with every taxon at this level as a tip.
Cladecombiner defaults to monophyletic groupings of these taxa where possible, but allows for non-monophyletic aggregated taxa.

### Ancestral versus tip taxa

When tracking ongoing evolution, the naming of new monophyletic taxa creates potential ambiguity when using the name of the now non-monophyletic ancestor.
For example, JN.1 can refer to the ancestral taxon, which has children such as KP.2.
It can also refer to JN.1 the taxon we can observe in data (a sequence which belongs to JN.1 but not any more specifically named descendant).

Cladecombiner distinguishes between these two senses in which a taxon name can be used by adding tips to the underlying taxonomy tree to represent the non-monophyletic remainder.
Thus, if needed, the ancestral `Taxon("JN.1", is_tip=False)` can have as a descendant tip `Taxon("JN.1", is_tip=True)`.

## Usage

If we have line list data like this:

| Lineage  | Location | Date       |
| -------- | -------- | ---------- |
| BA.2     | IA       | 2233-03-22 |
| BA.2.86  | IA       | 2233-03-22 |
| JN.1.9.1 | IA       | 2233-03-22 |
| KP.1.1   | IA       | 2233-03-22 |
| FW.1.1.1 | IA       | 2233-03-22 |

cladecombiner is a tool that provides mappings so that we can make it look like this:

| Lineage | Location | Date       |
| ------- | -------- | ---------- |
| BA.2    | IA       | 2233-03-22 |
| BA.2    | IA       | 2233-03-22 |
| JN.1    | IA       | 2233-03-22 |
| JN.1    | IA       | 2233-03-22 |
| JN.1    | IA       | 2233-03-22 |

Note that cladecombiner _does not_ manipulate data tables itself.
Instead, it provides taxon to taxon mappings which are generally suitable for this purpose in the form of python dictionaries.
Thus, cladecombiner strives to be useful anywhere lineages need to be aggregated.

At the moment, cladecombiner only supports aggregation when the final desired set of lineages is known.
That is, it cannot currently choose what to aggregate to on its own.
For more on doing this, see the documentation on [phylogenetic aggregation of Pango lineages](fixed_agg_workflow.md).

## Taxonomy and nomenclature

In cladecombiner parlance, nomenclature refers to rules for naming taxa while taxonomy specifies how taxa are related.
Thus, if we want to check that all of the taxa in our data are valid Pango lineages for SARS-CoV-2, we use a `Nomenclature` object, while if we want to know whether JN.1 is ancestral to KP.2, we need to use a `TaxonomyScheme` object.
We make this split, despite the fact that Pango names encode evolutionary history, for clean code separation.

### Nomenclature

Working with `Nomenclature` objects, we can check lineage name validity, get longer form names, and such.

```
import cladecombiner
from cladecombiner import pango_sc2_nomenclature as pn

pn.setup_alias_map()

my_lineages = ["JN.1", "KP.2"]

pn.validate(my_lineages)
# returns None, indicating validity

pn.longer_name("EG.1")
# 'XBB.1.9.2.1'
```

Cladecombiner distinguishes the general Pango nomenclature from the specific instance of Pango nomenclature for SARS-CoV-2.
Other instances for other outbreaks, have different [special taxa](api.md#cladecombiner.nomenclature.PangoLikeNomenclature.is_special), different numbers of sublevels, and alias maps stored in different locations.
For example, [MPXV](https://github.com/mpxv-lineages/lineage-designation).
Instances of other Pango nomenclature schemes can be made from the [`PangoNomenclature` class](api.md#cladecombiner.nomenclature.PangoNomenclature) by specificying this (and some additional) information.

### Taxonomy

To ask questions about taxonomy, we need a `TaxonomyScheme`.
We can obtain one for the taxa of interest by getting the taxonomy tree describing their relationships, and creating a `PhylogeneticTaxonomyScheme` from it.

Note that where `Nomenclature`s work on `str`s, `TaxonomyScheme`s work on `Taxon` objects.
Every `Taxon` must either be declared as a tip (something observable directly) or not.
Failure to appropriately specify this can lead to answers which seem nonsensical but which are in fact pedantically correct.
For more, see [above](#ancestral-versus-tip-taxa).

Note also that `PhylogeneticTaxonomyScheme`s only know about the taxa from which they were created and their ancestors.
That is, cladecombiner does not create the tree of all Pango SARS-CoV-2 lineages at any point unless we explicitly feed in an exhaustive list of all known lineages.

```
tree = pn.taxonomy_tree(my_taxa)
taxonomy_scheme = cladecombiner.PhylogeneticTaxonomyScheme(tree)

taxonomy_scheme.contains(
    cladecombiner.Taxon("JN.1", is_tip=False),
    cladecombiner.Taxon("KP.2", is_tip=True)
)
# True, ancestral JN.1 contains KP.1

taxonomy_scheme.contains(
    cladecombiner.Taxon("JN.1", is_tip=True),
    cladecombiner.Taxon("KP.2", is_tip=True)
)
# False, not-more-specifically-named JN.1 does not contain KP.1

taxonomy_scheme.contains(
    cladecombiner.Taxon("JN.1", is_tip=False),
    cladecombiner.Taxon("KP.1", is_tip=True)
)
# False, because KP.1 is unknown to taxonomy_scheme
```
