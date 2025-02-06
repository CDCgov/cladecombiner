# Other nomenclatures and taxonomies

## Nextstrain SARS-CoV-2 clades

In addition to user-ready Pango nomenclature and taxonomy tools for SARS-CoV-2 (as featured throughout the documentation), cladecombiner provides user-ready tools for Nextstrain clades of SARS-CoV-2.
The workflow is parallel to that for Pango lineages for SARS-CoV-2, starting with a `Nomenclature` and then defining a `PhylogeneticTaxonomyScheme` for the relevant taxa.

Unlike `cladecombiner.pango_sc2_nomenclature` which has some useful functionality without an internet connection, `cladecombiner.nextstrain_sc2_nomenclature` must be able to access the `nextstrain/ncov` repo (via [PyGithub](https://pygithub.readthedocs.io/en/stable/#)) in order to be used for essentially all tasks.
There will be some delay the first time this object is used as files are read from the repo.

We can use the `Nomenclature` for basic tasks, like validating names.

```
import cladecombiner

ns = cladecombiner.nextstrain_sc2_nomenclature

ns.is_valid_name("JN.1")
# False
```

We can also get taxonomy trees from the nomenclature, such as for all the 2024 clades.
```
taxa = [
    cladecombiner.Taxon(taxon, True)
    for taxon in ['24A', '24B', '24C', '24D', '24E', '24F', '24G', '24H', '24I']
]

tree = ns.taxonomy_tree(taxa)
```

We can view the tree with `print(tree.as_ascii_plot(plot_metric="level", show_internal_node_labels=True))`, yielding:
```
                                                                        /-------- 24E
                                                               /--------24C
                                                               |        \-------- 24C
                                                      /--------24B
                                                      |        |-------- 24G
                                                      |        |
                                                      |        \-------- 24B
                                                      |
                                    /--------23I------24A------ 24F
                                    |                 |
                                    |                 |-------- 24H
                                    |                 |
19A------20A------20B------21M------21L               |-------- 24I
                                    |                 |
                                    |                 \-------- 24A
                                    |
                                    \-------- 24D
```
Note that tips for 24A, 24B, and 24C have been added to distinguish between the ancestral taxon and observed taxa.

We can use this tree to construct a `PhylogeneticTaxonomyScheme` which we can then use in aggregation.

```
import datetime

scheme = cladecombiner.PhylogeneticTaxonomyScheme(tree)

agg = cladecombiner.AsOfAggregator(scheme, ns, datetime.date(2024, 6, 1))

res = agg.aggregate(taxa)
```

This yields:

```
Taxon(24A, tip=True) : Taxon(24A, tip=False)
Taxon(24B, tip=True) : Taxon(24B, tip=False)
Taxon(24C, tip=True) : Taxon(24B, tip=False)
Taxon(24D, tip=True) : Taxon(21L, tip=False)
Taxon(24E, tip=True) : Taxon(24B, tip=False)
Taxon(24F, tip=True) : Taxon(24A, tip=False)
Taxon(24G, tip=True) : Taxon(24B, tip=False)
Taxon(24I, tip=True) : Taxon(24A, tip=False)
Taxon(24H, tip=True) : Taxon(24A, tip=False)
```
