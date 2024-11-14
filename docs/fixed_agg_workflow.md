# Phylogenetic aggregation to fixed targets

Using the `BasicPhylogeneticAggregator`, we can phylogenetically aggregate a set of observed taxa into a pre-defined set of aggregated taxa.

## Observed taxa

First, we need to get the observed taxa as `cladecombiner.Taxon` objects.
If these are in a file, we can read them in:

```
import cladecombiner

taxa = cladecombiner.read_taxa("path/to/file")
```

If we already have these as a list of strings, we can convert them readily:

```
import cladecombiner

lineages = [
    "BA.2",
    "BA.2.86",
    "XCU",
    "XDQ",
    "BQ.1.1.23",
    "JN.6",
    "JN.1.9.1",
    "JN.9",
    "KP.1.1",
    "FW.1.1.1",
    "BA.5.2.6",
]

taxa = [cladecombiner.Taxon(lin, True) for lin in lineages]
```

Note that we are declaring these to all be _tip_ taxa.
In cladecombiner parlance, we use this to distinguish taxa which we have observed directly from those we have not.
Observed taxa are tips in the tree of all relevant taxa, while all other taxa are internal nodes.
The `read_taxa` utility function defaults to the assumption that all taxa are tips, while creating `Taxon` objects requires specifying either way.

## Creating the taxonomy scheme

Aggregating phylogenetically requires us to create a `PhylogeneticTaxonomyScheme`.
We do this by using Pango nomenclature's rules to get the taxonomic tree relating the observed lineages, and creating a taxonomy scheme from it.
We need access the alias json for Pango SARS-CoV-2 lineages so that we know that, for example, BA is B.1.1.529.
By default, cladecombiner automatically retrieves the latest version from the [pango-designation](https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json) git repo.

```
from cladecombiner import pango_sc2_nomenclature as pn

pn.setup_alias_map()
tree = pn.taxonomy_tree(taxa)

scheme = cladecombiner.PhylogeneticTaxonomyScheme(tree)
```

These trees are wrappers around `dendropy.Tree` objects, and can be printed to screen for inspection.
The command `print(tree.as_ascii_plot(plot_metric="level", show_internal_node_labels=False))` will, for these taxa, yield

```
                         /----+---- BA.5.2.6
                    /----+
                    |    \----+----+----+----+----+----+----+---- BQ.1.1.23
                    |
                    |                   /----+----+----+---- KP.1.1
                    |              /----+
                    |              |    \----+---- JN.1.9.1
                    |              |
+----+----+----+----+         /----+---- XDQ
                    |         |    |
                    |         |    |---- JN.9
                    |    /----+    |
                    |    |    |    \---- JN.6
                    |    |    |
                    |    |    \---- BA.2.86
                    \----+
                         |----+----+----+----+----+----+----+----+----+---- FW.1.1.1
                         |
                         |----+----+----+----+---- XCU
                         |
                         \---- BA.2
```

The observed taxa are all tips, and the internal nodes are represented by `+` (since we have set `show_internal_node_labels=False`, otherwise they would be named and the tree would get quite wide).
Some of these nodes here have only one child, while others have many, because cladecombiner makes this tree to include every named Pango lineage in the complete histories of all lineages.
Thus, the first four nodes are the root (in cladecombiner parlance, an empty string, `""`), B, B.1, B.1.1, and B.1.1.529 (aka BA, but note that cladecombiner is pedantic about names and does not accept naked aliases like BA).
Of these, only BA has multiple children among the observed lineages, so only it is a branching point in the tree.

## Aggregation targets

We also need the desired aggregation targets, the taxa to which we want to aggregate.
As with observed taxa, we can read these in or create them from strings.
Let us assume we want to aggregate to BA.2.86.1 (aka JN), BA.2, BA.5, and KP.1.1.

Aggregation requires us to be explicit about whether we are truly aggregating or not.
Three of these represent actual aggregations of the observed taxa, as we have observed lineages which are children of BA.2, BA.5, and JN.
However, as we have no children of KP.1.1 in our observed lineages, it is a tip.

```
target_taxa = [
    cladecombiner.Taxon("BA.2", False),
    cladecombiner.Taxon("BA.2.86.1", False),
    cladecombiner.Taxon("BA.5", False),
    cladecombiner.Taxon("KP.1.1", True),
]
```

## Aggregate

We are now ready to create an aggregation object and aggregate.

```
agg = cladecombiner.BasicPhylogeneticAggregator(targets=target_taxa, taxonomy_scheme=scheme)

res = agg.aggregate(input_taxa)
```

The `res` object is a python dictionary (a `dict[Taxon, Taxon]`) which maps each of the observed lineages in `taxa` to some aggregated taxon.
For the above input taxa and targets, the resulting mapping is:

```
Taxon(KP.1.1, tip=True)    : Taxon(KP.1.1, tip=True)
Taxon(BA.5.2.6, tip=True)  : Taxon(BA.5, tip=False)
Taxon(BQ.1.1.23, tip=True) : Taxon(BA.5, tip=False)
Taxon(JN.1.9.1, tip=True)  : Taxon(BA.2.86.1, tip=False)
Taxon(XDQ, tip=True)       : Taxon(BA.2.86.1, tip=False)
Taxon(JN.9, tip=True)      : Taxon(BA.2.86.1, tip=False)
Taxon(JN.6, tip=True)      : Taxon(BA.2.86.1, tip=False)
Taxon(XCU, tip=True)       : Taxon(BA.2, tip=False)
Taxon(BA.2.86, tip=True)   : Taxon(BA.2, tip=False)
Taxon(FW.1.1.1, tip=True)  : Taxon(BA.2, tip=False)
Taxon(BA.2, tip=True)      : Taxon(BA.2, tip=False)
```

## Important settings

When creating the `BasicPhylogeneticAggregator`, there are two important arguments: `off_target` and `sort_clades`.

### `off_target`

The mapping for a taxon which is not in the `targets` is determined by `off_target`, which can be:

- "other" to put all such taxa into `Taxon(other)`. **This is the default.**
- "self" to map all such taxa to themselves.

For example, if we hadn't included BA.2 in the targets, then the resulting map would have been, by default,

```
Taxon(KP.1.1, tip=True)    : Taxon(KP.1.1, tip=True)
Taxon(BA.5.2.6, tip=True)  : Taxon(BA.5, tip=False)
Taxon(BQ.1.1.23, tip=True) : Taxon(BA.5, tip=False)
Taxon(JN.1.9.1, tip=True)  : Taxon(BA.2.86.1, tip=False)
Taxon(XDQ, tip=True)       : Taxon(BA.2.86.1, tip=False)
Taxon(JN.9, tip=True)      : Taxon(BA.2.86.1, tip=False)
Taxon(JN.6, tip=True)      : Taxon(BA.2.86.1, tip=False)
Taxon(XCU, tip=True)       : Taxon(other, tip=False)
Taxon(BA.2.86, tip=True)   : Taxon(other, tip=False)
Taxon(FW.1.1.1, tip=True)  : Taxon(other, tip=False)
Taxon(BA.2, tip=True)      : Taxon(other, tip=False)
```

### `sort_clades`

The order in which clades are processed is important and is determined by `sort_clades`. When we process target a clade, we map all of its children to it, and then remove them from the pool of unmapped taxa. Thus, when nested target taxa are present, the order in which taxa are processed determines the mapping.

- If `sort_clades` is True, then we process clades so that, when clade X contains clade Y, we always aggregate to Y before X. Thus, we process KP.1.1 before BA.2.86.1, and BA.2.86.1 before BA.1. **This is the default.**
- If `sort_clades` is False, then we process clades in the order they are listed. As BA.2 is listed before BA.2.86.1 in the above declaration of `target_taxa`, all children of BA.2.86.1 would end up mapped to BA.2, and the aggregation would be
  ```
  Taxon(XCU, tip=True)       : Taxon(BA.2, tip=False)
  Taxon(BA.2.86, tip=True)   : Taxon(BA.2, tip=False)
  Taxon(FW.1.1.1, tip=True)  : Taxon(BA.2, tip=False)
  Taxon(JN.1.9.1, tip=True)  : Taxon(BA.2, tip=False)
  Taxon(BA.2, tip=True)      : Taxon(BA.2, tip=False)
  Taxon(XDQ, tip=True)       : Taxon(BA.2, tip=False)
  Taxon(JN.9, tip=True)      : Taxon(BA.2, tip=False)
  Taxon(JN.6, tip=True)      : Taxon(BA.2, tip=False)
  Taxon(KP.1.1, tip=True)    : Taxon(BA.2, tip=False)
  Taxon(BA.5.2.6, tip=True)  : Taxon(BA.5, tip=False)
  Taxon(BQ.1.1.23, tip=True) : Taxon(BA.5, tip=False)
  ```
