# As-of aggregation

As-of aggregation is a tool that can be used to approximate the taxon labels that would have been assigned to sequences on a prior (as-of) date given the taxon labels assigned on a later date.
Without loss of generality, take the later date to be the present.
Historical aggregation approximation takes in
- current taxa (whose labels were assigned with the current version of the assignment tool) as `input_taxa`
- the current taxonomy tree (in the form of a `PhylogeneticTaxonomyScheme`)
- information regarding what taxon names were recognized on the as-of date (this comes from a `HistoryAwareNomenclature`)

It outputs an `Aggregation` mapping each current taxon to the most recent ancestor which was recognized on the as-of date.

⚠️ If it is feasible to re-run sequence assignment on sequence data from the as-of date using the appropriate version of the assignment tool for that date, it is recommended that you do so.
(See for example [Cladetime](https://cladetime.readthedocs.io/en/latest/).) ⚠️

## In practice

We will continue with the taxa from [](fixed_agg_workflow/#observed-taxa)

```
import cladecombiner
from cladecombiner import pango_sc2_nomenclature as pn

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

pn.setup_alias_map()
tree = pn.taxonomy_tree(taxa)

scheme = cladecombiner.PhylogeneticTaxonomyScheme(tree)
```

Let us take our as-of date to be 2023-01-01.
In fact, we really need to specify a date and time (historical information is obtained through git histories using time stamps), so let us use midnight 2023-01-02.
Note that the call to set up the `AsOfAggregator` requires an internet connection (so that we can search the history of the public repository containing the data).
```
import datetime
agg = cladecombiner.AsOfAggregator(taxonomy_scheme=scheme, versioning_provider=pn, as_of=datetime.datetime(2023, 1, 2))
```

We can now aggregate our input taxa.

```
res = agg.aggregate(taxa)
```

This yields the following `Aggregation`:

```
Taxon(BA.2, tip=True) : Taxon(BA.2, tip=True)
Taxon(FW.1.1.1, tip=True) : Taxon(XBB.1, tip=False)
Taxon(XCU, tip=True) : Taxon(XBC.1, tip=False)
Taxon(BA.5.2.6, tip=True) : Taxon(BA.5.2.6, tip=True)
Taxon(BQ.1.1.23, tip=True) : Taxon(BQ.1.1.23, tip=True)
Taxon(KP.1.1, tip=True) : Taxon(B.1.1.529, tip=False)
Taxon(JN.6, tip=True) : Taxon(B.1.1.529, tip=False)
Taxon(JN.1.9.1, tip=True) : Taxon(B.1.1.529, tip=False)
Taxon(JN.9, tip=True) : Taxon(B.1.1.529, tip=False)
Taxon(XDQ, tip=True) : Taxon(B.1.1.529, tip=False)
Taxon(BA.2.86, tip=True) : Taxon(B.1.1.529, tip=False)
```
