from cladecombiner import Taxon
from cladecombiner.agg_utils import get_versioned_tip_taxa


def test_tip_clades(pango_historical_bundle):
    versioned_taxa, pango_historical = pango_historical_bundle

    current_taxa = [
        "MONTY.42",
        "MONTY.25.25.25",
        "OF.9.9.9",
        "FLYING.1.1",
        "FLYING.1.1.1",
        "FLYING.8472",
    ]

    taxa_tips = get_versioned_tip_taxa(
        pango_historical.taxonomy_tree(
            [Taxon(taxon, True) for taxon in current_taxa]
        ),
        pango_historical.get_versioner(None),
    )
    print(f">>>>>> {taxa_tips} <<<<<<")

    expected_taxa_tips = [
        ("MONTY.25.25.25", True),
        ("PYTHONS.0.0", False),
        ("PYTHONS.47.47.47", False),
    ]
    print(f">>>>>> {expected_taxa_tips} <<<<<<")

    # assert all(tt in expected_taxa_tips for tt in taxa_tips)
    for tt in taxa_tips:
        print(tt)
        assert tt in expected_taxa_tips
