import pytest

from cladecombiner import (
    BasicPhylogeneticAggregator,
    PhylogeneticTaxonomyScheme,
    Taxon,
    read_taxa,
)


@pytest.fixture
def expected_one_agg():
    expected = {}
    expected["LIFE"] = ["LIFE.1", "LIFE.1.2", "LIFE.123", "LIFE.7"]
    expected["FLYING.8472"] = [
        "FLYING.8472.1",
        "FLYING.8472.1.2",
        "FLYING.8472.123",
        "FLYING.8472.7",
    ]
    expected["FLYING.123"] = [
        "FLYING.123.1",
        "FLYING.123.1.2",
        "FLYING.123.123",
        "FLYING.123.7",
    ]
    return expected


def test_no_drop_no_overlap(pango_with_toy_alias):
    targets = [
        Taxon("FLYING.8472", False),
        Taxon("FLYING.123", False),
        Taxon("PYTHONS.0.0.0", False),
    ]

    input_taxa = read_taxa(
        "tests/toy_lineages.txt",
        is_tip=True,
        nomenclature=pango_with_toy_alias,
    )
    tree = pango_with_toy_alias.taxonomy_tree(input_taxa)
    taxonomy_scheme = PhylogeneticTaxonomyScheme(tree)
    agg = BasicPhylogeneticAggregator(input_taxa, targets, taxonomy_scheme)

    expected = {
        "LIFE.1": "PYTHONS.0.0.0",
        "LIFE.1.2": "PYTHONS.0.0.0",
        "LIFE.123": "PYTHONS.0.0.0",
        "LIFE.7": "PYTHONS.0.0.0",
        "FLYING.8472.1": "FLYING.8472",
        "FLYING.8472.1.2": "FLYING.8472",
        "FLYING.8472.123": "FLYING.8472",
        "FLYING.8472.7": "FLYING.8472",
        "FLYING.123.1": "FLYING.123",
        "FLYING.123.1.2": "FLYING.123",
        "FLYING.123.123": "FLYING.123",
        "FLYING.123.7": "FLYING.123",
    }
    print("+++++++++++++++++++++++++++++++++++++++++")
    print(agg.map())
    print("+++++++++++++++++++++++++++++++++++++++++")
    assert agg.map() == expected
