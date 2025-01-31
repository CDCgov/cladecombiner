import pytest

from cladecombiner import (  # HistoricalAggregator,
    ArbitraryAggregator,
    BasicPhylogeneticAggregator,
    PhylogeneticTaxonomyScheme,
    SerialAggregator,
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


def test_nocombiner(pango_with_toy_alias):
    expected = {
        "LIFE.1": "ARBITRARYTAXON.0",
        "LIFE.1.2": "ARBITRARYTAXON.1",
        "LIFE.123": "ARBITRARYTAXON.1",
        "LIFE.7": "ARBITRARYTAXON.1",
        "FLYING.8472.1": "ARBITRARYTAXON.0",
        "FLYING.8472.1.2": "ARBITRARYTAXON.0",
        "FLYING.8472.123": "ARBITRARYTAXON.0",
        "FLYING.8472.7": "ARBITRARYTAXON.1",
        "FLYING.123.1": "ARBITRARYTAXON.1",
        "FLYING.123.1.2": "ARBITRARYTAXON.1",
        "FLYING.123.123": "ARBITRARYTAXON.0",
        "FLYING.123.7": "ARBITRARYTAXON.0",
    }

    input_map = {Taxon(k, True): Taxon(v, False) for k, v in expected.items()}
    input_taxa = read_taxa(
        "tests/toy_lineages.txt",
        is_tip=True,
        nomenclature=pango_with_toy_alias,
    )
    agg = ArbitraryAggregator(input_map)

    assert agg.aggregate([taxon for taxon in input_taxa]).to_str() == expected


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
    agg = BasicPhylogeneticAggregator(targets, taxonomy_scheme)

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
    assert agg.aggregate(input_taxa).to_str() == expected


def test_drop_other(pango_with_toy_alias):
    targets = [
        Taxon("FLYING.8472", False),
        Taxon("FLYING.123", False),
    ]

    input_taxa = read_taxa(
        "tests/toy_lineages.txt",
        is_tip=True,
        nomenclature=pango_with_toy_alias,
    )
    tree = pango_with_toy_alias.taxonomy_tree(input_taxa)
    taxonomy_scheme = PhylogeneticTaxonomyScheme(tree)
    agg = BasicPhylogeneticAggregator(
        targets, taxonomy_scheme, off_target="other"
    )

    expected = {
        "LIFE.1": "other",
        "LIFE.1.2": "other",
        "LIFE.123": "other",
        "LIFE.7": "other",
        "FLYING.8472.1": "FLYING.8472",
        "FLYING.8472.1.2": "FLYING.8472",
        "FLYING.8472.123": "FLYING.8472",
        "FLYING.8472.7": "FLYING.8472",
        "FLYING.123.1": "FLYING.123",
        "FLYING.123.1.2": "FLYING.123",
        "FLYING.123.123": "FLYING.123",
        "FLYING.123.7": "FLYING.123",
    }
    assert agg.aggregate(input_taxa).to_str() == expected


def test_drop_self(pango_with_toy_alias):
    targets = [
        Taxon("FLYING.8472", False),
        Taxon("FLYING.123", False),
    ]

    input_taxa = read_taxa(
        "tests/toy_lineages.txt",
        is_tip=True,
        nomenclature=pango_with_toy_alias,
    )
    tree = pango_with_toy_alias.taxonomy_tree(input_taxa)
    taxonomy_scheme = PhylogeneticTaxonomyScheme(tree)
    agg = BasicPhylogeneticAggregator(
        targets, taxonomy_scheme, off_target="self"
    )

    expected = {
        "LIFE.1": "LIFE.1",
        "LIFE.1.2": "LIFE.1.2",
        "LIFE.123": "LIFE.123",
        "LIFE.7": "LIFE.7",
        "FLYING.8472.1": "FLYING.8472",
        "FLYING.8472.1.2": "FLYING.8472",
        "FLYING.8472.123": "FLYING.8472",
        "FLYING.8472.7": "FLYING.8472",
        "FLYING.123.1": "FLYING.123",
        "FLYING.123.1.2": "FLYING.123",
        "FLYING.123.123": "FLYING.123",
        "FLYING.123.7": "FLYING.123",
    }
    assert agg.aggregate(input_taxa).to_str() == expected


def test_overlap_sorted(pango_with_toy_alias):
    targets = [
        Taxon("PYTHONS.0.0.0", False),
        Taxon("LIFE.1", False),
        Taxon("FLYING.8472", False),
        Taxon("FLYING.123", False),
    ]

    input_taxa = read_taxa(
        "tests/toy_lineages.txt",
        is_tip=True,
        nomenclature=pango_with_toy_alias,
    )
    tree = pango_with_toy_alias.taxonomy_tree(input_taxa)
    taxonomy_scheme = PhylogeneticTaxonomyScheme(tree)
    agg = BasicPhylogeneticAggregator(targets, taxonomy_scheme)

    expected = {
        "LIFE.1": "LIFE.1",
        "LIFE.1.2": "LIFE.1",
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
    assert agg.aggregate(input_taxa).to_str() == expected


def test_overlap_unsorted(pango_with_toy_alias):
    targets = [
        Taxon("PYTHONS.0.0.0", False),
        Taxon("LIFE.1", False),
        Taxon("FLYING.8472", False),
        Taxon("FLYING.123", False),
    ]

    input_taxa = read_taxa(
        "tests/toy_lineages.txt",
        is_tip=True,
        nomenclature=pango_with_toy_alias,
    )
    tree = pango_with_toy_alias.taxonomy_tree(input_taxa)
    taxonomy_scheme = PhylogeneticTaxonomyScheme(tree)
    agg = BasicPhylogeneticAggregator(
        targets, taxonomy_scheme, sort_clades=False, warn=False
    )

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
    assert agg.aggregate(input_taxa).to_str() == expected


def test_serial_basic_arbitrary(pango_with_toy_alias):
    basic_targets = [
        Taxon("PYTHONS.0.0.0", False),
        Taxon("LIFE.1", False),
        Taxon("FLYING.8472", False),
        Taxon("FLYING.123", False),
    ]

    arbitrary_dict = {
        Taxon("PYTHONS.0.0.0", False): Taxon("ARBITRARY.1", False),
        Taxon("LIFE.1", False): Taxon("ARBITRARY.2", False),
        Taxon("FLYING.8472", False): Taxon("ARBITRARY.3", False),
        Taxon("FLYING.123", False): Taxon("ARBITRARY.1", False),
    }

    input_taxa = read_taxa(
        "tests/toy_lineages.txt",
        is_tip=True,
        nomenclature=pango_with_toy_alias,
    )
    tree = pango_with_toy_alias.taxonomy_tree(input_taxa)
    taxonomy_scheme = PhylogeneticTaxonomyScheme(tree)
    basic_agg = BasicPhylogeneticAggregator(
        basic_targets, taxonomy_scheme, sort_clades=False, warn=False
    )

    arbitrary_agg = ArbitraryAggregator(arbitrary_dict)

    agg = SerialAggregator([basic_agg, arbitrary_agg])

    expected = {
        "LIFE.1": "ARBITRARY.1",
        "LIFE.1.2": "ARBITRARY.1",
        "LIFE.123": "ARBITRARY.1",
        "LIFE.7": "ARBITRARY.1",
        "FLYING.8472.1": "ARBITRARY.3",
        "FLYING.8472.1.2": "ARBITRARY.3",
        "FLYING.8472.123": "ARBITRARY.3",
        "FLYING.8472.7": "ARBITRARY.3",
        "FLYING.123.1": "ARBITRARY.1",
        "FLYING.123.1.2": "ARBITRARY.1",
        "FLYING.123.123": "ARBITRARY.1",
        "FLYING.123.7": "ARBITRARY.1",
    }
    assert agg.aggregate(input_taxa).to_str() == expected


# def test_historical(pango_with_toy_alias):
#     expected_map = {
#         Taxon("MONTY.42", True): Taxon("MONTY.42", True),
#         Taxon("MONTY.25.25.25", True): Taxon("MONTY.25.25.25", True),
#         Taxon("OF.9.9.9", True): Taxon("PYTHONS.0.0", False),
#         Taxon("FLYING.1.1", True): Taxon("PYTHONS.47.47.47", False),
#         Taxon("FLYING.1.1.1", True): Taxon("PYTHONS.47.47.47", False),
#         Taxon("FLYING.8472", True): Taxon("PYTHONS.47.47.47", False),
#     }

#     observed_taxa_str = [taxon.name for taxon in expected_map]


#     tree = pango_with_toy_alias.taxonomy_tree(observed_taxa_str)
#     taxonomy_scheme = PhylogeneticTaxonomyScheme(tree)

#     arbitrary_unused_date = "1000-1-1"
#     aggregator = HistoricalAggregator(
#         taxonomy_scheme, pango_with_toy_alias, as_of=arbitrary_unused_date
#     )

#     aggregation = aggregator.aggregate(expected_map.keys())

#     assert aggregation == expected_map
