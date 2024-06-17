import pytest

from cladecombiner import PangoSc2Nomenclature


@pytest.fixture
def pango_with_toy_alias():
    pango = PangoSc2Nomenclature()
    # Toy alias is set up to look a lot like real Pango,
    # that includes "MONTY" which dealiases to the root ("")
    # which means it is a special taxon
    pango.special = ["MONTY"]
    pango.setup_alias_map(fp="tests/toy_alias.json")
    return pango


class Test_pango_conversions:
    def test_split(self):
        pn = PangoSc2Nomenclature()
        assert pn.split("CIRCUS.1001.1002.1003") == [
            "CIRCUS",
            "1001",
            "1002",
            "1003",
        ]

    def test_join(self):
        pn = PangoSc2Nomenclature()
        assert (
            pn.join(["CIRCUS", "1001", "1002", "1003"])
            == "CIRCUS.1001.1002.1003"
        )


@pytest.mark.usefixtures("pango_with_toy_alias")
class Test_pango_name_types:
    def test_special(self, pango_with_toy_alias):
        assert pango_with_toy_alias.is_special("MONTY")
        assert not pango_with_toy_alias.is_special("PYTHONS")
        assert pango_with_toy_alias.is_special("XPYTHONS")

    def test_recomb(self, pango_with_toy_alias):
        assert pango_with_toy_alias.is_hybrid("XPYTHONS")
        assert not pango_with_toy_alias.is_hybrid("PYTHONS")

    def test_root(self, pango_with_toy_alias):
        assert pango_with_toy_alias.is_root("")
        assert not pango_with_toy_alias.is_root("PYTHONS")

    def test_amb(self, pango_with_toy_alias):
        assert pango_with_toy_alias.is_ambiguous("PYTHONS.1.2*")
        assert not pango_with_toy_alias.is_ambiguous("PYTHONS.1.2")

    def test_valid(self, pango_with_toy_alias):
        assert pango_with_toy_alias.is_valid_alias("PYTHONS")
        assert not pango_with_toy_alias.is_valid_name("PYTHONS")
        assert pango_with_toy_alias.is_valid_name("PYTHONS", min_sublevels=0)
        assert not pango_with_toy_alias.is_valid_name("PYTHONS.1.1.2.3.4.5.2")
        assert pango_with_toy_alias.is_valid_name(
            "PYTHONS.1.1.2.3.4.5.2", max_sublevels=10000
        )
        pango_with_toy_alias.max_sublevels = 10000
        assert pango_with_toy_alias.is_valid_name("PYTHONS.1.1.2.3.4.5.2")


@pytest.mark.usefixtures("pango_with_toy_alias")
class Test_pango_alias:
    def test_longer_name_nolist(self):
        pn = PangoSc2Nomenclature()
        with pytest.raises(RuntimeError):
            pn.longer_name("CIRCUS.1001.1002.1003")

    def test_alias(self, pango_with_toy_alias):
        observed = pango_with_toy_alias.longer_name("CIRCUS.1001.1002.1003")
        expected = "MONTY.42.42.42.47.47.47.8472.8472.8472.1001.1002.1003"
        assert observed == expected

    def test_alias_bottom(self, pango_with_toy_alias):
        observed = pango_with_toy_alias.longer_name("MONTY.1001.1002.1003")
        expected = "MONTY.1001.1002.1003"
        assert observed == expected

    def test_shorter_name(self, pango_with_toy_alias):
        observed = pango_with_toy_alias.shorter_name(
            "MONTY.42.42.42.47.47.47.8472.8472.8472.1001.1002.1003",
        )
        expected = "CIRCUS.1001.1002.1003"
        assert observed == expected

        observed = pango_with_toy_alias.shorter_name(
            "MONTY.42.42.42.47.47.47.1.1.1.1001.1002.1003"
        )
        expected = "CARNIVAL.1001.1002.1003"
        assert observed == expected

        observed = pango_with_toy_alias.shorter_name(
            "MONTY.42.42.42.47.47.47.1.1"
        )
        expected = "FLYING.1.1"
        assert observed == expected
