import pytest

from cladecombiner import PangoSc2Nomenclature


def test_split():
    pn = PangoSc2Nomenclature()
    assert pn.split("CIRCUS.1001.1002.1003") == [
        "CIRCUS",
        "1001",
        "1002",
        "1003",
    ]


def test_join():
    pn = PangoSc2Nomenclature()
    assert (
        pn.join(["CIRCUS", "1001", "1002", "1003"]) == "CIRCUS.1001.1002.1003"
    )


def test_special_norecomb(pango_with_toy_alias):
    assert pango_with_toy_alias.is_special("MONTY")
    assert not pango_with_toy_alias.is_special("PYTHONS")
    assert pango_with_toy_alias.is_special("XPYTHONS")


def test_special_recomb(pango_with_recomb_alias):
    assert pango_with_recomb_alias.special == ["A"]
    assert pango_with_recomb_alias.is_special("A")
    assert not pango_with_recomb_alias.is_special("B")
    assert pango_with_recomb_alias.is_special("XE")


def test_recomb(pango_with_toy_alias):
    assert pango_with_toy_alias.is_hybrid("XPYTHONS")
    assert not pango_with_toy_alias.is_hybrid("PYTHONS")


def test_root(pango_with_toy_alias):
    assert pango_with_toy_alias.is_root("")
    assert not pango_with_toy_alias.is_root("PYTHONS")


def test_ambiguous(pango_with_toy_alias):
    assert pango_with_toy_alias.is_ambiguous("PYTHONS.1.2*")
    assert not pango_with_toy_alias.is_ambiguous("PYTHONS.1.2")


def test_valid(pango_with_toy_alias):
    assert pango_with_toy_alias.is_valid_alias("PYTHONS")
    assert not pango_with_toy_alias.is_valid_name("PYTHONS")
    assert pango_with_toy_alias.is_valid_name("PYTHONS", min_sublevels=0)
    assert not pango_with_toy_alias.is_valid_name("PYTHONS.1.1.2.3.4.5.2")
    assert pango_with_toy_alias.is_valid_name(
        "PYTHONS.1.1.2.3.4.5.2", max_sublevels=10000
    )
    pango_with_toy_alias.max_sublevels = 10000
    assert pango_with_toy_alias.is_valid_name("PYTHONS.1.1.2.3.4.5.2")
    pango_with_toy_alias.max_sublevels = 3


def test_longer_name_nolist():
    pn = PangoSc2Nomenclature()
    with pytest.raises(RuntimeError):
        pn.longer_name("CIRCUS.1001.1002.1003")


def test_longer_name(pango_with_toy_alias):
    observed = pango_with_toy_alias.longer_name("CIRCUS.1001.1002.1003")
    expected = "MONTY.42.42.42.47.47.47.8472.8472.8472.1001.1002.1003"
    assert observed == expected


def test_alias_bottom(pango_with_toy_alias):
    observed = pango_with_toy_alias.longer_name("MONTY.1001.1002.1003")
    expected = "MONTY.1001.1002.1003"
    assert observed == expected


def test_shorter_name(pango_with_toy_alias):
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

    observed = pango_with_toy_alias.shorter_name("MONTY.42.42.42.47.47.47.1.1")
    expected = "FLYING.1.1"
    assert observed == expected


def test_history_norecomb(pango_with_toy_alias):
    observed = pango_with_toy_alias.get_history(
        "CARNIVAL.1001.1002.1003", stop_at_hybrid=True
    )
    expected = [
        "",
        "MONTY",
        "MONTY.42",
        "MONTY.42.42",
        "MONTY.42.42.42",
        "MONTY.42.42.42.47",
        "MONTY.42.42.42.47.47",
        "MONTY.42.42.42.47.47.47",
        "MONTY.42.42.42.47.47.47.1",
        "MONTY.42.42.42.47.47.47.1.1",
        "MONTY.42.42.42.47.47.47.1.1.1",
        "MONTY.42.42.42.47.47.47.1.1.1.1001",
        "MONTY.42.42.42.47.47.47.1.1.1.1001.1002",
        "MONTY.42.42.42.47.47.47.1.1.1.1001.1002.1003",
    ]
    for e, o in zip(expected[1:], observed[1:]):
        assert pango_with_toy_alias.equals_ignore_alias(e, o)
    assert observed[0] == ""
    for o in observed:
        assert o == pango_with_toy_alias.longer_name(o)

    observed = pango_with_toy_alias.get_history(
        "CIRCUS.1001.1002.1003", stop_at_hybrid=True
    )
    expected = [
        "",
        "MONTY",
        "MONTY.42",
        "MONTY.42.42",
        "MONTY.42.42.42",
        "MONTY.42.42.42.47",
        "MONTY.42.42.42.47.47",
        "MONTY.42.42.42.47.47.47",
        "MONTY.42.42.42.47.47.47.8472",
        "MONTY.42.42.42.47.47.47.8472.8472",
        "MONTY.42.42.42.47.47.47.8472.8472.8472",
        "MONTY.42.42.42.47.47.47.8472.8472.8472.1001",
        "MONTY.42.42.42.47.47.47.8472.8472.8472.1001.1002",
        "MONTY.42.42.42.47.47.47.8472.8472.8472.1001.1002.1003",
    ]
    for e, o in zip(expected[1:], observed[1:]):
        assert pango_with_toy_alias.equals_ignore_alias(e, o)
    assert observed[0] == ""
    for o in observed:
        assert o == pango_with_toy_alias.longer_name(o)


def test_history_recomb(pango_with_recomb_alias):
    observed = pango_with_recomb_alias.get_history("XB", stop_at_hybrid=True)
    assert observed == ["XB"]

    observed = pango_with_recomb_alias.get_history("XB", stop_at_hybrid=False)
    expected = [
        "",
        "A",
        "A.2",
        "A.2.1",
        "A.2.1.2",
        "A.2.1.2.2",
        "XA",
        "XA.1",
        "XA.1.1",
        "XA.1.1.1",
        "XB",
    ]
    for e, o in zip(expected[1:], observed[1:]):
        print("Asserting " + o + " == " + e)
        assert pango_with_recomb_alias.equals_ignore_alias(e, o)
    assert observed[0] == ""
