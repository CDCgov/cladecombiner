import pytest

import cladecombiner


@pytest.fixture(scope="session")
def pango_with_toy_alias():
    r"""
    The alias file toy_alias.json is based on this tree, which follows the
        Pango nomenclatural rules, with 3 sublevels, but is not a real tree.
    The aliases are not all complete, in that to get the fully long name of
        some names, you must use the alias key repeatedly.
    The only "special" lineage descending directly from the root is MONTY.

      /---MONTY.25---MONTY.25.25---MONTY.25.25.25=PATHG
      |
      |
    MONTY                                            /---PYTHONS.0---PYTHONS.0.0---PYTHONS.0.0.0=LIFE---LIFE.7---LIFE.7.7---LIFE.7.7.7=OF---OF.9---OF.9.9---OF.9.9.9=BRIAN
      |                                              |
      |                                              |
      \---MONTY.42---MONTY.42.42---MONTY.42.42.42=PYTHONS                                         /---FLYING.1---FLYING.1.1---FLYING.1.1.1=CARNIVAL
                                                     |                                            |
                                                     |                                            |
                                                     \---PYTHONS.47---PYTHONS.47.47---PYTHONS.47.47.47=FLYING
                                                                                                  |
                                                                                                  |
                                                                                                  \---FLYING.8472---FLYING.8472.8472---FLYING.8472.8472.8472=CIRCUS
    """
    pango = cladecombiner.PangoSc2Nomenclature()
    pango.special = ["MONTY"]
    pango.setup_alias_map(fp="tests/toy_alias.json")
    return pango


@pytest.fixture(scope="session")
def pango_with_recomb_alias():
    r"""
    The alias file recomb_alias.json is based on this tree, which follows the
        Pango nomenclatural rules, with 3 sublevels, but is not a real tree.
    Internal bifurcating nodes are labeled directly, where they map to an aliased
        shorter name, this is noted as <long>=<short>.
    Internal multifurcating nodes are labeled once and other descendants are marked
        with a +.
    Recombination events are denoted with &, which is followed by a branch and then
        the label of the recombinant taxon.
    The only "special" lineage descending directly from the root is A.

           /-----A.1.1
           |
    /-----A.1
    |      |
    |      \-----A.1.2
    |
    |              /-----A.2.1.1
    |              |
    A      /-----A.2.1       /-----B.1                      /-----XA.1.1.1-----\
    |      |       |         |                              |                  |
    \-----A.2      \-----A.2.1.2=B                 /-----XA.1.1                &-----XB
           |                 |                     |        |                  |
           |                 \-----B.2     /-----XA.1       \-----XA.1.1.2-----/
           |                        |      |       |
           |                        |      |       \-----XA.1.2
           |                        |      |                                  /-----XC
           \-----A.2.2--------------&-----XA                                  |
                                           |       /-----XA.2.1               +-----XD
                                           |       |                          |
                                           \-----XA.2       /-----XA.2.2.1----&-----XE
                                                   |        |                 |
                                                   \-----XA.2.2        /-----C.1
                                                            |          |
                                                            \-----XA.2.2.2=C
                                                                       |
                                                                       +-----C.2-----\
                                                                       |             |
                                                                       |             |
                                                                       +-----C.3-----&-----XF
                                                                       |             |
                                                                       |             |
                                                                       \-----C.4-----/
    """

    pango = cladecombiner.PangoSc2Nomenclature()
    pango.special = ["A"]
    pango.setup_alias_map(fp="tests/recomb_alias.json")
    return pango
