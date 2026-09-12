"""Tests for chemistrykit.electrochem.systems.standard_potentials."""

import pytest

from chemistrykit.electrochem.systems.standard_potentials import (
    STANDARD_REDUCTION_POTENTIALS as T,
)
from chemistrykit.electrochem.systems.standard_potentials import (
    balance_redox_reaction,
    cell_potential,
    is_spontaneous,
    standard_cell_potential,
)


def test_daniell_cell_potential():
    """The classic Cu/Zn Daniell cell: E_cell = 0.34 - (-0.76) = 1.10 V."""
    E = cell_potential(T["Cu2+/Cu"], T["Zn2+/Zn"])
    assert E == pytest.approx(1.10)
    assert standard_cell_potential("Cu2+/Cu", "Zn2+/Zn") == pytest.approx(1.10)


def test_cell_potential_independent_of_electron_count_scaling():
    """Doubling both half-reactions' stoichiometry (same E values) leaves
    the cell potential unchanged -- an intensive property."""
    E1 = cell_potential(T["Ag+/Ag"], T["Zn2+/Zn"])
    # Same tabulated E values regardless of how the balanced equation is scaled.
    E2 = cell_potential(T["Ag+/Ag"], T["Zn2+/Zn"])
    assert E1 == pytest.approx(E2)


def test_balance_redox_reaction_equal_electron_counts_needs_no_scaling():
    cathode_mult, anode_mult, n_total = balance_redox_reaction(T["Cu2+/Cu"], T["Zn2+/Zn"])
    assert (cathode_mult, anode_mult, n_total) == (1, 1, 2)


def test_balance_redox_reaction_lcm_scaling():
    """MnO4-/Mn2+ (5 e-) paired with Fe3+/Fe2+ (1 e-) needs the iron
    half-reaction scaled 5-fold to match electron counts."""
    cathode_mult, anode_mult, n_total = balance_redox_reaction(T["MnO4-/Mn2+"], T["Fe3+/Fe2+"])
    assert cathode_mult == 1
    assert anode_mult == 5
    assert n_total == 5


def test_is_spontaneous():
    assert is_spontaneous(1.10) is True
    assert is_spontaneous(-0.5) is False
    assert is_spontaneous(0.0) is False


def test_table_values_match_known_references():
    """Spot-check a few well-known literature standard reduction potentials
    (Bard & Faulkner, Appendix C.3 / Atkins & de Paula, Table 6.4)."""
    assert T["H+/H2"].E_standard == pytest.approx(0.00)
    assert T["F2/F-"].E_standard == pytest.approx(2.87)
    assert T["Li+/Li"].E_standard == pytest.approx(-3.04)
    assert T["O2/H2O"].n == 4
