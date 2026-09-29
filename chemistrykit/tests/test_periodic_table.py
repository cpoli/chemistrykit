"""Tests for chemistrykit.periodic_table: element coverage and formula parsing."""

import pytest

from chemistrykit.periodic_table import PERIODIC_TABLE, electronegativity, formula_charge, get_element, molar_mass, parse_formula, valence_electrons


def test_table_covers_every_element_through_oganesson():
    assert len(PERIODIC_TABLE) == 118
    assert sorted(e.atomic_number for e in PERIODIC_TABLE.values()) == list(range(1, 119))


@pytest.mark.parametrize(
    "symbol, z, mass", [("Au", 79, 196.97), ("Pt", 78, 195.08), ("Hg", 80, 200.59), ("Pb", 82, 207.2), ("U", 92, 238.03), ("Gd", 64, 157.25)]
)
def test_heavy_elements_are_present(symbol, z, mass):
    element = get_element(symbol)
    assert element.atomic_number == z
    assert element.atomic_mass == pytest.approx(mass)
    assert get_element(z) is element


def test_heavy_main_group_lookups():
    assert valence_electrons("Pb") == 4
    assert electronegativity("Au") == pytest.approx(2.54)


@pytest.mark.parametrize(
    "formula, expected",
    [
        ("H2O", {"H": 2, "O": 1}),
        ("Ca(OH)2", {"Ca": 1, "O": 2, "H": 2}),
        ("Al2(SO4)3", {"Al": 2, "S": 3, "O": 12}),
        ("K4[Fe(CN)6]", {"K": 4, "Fe": 1, "C": 6, "N": 6}),
        ("CuSO4*5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
        ("CuSO4·5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
        ("Fe0.95O", {"Fe": 0.95, "O": 1}),
        ("e-", {}),
    ],
)
def test_parse_formula(formula, expected):
    assert parse_formula(formula) == pytest.approx(expected)


@pytest.mark.parametrize("formula, charge", [("Fe3+", 3), ("O2-", -2), ("MnO4-", -1), ("NH4+", 1), ("SO4^2-", -2), ("Hg2^2+", 2), ("O2^-", -1), ("H2O", 0)])
def test_formula_charge(formula, charge):
    assert formula_charge(formula) == charge


def test_polyatomic_ion_count_is_not_read_as_charge():
    assert parse_formula("MnO4-") == {"Mn": 1, "O": 4}
    assert parse_formula("Hg2^2+") == {"Hg": 2}


@pytest.mark.parametrize("bad", ["Ca(OH", "H2)", "Xx", "2", "", "H2O$"])
def test_malformed_formulas_raise(bad):
    with pytest.raises(ValueError):
        parse_formula(bad)


def test_molar_mass_string_matches_dict():
    assert molar_mass("Ca(OH)2") == pytest.approx(molar_mass({"Ca": 1, "O": 2, "H": 2}))
    assert molar_mass("Ca(OH)2") == pytest.approx(40.078 + 2 * (15.999 + 1.008))
