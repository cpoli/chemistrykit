"""Tests for chemistrykit.stoichiometry: balancing conserves atoms and charge; yield bookkeeping."""

import pytest

from chemistrykit.periodic_table import formula_charge, parse_formula
from chemistrykit.stoichiometry import (
    BalancedEquation,
    balance_equation,
    empirical_formula,
    limiting_reagent,
    mass_percent,
    percent_yield,
    reaction_extent,
    theoretical_yield,
)


def _conserved(eq: BalancedEquation) -> bool:
    totals: dict = {}
    for species, nu in eq.coefficients.items():
        for element, count in parse_formula(species).items():
            totals[element] = totals.get(element, 0) + nu * count
        totals["charge"] = totals.get("charge", 0) + nu * formula_charge(species)
    return all(v == 0 for v in totals.values())


@pytest.mark.parametrize(
    "equation, reactants, products",
    [
        ("Fe + O2 -> Fe2O3", {"Fe": 4, "O2": 3}, {"Fe2O3": 2}),
        ("C3H8 + O2 -> CO2 + H2O", {"C3H8": 1, "O2": 5}, {"CO2": 3, "H2O": 4}),
        ("KMnO4 + HCl = KCl + MnCl2 + H2O + Cl2", {"KMnO4": 2, "HCl": 16}, {"KCl": 2, "MnCl2": 2, "H2O": 8, "Cl2": 5}),
        ("Ca(OH)2 + H3PO4 -> Ca3(PO4)2 + H2O", {"Ca(OH)2": 3, "H3PO4": 2}, {"Ca3(PO4)2": 1, "H2O": 6}),
        ("MnO4- + Fe2+ + H+ -> Mn2+ + Fe3+ + H2O", {"MnO4-": 1, "Fe2+": 5, "H+": 8}, {"Mn2+": 1, "Fe3+": 5, "H2O": 4}),
        ("Cr2O7^2- + H+ + e- -> Cr3+ + H2O", {"Cr2O7^2-": 1, "H+": 14, "e-": 6}, {"Cr3+": 2, "H2O": 7}),
    ],
)
def test_balance_equation_known_results(equation, reactants, products):
    eq = balance_equation(equation)
    assert eq.reactants == reactants
    assert eq.products == products
    assert _conserved(eq)


def test_existing_coefficients_are_recomputed():
    assert balance_equation("5H2 + 7O2 -> 3H2O") == balance_equation("H2 + O2 -> H2O")


@pytest.mark.parametrize("equation", ["H2 -> O2", "H2 + O2 -> H2O + H2O2", "H2O -> H2 + O2 + H2O2", "H2 + O2", "H2O + O2 -> H2"])
def test_unbalanceable_or_ambiguous_equations_raise(equation):
    with pytest.raises(ValueError):
        balance_equation(equation)


def test_limiting_reagent_and_yield():
    moles = {"N2": 1.0, "H2": 2.0}
    eq = "N2 + H2 -> NH3"
    assert limiting_reagent(eq, moles) == "H2"
    assert reaction_extent(eq, moles) == pytest.approx(2.0 / 3.0)
    assert theoretical_yield(eq, moles, "NH3") == pytest.approx(4.0 / 3.0)
    assert theoretical_yield(eq, moles, "NH3", in_grams=True) == pytest.approx(4.0 / 3.0 * (14.007 + 3 * 1.008))
    assert percent_yield(1.0, 4.0 / 3.0) == pytest.approx(75.0)


def test_yield_errors():
    with pytest.raises(ValueError):
        theoretical_yield("N2 + H2 -> NH3", {"N2": 1.0}, "N2")
    with pytest.raises(ValueError):
        limiting_reagent("N2 + H2 -> NH3", {"O2": 1.0})
    with pytest.raises(ValueError):
        percent_yield(1.0, 0.0)


def test_mass_percent_sums_to_100():
    assert sum(mass_percent("Ca3(PO4)2").values()) == pytest.approx(100.0)


@pytest.mark.parametrize(
    "formula, expected", [("C6H12O6", {"C": 1, "H": 2, "O": 1}), ("Fe3O4", {"Fe": 3, "O": 4}), ("C4H10", {"C": 2, "H": 5}), ("P4O10", {"P": 2, "O": 5})]
)
def test_empirical_formula_inverts_mass_percent(formula, expected):
    assert empirical_formula(mass_percent(formula)) == expected
