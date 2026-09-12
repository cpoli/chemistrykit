"""Tests for chemistrykit.structure.systems.lewis against known formal charges and oxidation states."""

import pytest

from chemistrykit.structure.systems.lewis import LewisStructure


def test_water_formal_charges_all_zero():
    water = LewisStructure(symbols=["O", "H", "H"], bond_orders={(0, 1): 1, (0, 2): 1}, lone_pairs={0: 2})
    charges = water.formal_charges()
    assert charges == pytest.approx({0: 0.0, 1: 0.0, 2: 0.0})
    assert water.total_formal_charge() == pytest.approx(0.0)


def test_water_oxidation_states():
    water = LewisStructure(symbols=["O", "H", "H"], bond_orders={(0, 1): 1, (0, 2): 1}, lone_pairs={0: 2})
    states = water.oxidation_states()
    assert states[0] == pytest.approx(-2.0)
    assert states[1] == pytest.approx(1.0)
    assert states[2] == pytest.approx(1.0)


def test_ammonium_nitrogen_has_plus_one_formal_charge():
    nh4_plus = LewisStructure(symbols=["N", "H", "H", "H", "H"], bond_orders={(0, k): 1 for k in (1, 2, 3, 4)})
    charges = nh4_plus.formal_charges()
    assert charges[0] == pytest.approx(1.0)
    # Consistency check: total formal charge equals the actual net charge of NH4+.
    assert nh4_plus.total_formal_charge() == pytest.approx(1.0)


def test_hydronium_oxygen_has_plus_one_formal_charge():
    h3o_plus = LewisStructure(symbols=["O", "H", "H", "H"], bond_orders={(0, k): 1 for k in (1, 2, 3)}, lone_pairs={0: 1})
    charges = h3o_plus.formal_charges()
    assert charges[0] == pytest.approx(1.0)
    assert h3o_plus.total_formal_charge() == pytest.approx(1.0)


def test_hydrogen_peroxide_oxidation_state_is_minus_one():
    h2o2 = LewisStructure(symbols=["H", "O", "O", "H"], bond_orders={(0, 1): 1, (1, 2): 1, (2, 3): 1}, lone_pairs={1: 2, 2: 2})
    states = h2o2.oxidation_states()
    assert states[1] == pytest.approx(-1.0)
    assert states[2] == pytest.approx(-1.0)
    assert states[0] == pytest.approx(1.0)
    assert states[3] == pytest.approx(1.0)


def test_carbon_dioxide_carbon_oxidation_state_is_plus_four():
    co2 = LewisStructure(symbols=["O", "C", "O"], bond_orders={(0, 1): 2, (1, 2): 2}, lone_pairs={0: 2, 2: 2})
    states = co2.oxidation_states()
    assert states[1] == pytest.approx(4.0)
    assert states[0] == pytest.approx(-2.0)
    assert states[2] == pytest.approx(-2.0)


def test_rejects_bond_index_out_of_range():
    with pytest.raises(ValueError):
        LewisStructure(symbols=["H", "H"], bond_orders={(0, 5): 1})


def test_rejects_lone_pair_index_out_of_range():
    with pytest.raises(ValueError):
        LewisStructure(symbols=["H", "H"], bond_orders={(0, 1): 1}, lone_pairs={5: 1})
