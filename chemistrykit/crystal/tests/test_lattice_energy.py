"""Tests for chemistrykit.crystal.systems.lattice_energy against known NaCl lattice-energy estimates."""

import pytest

from chemistrykit.crystal.systems.lattice_energy import BornLandeLatticeEnergy, KapustinskiiLatticeEnergy


def test_born_lande_nacl_matches_hand_calculation():
    model = BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=282e-12, born_exponent=8.0)
    assert model.lattice_energy() / 1000.0 == pytest.approx(-753.4, abs=0.1)


def test_born_lande_lattice_energy_is_negative_and_magnitude_increases_with_madelung_constant():
    low_M = BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=282e-12, born_exponent=8.0, madelung_constant=1.0)
    high_M = BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=282e-12, born_exponent=8.0, madelung_constant=2.0)
    assert low_M.lattice_energy() < 0.0
    assert abs(high_M.lattice_energy()) > abs(low_M.lattice_energy())


def test_born_lande_magnitude_increases_with_charge():
    plus1 = BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=282e-12, born_exponent=8.0)
    plus2 = BornLandeLatticeEnergy(cation_charge=2, anion_charge=2, r0=282e-12, born_exponent=8.0)
    assert abs(plus2.lattice_energy()) > abs(plus1.lattice_energy())


def test_born_lande_rejects_invalid_parameters():
    with pytest.raises(ValueError):
        BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=-1.0, born_exponent=8.0)
    with pytest.raises(ValueError):
        BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=1.0, born_exponent=1.0)


def test_kapustinskii_nacl_matches_hand_calculation():
    model = KapustinskiiLatticeEnergy(n_ions=2, cation_charge=1, anion_charge=1, r_cation_pm=102.0, r_anion_pm=181.0)
    assert model.lattice_energy() / 1000.0 == pytest.approx(-746.2, abs=0.1)


def test_kapustinskii_agrees_with_born_lande_within_10_percent_for_nacl():
    born_lande = BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=282e-12, born_exponent=8.0).lattice_energy()
    kapustinskii = KapustinskiiLatticeEnergy(n_ions=2, cation_charge=1, anion_charge=1, r_cation_pm=102.0, r_anion_pm=181.0).lattice_energy()
    assert abs(born_lande - kapustinskii) / abs(born_lande) < 0.10


def test_kapustinskii_rejects_invalid_parameters():
    with pytest.raises(ValueError):
        KapustinskiiLatticeEnergy(n_ions=2, cation_charge=1, anion_charge=1, r_cation_pm=-1.0, r_anion_pm=181.0)
    with pytest.raises(ValueError):
        KapustinskiiLatticeEnergy(n_ions=1, cation_charge=1, anion_charge=1, r_cation_pm=100.0, r_anion_pm=100.0)


def test_kapustinskii_magnitude_increases_with_number_of_ions():
    two_ions = KapustinskiiLatticeEnergy(n_ions=2, cation_charge=1, anion_charge=1, r_cation_pm=100.0, r_anion_pm=180.0)
    three_ions = KapustinskiiLatticeEnergy(n_ions=3, cation_charge=1, anion_charge=1, r_cation_pm=100.0, r_anion_pm=180.0)
    assert abs(three_ions.lattice_energy()) > abs(two_ions.lattice_energy())
