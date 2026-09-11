"""Tests for chemistrykit.statmech.systems.partition_functions against closed-form results."""

import numpy as np
import pytest
import scipy.constants as sc

from chemistrykit.constants import K_B, R
from chemistrykit.statmech.systems.partition_functions import (
    IdealGasMolecule,
    RotationalPartitionFunctionLinear,
    TranslationalPartitionFunction,
    VibrationalPartitionFunctionHarmonic,
    sackur_tetrode_entropy,
)

ARGON_MASS = 39.948 * sc.atomic_mass


def test_translational_internal_energy_is_equipartition_value():
    q = TranslationalPartitionFunction(mass=ARGON_MASS, volume=1.0e-3)
    U = q.internal_energy(300.0, N=1.0)
    assert U == pytest.approx(1.5 * K_B * 300.0)


def test_translational_heat_capacity_is_three_halves_nkb():
    q = TranslationalPartitionFunction(mass=ARGON_MASS, volume=1.0e-3)
    assert q.heat_capacity_v(300.0, N=1.0) == pytest.approx(1.5 * K_B)


def test_translational_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        TranslationalPartitionFunction(mass=0.0, volume=1.0)
    with pytest.raises(ValueError):
        TranslationalPartitionFunction(mass=1.0, volume=-1.0)


def test_sackur_tetrode_matches_tabulated_argon_entropy():
    """Textbook value: S_m(Ar, 298.15 K, 1 bar) ~ 154.8 J/(mol K)."""
    S = sackur_tetrode_entropy(mass=ARGON_MASS, T=298.15, P=1.0e5)
    assert S == pytest.approx(154.8, abs=0.5)


def test_sackur_tetrode_entropy_increases_with_temperature():
    S_low = sackur_tetrode_entropy(mass=ARGON_MASS, T=200.0, P=1.0e5)
    S_high = sackur_tetrode_entropy(mass=ARGON_MASS, T=400.0, P=1.0e5)
    assert S_high > S_low


def test_sackur_tetrode_entropy_decreases_with_pressure():
    S_low_p = sackur_tetrode_entropy(mass=ARGON_MASS, T=298.15, P=1.0e4)
    S_high_p = sackur_tetrode_entropy(mass=ARGON_MASS, T=298.15, P=1.0e6)
    assert S_low_p > S_high_p


def test_rotational_classical_limit_heat_capacity_is_nkb():
    q = RotationalPartitionFunctionLinear(moment_of_inertia=1.45e-46, symmetry_number=1)
    assert q.heat_capacity_v(300.0, N=1.0) == pytest.approx(K_B)


def test_rotational_symmetry_number_halves_partition_function():
    I = 1.45e-46
    q1 = RotationalPartitionFunctionLinear(moment_of_inertia=I, symmetry_number=1)
    q2 = RotationalPartitionFunctionLinear(moment_of_inertia=I, symmetry_number=2)
    assert q2.value(300.0) == pytest.approx(q1.value(300.0) / 2.0)


def test_rotational_rejects_invalid_parameters():
    with pytest.raises(ValueError):
        RotationalPartitionFunctionLinear(moment_of_inertia=-1.0)
    with pytest.raises(ValueError):
        RotationalPartitionFunctionLinear(moment_of_inertia=1.0, symmetry_number=0)


def test_vibrational_heat_capacity_vanishes_at_low_temperature():
    q = VibrationalPartitionFunctionHarmonic(frequency=8.7e13)
    assert q.heat_capacity_v(1.0, N=1.0) == pytest.approx(0.0, abs=1e-40)


def test_vibrational_heat_capacity_approaches_nkb_at_high_temperature():
    q = VibrationalPartitionFunctionHarmonic(frequency=8.7e13)
    ratio = q.heat_capacity_v(1.0e7, N=1.0) / K_B
    assert ratio == pytest.approx(1.0, abs=1e-3)


def test_vibrational_from_wavenumber():
    q = VibrationalPartitionFunctionHarmonic.from_wavenumber(2886.0)
    expected_frequency = sc.c * 2886.0 * 100.0
    assert q.frequency == pytest.approx(expected_frequency)


def test_vibrational_rejects_nonpositive_frequency():
    with pytest.raises(ValueError):
        VibrationalPartitionFunctionHarmonic(frequency=0.0)


def test_ideal_gas_molecule_high_temperature_heat_capacity_matches_equipartition():
    molecule = IdealGasMolecule(mass=ARGON_MASS, volume=1.0e-3, moment_of_inertia=1.45e-46, vibrational_frequencies=[8.7e13])
    Cv_over_R = molecule.heat_capacity_v(1.0e7) / R
    assert Cv_over_R == pytest.approx(3.5, abs=1e-3)


def test_ideal_gas_molecule_atom_has_no_rotational_contribution():
    molecule = IdealGasMolecule(mass=ARGON_MASS, volume=1.0e-3)
    Cv_over_R = molecule.heat_capacity_v(300.0) / R
    assert Cv_over_R == pytest.approx(1.5)


def test_ideal_gas_molecule_thermodynamic_functions_bundle():
    molecule = IdealGasMolecule(mass=ARGON_MASS, volume=1.0e-3)
    functions = molecule.thermodynamic_functions(300.0, N=1.0)
    assert functions.T == 300.0
    assert functions.U == pytest.approx(molecule.internal_energy(300.0, N=1.0))
    assert functions.A == pytest.approx(functions.U - functions.T * functions.S)


def test_helmholtz_free_energy_consistent_with_u_and_s():
    q = TranslationalPartitionFunction(mass=ARGON_MASS, volume=1.0e-3)
    T = 300.0
    A = q.helmholtz_free_energy(T, N=1.0)
    assert A == pytest.approx(q.internal_energy(T, N=1.0) - T * q.entropy(T, N=1.0))


def test_translational_partition_function_scales_with_volume():
    q1 = TranslationalPartitionFunction(mass=ARGON_MASS, volume=1.0e-3)
    q2 = TranslationalPartitionFunction(mass=ARGON_MASS, volume=2.0e-3)
    assert q2.value(300.0) == pytest.approx(2.0 * q1.value(300.0))


def test_all_partition_functions_are_positive():
    T = np.linspace(50.0, 500.0, 5)
    trans = TranslationalPartitionFunction(mass=ARGON_MASS, volume=1e-3)
    rot = RotationalPartitionFunctionLinear(moment_of_inertia=1.45e-46)
    vib = VibrationalPartitionFunctionHarmonic(frequency=8.7e13)
    for t in T:
        assert trans.value(t) > 0
        assert rot.value(t) > 0
        assert vib.value(t) > 0
