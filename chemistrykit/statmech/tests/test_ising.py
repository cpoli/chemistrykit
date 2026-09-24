"""Tests for chemistrykit.statmech.systems.ising against exact enumeration and closed forms."""

import itertools

import numpy as np
import pytest

from chemistrykit.constants import K_B
from chemistrykit.statmech.systems.ising import Ising1D, Ising2DOnsager, ising_2d_critical_temperature, kramers_wannier_dual_coupling

J = 100.0 * K_B


def brute_force_ring(J, h, T, N):
    beta = 1.0 / (K_B * T)
    Z = 0.0
    for spins in itertools.product((-1, 1), repeat=N):
        s = np.array(spins)
        E = -J * np.sum(s * np.roll(s, 1)) - h * np.sum(s)
        Z += np.exp(-beta * E)
    return Z


@pytest.mark.parametrize("h", [0.0, 30.0 * K_B])
def test_transfer_matrix_matches_brute_force_enumeration(h):
    chain = Ising1D(coupling=J, field=h)
    assert chain.partition_function(150.0, 8) == pytest.approx(brute_force_ring(J, h, 150.0, 8), rel=1e-10)


def test_1d_magnetization_is_derivative_of_free_energy():
    h, dh, T = 20.0 * K_B, 1e-4 * K_B, 120.0
    m = Ising1D(J, h).magnetization(T)
    numeric = -(Ising1D(J, h + dh).free_energy_per_spin(T) - Ising1D(J, h - dh).free_energy_per_spin(T)) / (2 * dh)
    assert m == pytest.approx(numeric, rel=1e-6)


def test_1d_has_no_spontaneous_magnetization_and_matches_heat_capacity_formula():
    chain = Ising1D(J)
    assert float(chain.magnetization(10.0)) == 0.0
    T, dT = 80.0, 1e-3
    f = chain.free_energy_per_spin
    # c = -T d^2f/dT^2
    numeric = -T * (f(T + dT) - 2 * f(T) + f(T - dT)) / dT**2
    assert numeric == pytest.approx(chain.heat_capacity_per_spin(T), rel=1e-4)


def test_kramers_wannier_duality_relations():
    K = np.array([0.1, 0.3, 0.8])
    Ks = kramers_wannier_dual_coupling(K)
    assert np.sinh(2 * K) * np.sinh(2 * Ks) == pytest.approx(np.ones(3))
    assert kramers_wannier_dual_coupling(Ks) == pytest.approx(K)
    assert ising_2d_critical_temperature(K_B) == pytest.approx(2.0 / np.log(1.0 + np.sqrt(2.0)))


def test_onsager_energy_and_heat_capacity_are_thermodynamically_consistent():
    model = Ising2DOnsager(J)
    for T in (150.0, 200.0, 260.0, 400.0):
        dT = 1e-3
        beta_f = lambda t: model.free_energy_per_spin(t) / (K_B * t)
        # u = d(beta f) / d(beta)
        u_numeric = (beta_f(T + dT) - beta_f(T - dT)) / (1 / (K_B * (T + dT)) - 1 / (K_B * (T - dT)))
        assert model.internal_energy_per_spin(T) == pytest.approx(u_numeric, rel=1e-6)
        c_numeric = (model.internal_energy_per_spin(T + dT) - model.internal_energy_per_spin(T - dT)) / (2 * dT)
        assert model.heat_capacity_per_spin(T) == pytest.approx(c_numeric, rel=1e-5)


def test_onsager_limits():
    model = Ising2DOnsager(J)
    Tc = model.critical_temperature
    assert model.internal_energy_per_spin(Tc) == pytest.approx(-np.sqrt(2.0) * J)
    assert model.free_energy_per_spin(1.0e6) / (K_B * 1.0e6) == pytest.approx(-np.log(2.0), rel=1e-4)
    assert model.internal_energy_per_spin(1.0) == pytest.approx(-2.0 * J)
    assert float(model.spontaneous_magnetization(0.5 * Tc)) > 0.99
    assert float(model.spontaneous_magnetization(1.01 * Tc)) == 0.0
    # logarithmic divergence: the peak sharpens as T -> Tc
    c = [model.heat_capacity_per_spin(Tc * (1 + d)) / K_B for d in (-1e-2, -1e-4, 1e-4, 1e-2)]
    assert c[1] > c[0] and c[2] > c[3] and min(c[1], c[2]) > 4.0
