"""Tests for chemistrykit.photochem.systems.chain_reaction."""

import numpy as np
import pytest

from chemistrykit.photochem.systems.chain_reaction import chain_quantum_yield, hydrogen_chlorine_chain_network


def test_chain_quantum_yield_closed_form():
    assert chain_quantum_yield(k2=10.0, kt=100.0, H2=1.0, I_abs=1e-4) == pytest.approx(200.0)


def test_chain_quantum_yield_scales_as_inverse_sqrt_intensity():
    I = np.array([1e-6, 1e-4, 1e-2])
    phi = chain_quantum_yield(k2=5.0, kt=10.0, H2=0.5, I_abs=I)
    assert phi * np.sqrt(I) == pytest.approx(np.full(3, 2 * 5.0 * 0.5 / np.sqrt(10.0)))


def test_network_conserves_atoms_and_matches_steady_state_yield():
    j, k2, k3, kt = 1e-4, 10.0, 100.0, 100.0
    res = hydrogen_chlorine_chain_network(j, k2, k3, kt).integrate((0.0, 20.0), dt=1e-3, method="rk4")
    H2, Cl2 = res.concentration("H2"), res.concentration("Cl2")
    Cl, H, HCl = res.concentration("Cl"), res.concentration("H"), res.concentration("HCl")
    assert 2 * H2 + H + HCl == pytest.approx(np.full_like(H2, 2.0))
    assert 2 * Cl2 + Cl + HCl == pytest.approx(np.full_like(H2, 2.0))
    phi_numeric = np.gradient(HCl, res.t)[-1] / (j * Cl2[-1])
    assert phi_numeric > 100.0
    assert phi_numeric == pytest.approx(chain_quantum_yield(k2, kt, H2[-1], j * Cl2[-1]), rel=0.03)
