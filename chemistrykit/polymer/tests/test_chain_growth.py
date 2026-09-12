"""Tests for chemistrykit.polymer.systems.chain_growth against closed-form/known results.

The key cross-check (mirroring
chemistrykit.kinetics.systems.networks.ssa_intermediate_concentration's
validation style) is that numerically integrating free_radical_network
reproduces the steady-state-approximation's closed-form radical
concentration once the fast radical-concentration transient has decayed.
"""

import numpy as np
import pytest

from chemistrykit.polymer.systems.chain_growth import (
    free_radical_network,
    kinetic_chain_length,
    steady_state_radical_concentration,
    steady_state_rate_of_polymerization,
)

KD, F, KP, KT = 1.0e-5, 0.5, 1.0e3, 1.0e7
I0, M0 = 0.01, 5.0


def test_steady_state_radical_concentration_known_value():
    R_ss = steady_state_radical_concentration(KD, F, I0, KT)
    assert R_ss == pytest.approx(np.sqrt(F * KD * I0 / KT), rel=1e-12)


def test_steady_state_radical_concentration_scales_as_sqrt_of_initiator():
    R1 = steady_state_radical_concentration(KD, F, I=0.01, kt=KT)
    R2 = steady_state_radical_concentration(KD, F, I=0.04, kt=KT)
    assert R2 / R1 == pytest.approx(2.0, rel=1e-9)  # sqrt(0.04/0.01) = 2


def test_free_radical_network_numerical_integration_matches_ssa():
    net = free_radical_network(KD, F, KP, KT, I0, M0)
    result = net.integrate((0.0, 50.0), dt=1e-2, method="rk4")
    R_ss = steady_state_radical_concentration(KD, F, I0, KT)
    # After the fast initial transient, R(t) should track the (slowly
    # drifting) SSA value closely.
    assert result.concentration("R")[-1] == pytest.approx(R_ss, rel=0.05)


def test_free_radical_network_conserves_no_mass_violations():
    net = free_radical_network(KD, F, KP, KT, I0, M0)
    result = net.integrate((0.0, 50.0), dt=1e-2, method="rk4")
    # Concentrations must stay non-negative throughout (no numerical blow-up / NaN).
    assert np.all(np.isfinite(result.y))
    assert np.all(result.y >= -1e-8)


def test_free_radical_network_initiator_and_monomer_are_depleted():
    net = free_radical_network(KD, F, KP, KT, I0, M0)
    result = net.integrate((0.0, 50.0), dt=1e-2, method="rk4")
    assert result.concentration("I")[-1] < I0
    assert result.concentration("M")[-1] < M0


def test_free_radical_network_rejects_bad_mode():
    with pytest.raises(ValueError):
        free_radical_network(KD, F, KP, KT, I0, M0, mode="bogus")


def test_free_radical_network_combination_vs_disproportionation_dead_chain_count():
    net_comb = free_radical_network(KD, F, KP, KT, I0, M0, mode="combination")
    result_comb = net_comb.integrate((0.0, 50.0), dt=1e-2, method="rk4")
    net_disp = free_radical_network(KD, F, KP, KT, I0, M0, mode="disproportionation")
    result_disp = net_disp.integrate((0.0, 50.0), dt=1e-2, method="rk4")
    # Same number of termination events, but disproportionation makes two
    # dead chains per event vs. combination's one.
    assert result_disp.concentration("D")[-1] == pytest.approx(2.0 * result_comb.concentration("D")[-1], rel=1e-6)


def test_steady_state_rate_of_polymerization_scales_as_sqrt_of_initiator():
    Rp1 = steady_state_rate_of_polymerization(KD, F, KP, KT, I=0.01, M=M0)
    Rp2 = steady_state_rate_of_polymerization(KD, F, KP, KT, I=0.04, M=M0)
    assert Rp2 / Rp1 == pytest.approx(2.0, rel=1e-9)


def test_steady_state_rate_of_polymerization_is_linear_in_monomer():
    Rp1 = steady_state_rate_of_polymerization(KD, F, KP, KT, I=I0, M=1.0)
    Rp2 = steady_state_rate_of_polymerization(KD, F, KP, KT, I=I0, M=3.0)
    assert Rp2 / Rp1 == pytest.approx(3.0, rel=1e-9)


def test_kinetic_chain_length_equals_Rp_over_Ri():
    Rp = steady_state_rate_of_polymerization(KD, F, KP, KT, I=I0, M=M0)
    Ri = 2.0 * F * KD * I0
    nu = kinetic_chain_length(KD, F, KP, KT, I=I0, M=M0)
    assert nu == pytest.approx(Rp / Ri, rel=1e-9)
