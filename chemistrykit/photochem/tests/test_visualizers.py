"""Smoke tests for chemistrykit.photochem.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.photochem.systems.jablonski import jablonski_network
from chemistrykit.photochem.systems.photostationary_state import photostationary_state, photoswitch_network
from chemistrykit.photochem.systems.stern_volmer import fit_stern_volmer, stern_volmer_ratio
from chemistrykit.photochem.visualizers.photochem_plots import (
    plot_photostationary_approach,
    plot_state_populations,
    plot_stern_volmer,
)


def test_plot_state_populations_returns_axes():
    net = jablonski_network(kf=2.0, kic=1.0, kisc=0.5, kp=0.3, kic_T=0.2)
    result = net.integrate((0.0, 10.0), dt=1e-2, method="rk4")
    ax = plot_state_populations(result)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_stern_volmer_returns_axes():
    Q = np.array([0.0, 0.01, 0.02, 0.03])
    ratio = stern_volmer_ratio(Ksv=30.0, Q=Q)
    fit = fit_stern_volmer(Q, ratio)
    ax = plot_stern_volmer(Q, ratio, fit=fit)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_photostationary_approach_returns_axes():
    k_AB, k_BA = 1.5, 0.4
    net = photoswitch_network(k_AB, k_BA, A0=1.0)
    result = net.integrate((0.0, 50.0), dt=1e-2, method="rk4")
    pss = photostationary_state(k_AB, k_BA)
    ax = plot_photostationary_approach(result, pss)
    assert isinstance(ax, matplotlib.axes.Axes)
