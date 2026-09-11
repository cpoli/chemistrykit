"""Smoke tests for chemistrykit.thermo.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.constants import R
from chemistrykit.thermo.systems.equations_of_state import IdealGas, RedlichKwong, VanDerWaals
from chemistrykit.thermo.systems.equilibrium import fit_van_t_hoff, van_t_hoff_equilibrium_constant
from chemistrykit.thermo.systems.mixtures import BinaryIdealSolution
from chemistrykit.thermo.systems.phase_equilibria import ClausiusClapeyron
from chemistrykit.thermo.visualizers.thermo_plots import (
    plot_isotherms,
    plot_phase_boundary,
    plot_pxy_diagram,
    plot_van_t_hoff,
)


def test_plot_isotherms_returns_axes():
    Tc, Pc = 304.13, 7.3773e6
    ideal = IdealGas()
    vdw = VanDerWaals.from_critical_constants(Tc, Pc)
    rk = RedlichKwong.from_critical_constants(Tc, Pc)
    Vm = np.linspace(2.0 * vdw.b, 5.0e-4, 200)
    ax = plot_isotherms([ideal, vdw, rk], T=Tc, Vm_range=Vm)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 3


def test_plot_phase_boundary_returns_axes():
    model = ClausiusClapeyron(delta_h_vap=40700.0, T_ref=373.15, P_ref=101325.0)
    T = np.linspace(300.0, 373.15, 50)
    ax = plot_phase_boundary(model, T)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_van_t_hoff_returns_axes():
    T = np.linspace(280.0, 360.0, 8)
    K = van_t_hoff_equilibrium_constant(T, T_ref=298.15, K_ref=1.0, delta_h=50e3, R_gas=R)
    fit = fit_van_t_hoff(T, K)
    ax = plot_van_t_hoff(T, K, fit=fit)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_pxy_diagram_returns_axes():
    solution = BinaryIdealSolution(P_A_star=300.0, P_B_star=50.0)
    ax = plot_pxy_diagram(solution)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 2
