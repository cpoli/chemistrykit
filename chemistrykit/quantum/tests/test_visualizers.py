"""Smoke tests for chemistrykit.quantum.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes

from chemistrykit.quantum.systems.harmonic_oscillator import compare_harmonic_vs_morse
from chemistrykit.quantum.systems.huckel import HuckelSystem
from chemistrykit.quantum.systems.hydrogenlike import HydrogenLikeAtom
from chemistrykit.quantum.systems.particle_in_box import ParticleInBox1D
from chemistrykit.quantum.visualizers.quantum_plots import (
    plot_energy_levels,
    plot_harmonic_vs_morse_levels,
    plot_huckel_levels,
    plot_particle_in_box_wavefunctions,
    plot_radial_distribution,
)


def test_plot_energy_levels_returns_axes():
    ax = plot_energy_levels([1.0, 1.0, 2.0, 3.0, 3.0, 3.0])
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_energy_levels_with_labels():
    ax = plot_energy_levels([0.0, 1.0, 2.0], labels=["ground", "first", "second"])
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_particle_in_box_wavefunctions_returns_axes():
    box = ParticleInBox1D(length=1.0e-9)
    ax = plot_particle_in_box_wavefunctions(box, [1, 2, 3])
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 6  # 3 wavefunction curves + 3 dotted energy-offset reference lines


def test_plot_harmonic_vs_morse_levels_returns_axes():
    v, E_harmonic, E_morse = compare_harmonic_vs_morse(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19, v_max=5)
    ax = plot_harmonic_vs_morse_levels(v, E_harmonic, E_morse)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 2


def test_plot_radial_distribution_returns_axes():
    h_atom = HydrogenLikeAtom(Z=1)
    ax = plot_radial_distribution(h_atom, n=2, l=1)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_huckel_levels_returns_axes():
    benzene = HuckelSystem.cyclic_polyene(6)
    result = benzene.solve()
    ax = plot_huckel_levels(result, n_pi_electrons=6)
    assert isinstance(ax, matplotlib.axes.Axes)
