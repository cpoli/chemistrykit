"""Quantum-chemistry plotting: energy-level diagrams, wavefunctions, and orbital radial distributions.

matplotlib-based, following the established chemistrykit visualizer
convention (:mod:`chemistrykit.statmech.visualizers.statmech_plots`):
thin functions that accept an optional ``ax`` (creating a new figure if
omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "plot_energy_levels",
    "plot_particle_in_box_wavefunctions",
    "plot_harmonic_vs_morse_levels",
    "plot_radial_distribution",
    "plot_huckel_levels",
]


def plot_energy_levels(energies, ax=None, degeneracy_tol: float = 1e-9, labels=None, **kwargs):
    """Draw a simple energy-level diagram: one horizontal line per (possibly degenerate) level.

    Parameters
    ----------
    energies : array-like of float
        Energy of every state (degenerate levels may repeat).
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    degeneracy_tol : float, default 1e-9
        Levels closer together than this are drawn as a single line with
        multiple short segments (one per degenerate state).
    labels : sequence of str, optional
        Text label placed to the right of each distinct level.
    **kwargs
        Forwarded to ``ax.hlines``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    energies = np.sort(np.asarray(energies, dtype=np.float64))

    groups = []
    for e in energies:
        if groups and abs(e - groups[-1][0]) < degeneracy_tol:
            groups[-1].append(e)
        else:
            groups.append([e])

    for level_idx, group in enumerate(groups):
        n_degenerate = len(group)
        width = 0.6
        starts = np.linspace(-0.5, 0.5, n_degenerate + 1)
        for k in range(n_degenerate):
            x0 = starts[k] + 0.05
            x1 = starts[k + 1] - 0.05 if n_degenerate > 1 else width / 2
            x0 = -width / 2 if n_degenerate == 1 else x0
            ax.hlines(group[k], x0, x1, **kwargs)
        if labels is not None and level_idx < len(labels):
            ax.text(width / 2 + 0.1, group[0], str(labels[level_idx]), va="center")

    ax.set_xticks([])
    ax.set_ylabel("energy")
    ax.set_title("Energy-level diagram")
    return ax


def plot_particle_in_box_wavefunctions(box, n_values, ax=None, n_points: int = 300, scale: float = 1.0, **kwargs):
    """Plot particle-in-a-box wavefunctions, each vertically offset by its own energy level.

    Parameters
    ----------
    box : chemistrykit.quantum.systems.particle_in_box.ParticleInBox1D
    n_values : sequence of int
        Quantum numbers to plot.
    ax : matplotlib.axes.Axes, optional
    n_points : int, default 300
    scale : float, default 1.0
        Vertical scale factor applied to each wavefunction before adding
        its energy offset (purely for visual clarity).
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    x = np.linspace(0.0, box.length, n_points)
    for n in n_values:
        e = box.energy(n)
        psi = box.wavefunction(x, n)
        ax.plot(x, e + scale * psi, **kwargs)
        ax.axhline(e, color="gray", linestyle=":", linewidth=0.6)
    ax.set_xlabel("x")
    ax.set_ylabel("energy + scaled wavefunction")
    ax.set_title("Particle-in-a-box wavefunctions")
    return ax


def plot_harmonic_vs_morse_levels(v, harmonic_energies, morse_energies, ax=None):
    """Plot harmonic-oscillator vs. Morse-potential vibrational energies against `v`, side by side.

    Parameters
    ----------
    v : array-like of int
    harmonic_energies, morse_energies : array-like of float
        As returned by :func:`chemistrykit.quantum.systems.harmonic_oscillator.compare_harmonic_vs_morse`.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(v, harmonic_energies, "o-", label="harmonic")
    ax.plot(v, morse_energies, "s-", label="Morse")
    ax.set_xlabel("v")
    ax.set_ylabel("energy")
    ax.set_title("Harmonic vs. Morse vibrational levels")
    ax.legend()
    return ax


def plot_radial_distribution(hydrogen_atom, n: int, l: int, ax=None, r_max_bohr_radii: float = 20.0, n_points: int = 400, **kwargs):
    """Plot a hydrogen-like atom's radial distribution function P(r) = r^2 R_nl(r)^2.

    Parameters
    ----------
    hydrogen_atom : chemistrykit.quantum.systems.hydrogenlike.HydrogenLikeAtom
    n, l : int
    ax : matplotlib.axes.Axes, optional
    r_max_bohr_radii : float, default 20.0
        Upper plot limit, in units of the atom's Bohr radius.
    n_points : int, default 400
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    a0 = hydrogen_atom.bohr_radius
    r = np.linspace(1.0e-6 * a0, r_max_bohr_radii * a0, n_points)
    P = hydrogen_atom.radial_distribution_function(r, n, l)
    ax.plot(r / a0, P * a0, **kwargs)
    ax.set_xlabel("r / a0")
    ax.set_ylabel("P(r) * a0")
    ax.set_title(f"Radial distribution function (n={n}, l={l})")
    return ax


def plot_huckel_levels(result, n_pi_electrons=None, ax=None, degeneracy_tol: float = 1e-6):
    """Plot a Huckel molecular-orbital energy-level diagram, with electron occupation if given.

    Parameters
    ----------
    result : chemistrykit.quantum.core.base_system.EigenstateResult
        As returned by :meth:`chemistrykit.quantum.systems.huckel.HuckelSystem.solve`.
    n_pi_electrons : int, optional
        If given, fill the lowest orbitals with up/down arrows (2 electrons
        per orbital, Aufbau order) to show the ground-state occupation.
    ax : matplotlib.axes.Axes, optional
    degeneracy_tol : float, default 1e-6

    Returns
    -------
    matplotlib.axes.Axes
    """
    ax = plot_energy_levels(result.energies, ax=ax, degeneracy_tol=degeneracy_tol, color="black")
    if n_pi_electrons is not None:
        energies = np.sort(result.energies)
        n_full = n_pi_electrons // 2
        for i in range(n_full):
            ax.text(0.0, energies[i], "↑↓", ha="center", va="center")
        if n_pi_electrons % 2:
            ax.text(0.0, energies[n_full], "↑", ha="center", va="center")
    ax.set_title("Huckel molecular-orbital energy levels")
    return ax
