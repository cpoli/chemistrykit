"""Smoke tests for chemistrykit.structure.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes

from chemistrykit.structure.systems.vsepr import build_vsepr_molecule
from chemistrykit.structure.visualizers.structure_plots import plot_bond_order_correlation, plot_molecule_3d, plot_vsepr_geometry


def test_plot_molecule_3d_returns_axes():
    methane = build_vsepr_molecule(4, 0, bond_length=1.09, central_symbol="C", ligand_symbol="H")
    ax = plot_molecule_3d(methane)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_molecule_3d_without_labels():
    water = build_vsepr_molecule(4, 2, bond_length=0.96, central_symbol="O", ligand_symbol="H")
    ax = plot_molecule_3d(water, show_labels=False)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_vsepr_geometry_returns_axes():
    ax = plot_vsepr_geometry(steric_number=5, lone_pairs=1, central_symbol="S", ligand_symbol="F")
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_bond_order_correlation_returns_axes():
    ax = plot_bond_order_correlation(single_bond_length=1.54)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 1
