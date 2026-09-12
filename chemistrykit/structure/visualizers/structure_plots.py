"""Molecular-structure plotting: 3D ball-and-stick geometry, and the bond-order/length correlation curve.

matplotlib-based, following the established chemistrykit visualizer
convention (:mod:`chemistrykit.quantum.visualizers.quantum_plots`): thin
functions that accept an optional ``ax`` (creating a new figure if
omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_molecule_3d", "plot_vsepr_geometry", "plot_bond_order_correlation"]

#: A small, deterministic color cycle for common elements, purely
#: cosmetic (CPK-inspired but not a faithful CPK palette).
_ELEMENT_COLORS = {
    "H": "whitesmoke",
    "C": "dimgray",
    "N": "royalblue",
    "O": "crimson",
    "F": "seagreen",
    "Cl": "limegreen",
    "S": "gold",
    "P": "darkorange",
    "B": "peru",
    "Xe": "orchid",
}


def plot_molecule_3d(molecule, ax=None, show_labels: bool = True, atom_size: float = 400.0, **kwargs):
    """Draw a molecule as a 3D ball-and-stick figure.

    Parameters
    ----------
    molecule : chemistrykit.structure.core.base_system.Molecule
    ax : matplotlib.axes.Axes3D, optional
        3D axes to draw on; a new 3D figure is created if omitted.
    show_labels : bool, default True
        If True, annotate each atom with its element symbol.
    atom_size : float, default 400.0
        Marker size (points^2) for the atom scatter.
    **kwargs
        Forwarded to the bond ``ax.plot`` calls.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
    coords = molecule.coordinates
    colors = [_ELEMENT_COLORS.get(s, "slateblue") for s in molecule.symbols]
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=atom_size, c=colors, edgecolors="black", depthshade=True)
    for i, j in molecule.bonds:
        ax.plot(*zip(coords[i], coords[j], strict=True), color="gray", linewidth=2, **kwargs)
    if show_labels:
        for symbol, (x, y, z) in zip(molecule.symbols, coords, strict=True):
            ax.text(x, y, z, symbol, fontsize=10, ha="center", va="center")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Molecular geometry")
    return ax


def plot_vsepr_geometry(steric_number: int, lone_pairs: int, ax=None, bond_length: float = 1.0, central_symbol: str = "C", ligand_symbol: str = "H"):
    """Build and plot a VSEPR-predicted geometry in one call.

    Parameters
    ----------
    steric_number : int
    lone_pairs : int
    ax : matplotlib.axes.Axes3D, optional
    bond_length : float, default 1.0
    central_symbol : str, default "C"
    ligand_symbol : str, default "H"

    Returns
    -------
    matplotlib.axes.Axes
    """
    from chemistrykit.structure.systems.vsepr import VSEPRGeometry, build_vsepr_molecule

    molecule = build_vsepr_molecule(steric_number, lone_pairs, bond_length, central_symbol, ligand_symbol)
    ax = plot_molecule_3d(molecule, ax=ax)
    shape = VSEPRGeometry(steric_number, lone_pairs).shape_name
    ax.set_title(f"VSEPR: {shape} (steric number {steric_number}, {lone_pairs} lone pair(s))")
    return ax


def plot_bond_order_correlation(single_bond_length: float, bond_orders=None, c: float = 0.71, ax=None, **kwargs):
    """Plot the Pauling bond-length-vs-bond-order correlation curve.

    Parameters
    ----------
    single_bond_length : float
        Reference single-bond length :math:`D(1)`.
    bond_orders : array-like of float, optional
        Bond orders to evaluate; defaults to ``np.linspace(0.5, 3.5, 100)``.
    c : float, default 0.71
        Pauling correlation constant.
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    from chemistrykit.structure.systems.bonding import bond_length_from_order

    if ax is None:
        _, ax = plt.subplots()
    if bond_orders is None:
        bond_orders = np.linspace(0.5, 3.5, 100)
    lengths = [bond_length_from_order(single_bond_length, n, c) for n in bond_orders]
    ax.plot(bond_orders, lengths, **kwargs)
    ax.set_xlabel("bond order n")
    ax.set_ylabel("bond length")
    ax.set_title("Pauling bond-order/bond-length correlation")
    return ax
