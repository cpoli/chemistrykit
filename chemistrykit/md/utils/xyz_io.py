"""Reading and writing multi-frame XYZ trajectory files, for VMD, OVITO, Jmol and similar viewers.

Each XYZ frame is an atom-count line, a free-form comment line, and one
``symbol x y z`` line per atom; a trajectory is frames concatenated in one
file. When a cubic periodic box is known, the comment line uses the
*extended XYZ* convention (``Lattice="L 0 0 0 L 0 0 0 L" Properties=...``),
which OVITO and ASE read as the simulation cell; plain-XYZ readers such as
VMD simply ignore it.

Viewers assume coordinates in angstroms, while :mod:`chemistrykit.md`
works in whatever length unit the system was built in -- Lennard-Jones
reduced units (sigma) for :class:`~chemistrykit.md.LJFluid`, metres for
SI-parameterized molecules. Pass `scale` to convert: sigma in angstroms
(3.405 for argon) for reduced units, or ``1e10`` for metres.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

__all__ = ["XYZTrajectory", "write_xyz", "read_xyz"]


@dataclass
class XYZTrajectory:
    """Frames read back by :func:`read_xyz`."""

    symbols: list
    """list of str: Element symbol (or label) of each atom."""

    positions: np.ndarray
    """ndarray, shape (n_frames, n_atoms, 3): Coordinates, divided by the `scale` passed to :func:`read_xyz`."""

    comments: list = field(default_factory=list)
    """list of str: Comment line of each frame."""


def write_xyz(path, trajectory, symbols="Ar", scale: float = 1.0, box_length=None, comment: str = "") -> None:
    """Write one or more frames to an XYZ file.

    Parameters
    ----------
    path : str or path-like
        Output file (overwritten).
    trajectory : array-like or MDResult
        Positions of shape ``(n_atoms, n_dim)`` or ``(n_frames, n_atoms,
        n_dim)`` with ``n_dim`` 2 or 3 (2D systems get ``z = 0``), or an
        :class:`chemistrykit.md.MDResult`, whose frame times and
        ``box_length`` are written too.
    symbols : str or sequence of str, default "Ar"
        One symbol for every atom, or one per atom.
    scale : float, default 1.0
        Factor applied to coordinates (and box length) on writing, e.g.
        to convert to angstroms.
    box_length : float, optional
        Cubic box side (unscaled), written as an extended-XYZ lattice;
        taken from an `MDResult` when not given.
    comment : str, default ""
        Extra text for every frame's comment line (written as
        ``Comment="..."`` when the line is extended XYZ).

    Examples
    --------
    >>> import os, tempfile
    >>> import numpy as np
    >>> frames = np.random.default_rng(0).random((3, 4, 3))
    >>> path = os.path.join(tempfile.mkdtemp(), "traj.xyz")
    >>> write_xyz(path, frames, symbols="Ar", box_length=1.0)
    >>> traj = read_xyz(path)
    >>> traj.positions.shape, bool(np.allclose(traj.positions, frames, atol=1e-8))
    ((3, 4, 3), True)
    """
    times = None
    if hasattr(trajectory, "positions"):
        times = np.asarray(trajectory.t)
        if box_length is None:
            box_length = trajectory.box_length
        trajectory = trajectory.positions
    frames = np.asarray(trajectory, dtype=np.float64)
    if frames.ndim == 2:
        frames = frames[None]
    if frames.ndim != 3 or frames.shape[2] not in (2, 3):
        raise ValueError("positions must have shape (n_atoms, 2|3) or (n_frames, n_atoms, 2|3)")
    if frames.shape[2] == 2:
        frames = np.concatenate([frames, np.zeros(frames.shape[:2] + (1,))], axis=2)
    n_atoms = frames.shape[1]
    labels = [symbols] * n_atoms if isinstance(symbols, str) else list(symbols)
    if len(labels) != n_atoms:
        raise ValueError(f"got {len(labels)} symbols for {n_atoms} atoms")

    header = []
    if box_length is not None:
        L = box_length * scale
        header.append(f'Lattice="{L:.10g} 0 0 0 {L:.10g} 0 0 0 {L:.10g}" Properties=species:S:1:pos:R:3')
    if comment:
        # extended-XYZ comment lines are key=value pairs only
        header.append(f'Comment="{comment.replace(chr(34), chr(39))}"' if header else comment)
    with open(os.fspath(path), "w") as fh:
        for i, frame in enumerate(frames * scale):
            line = header + ([f"Time={times[i]:.10g}"] if times is not None else [])
            fh.write(f"{n_atoms}\n{' '.join(line)}\n")
            for label, (x, y, z) in zip(labels, frame, strict=True):
                fh.write(f"{label} {x:.10f} {y:.10f} {z:.10f}\n")


def read_xyz(path, scale: float = 1.0) -> XYZTrajectory:
    """Read every frame of an XYZ (or extended XYZ) file.

    Parameters
    ----------
    path : str or path-like
    scale : float, default 1.0
        Coordinates are divided by this on reading (the inverse of
        :func:`write_xyz`'s `scale`).

    Returns
    -------
    XYZTrajectory

    Raises
    ------
    ValueError
        If the file is truncated or frames have different atom counts or
        symbols.
    """
    with open(os.fspath(path)) as fh:
        lines = fh.read().splitlines()
    symbols: list = []
    frames: list[list[list[float]]] = []
    comments: list[str] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        if i + 2 + n > len(lines):
            raise ValueError(f"truncated frame starting at line {i + 1}")
        comments.append(lines[i + 1])
        rows = [lines[j].split() for j in range(i + 2, i + 2 + n)]
        frame_symbols = [r[0] for r in rows]
        if not frames:
            symbols = frame_symbols
        elif frame_symbols != symbols:
            raise ValueError(f"frame {len(frames)} has different atoms from frame 0")
        frames.append([[float(v) for v in r[1:4]] for r in rows])
        i += 2 + n
    return XYZTrajectory(symbols=symbols, positions=np.array(frames, dtype=np.float64).reshape(len(frames), -1, 3) / scale, comments=comments)
