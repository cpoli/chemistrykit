Examples
========

This gallery walks through every public feature of ``chemistrykit.quantum``:
particle-in-a-box models (with the free-electron model of conjugated-dye
color); the quantum harmonic oscillator compared against the exact Morse
potential; the rigid rotor; hydrogen-like orbitals; Huckel molecular-orbital
theory and its aromaticity rule; a minimal variational treatment of H2+;
and Rayleigh-Schrodinger perturbation theory for the anharmonic oscillator.

Each script in this gallery is self-contained and can be run directly with
``python examples/quantum/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **particle_in_box** -- 1D/3D particle-in-a-box energy levels and
  wavefunctions, and Kuhn's free-electron model of a conjugated dye's
  UV-Vis absorption wavelength.
- **harmonic_oscillator** -- the quantum harmonic oscillator's evenly
  spaced vibrational levels, compared against the exact (anharmonic)
  Morse-potential levels for the same force constant.
- **rigid_rotor** -- rigid-rotor rotational energy levels, degeneracies,
  and the evenly spaced microwave absorption spectrum they predict.
- **hydrogenlike** -- hydrogen-like radial wavefunctions, radial
  distribution functions, and orbital energies for 1s/2s/2p/3d.
- **huckel** -- Huckel molecular-orbital theory for butadiene and
  benzene: building and diagonalizing the secular matrix, and checking
  Huckel's 4n+2 aromaticity rule against the computed spectrum.
- **hartree_fock** -- a minimal 2-Gaussian LCAO variational treatment of
  H2+: solving ``HC=SCE`` and variationally optimizing the orbital
  exponent to improve on a naive guess.
- **perturbation** -- Rayleigh-Schrodinger perturbation theory for the
  quartic anharmonic oscillator, checked against exact numerical
  diagonalization in a truncated basis.
