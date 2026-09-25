# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Initial repository scaffold: `pyproject.toml`, shared
  `chemistrykit.integrators` (RK4, leapfrog/velocity-Verlet, Yoshida4,
  adaptive Dormand-Prince, ported from physicskit), `chemistrykit.constants`
  (chemical constants from `scipy.constants` plus a built-in periodic-table
  data table), CI/pre-commit/readthedocs config, and root docs.
- `chemistrykit.kinetics` — the first complete domain subpackage, and the
  template the remaining 13 domains will follow: elementary integrated
  rate laws and half-lives (zero/first/second order); the Arrhenius
  equation and activation-energy fitting; Michaelis-Menten enzyme
  kinetics with Lineweaver-Burk linearization and competitive/
  noncompetitive inhibition; a general stoichiometry-matrix
  reaction-network engine (parallel, consecutive, reversible, and
  steady-state-approximation chains) integrated via
  `chemistrykit.integrators`; and the Brusselator oscillating reaction
  network. Includes a matching `examples/kinetics/` sphinx-gallery and
  `docs/source/api/kinetics.rst` / `docs/source/examples/kinetics.rst`.
- `chemistrykit.thermo` — equations of state (ideal gas, van der Waals,
  Redlich-Kwong) sharing an `EquationOfState` interface; Clausius-Clapeyron
  phase boundaries; reaction equilibrium (Kp/Kc, van't Hoff, and a
  Gibbs-energy-minimization equilibrium-composition solver); and
  Raoult's/Henry's law mixtures with colligative properties.
- `chemistrykit.solutions` — pH/pOH and weak acid/base equilibria with
  Henderson-Hasselbalch buffers; strong/weak acid-base titration curves
  with numerical equivalence-point detection; Ksp solubility equilibria
  and the common-ion effect; and Debye-Huckel activity coefficients.
- `chemistrykit.md` — the Lennard-Jones fluid in reduced units (periodic
  boundary conditions, a Verlet neighbor list, pressure, radial
  distribution function), Numba-accelerated force evaluation;
  Morse/Buckingham/harmonic bonded potentials; and velocity-rescaling/
  Nose-Hoover thermostats.
- `chemistrykit.statmech` — translational/rotational/vibrational
  partition functions and their derived thermodynamic functions; the
  Maxwell-Boltzmann speed distribution; and a canonical-ensemble
  lattice-gas adsorption model.
- `chemistrykit.quantum` — particle-in-a-box models (with the
  free-electron model of conjugated-dye color); the quantum harmonic
  oscillator vs. the exact Morse potential; the rigid rotor;
  hydrogen-like orbitals; Huckel molecular-orbital theory and its 4n+2
  aromaticity rule; a minimal variational treatment of H2+; and
  Rayleigh-Schrodinger perturbation theory for the anharmonic oscillator.
- `chemistrykit.spectro` — the Beer-Lambert absorbance law and its
  stray-light deviation from linearity; rigid-rotor rotational spectra
  with isotope shifts; harmonic vs. Morse vibrational band positions plus
  a Wilson GF-matrix triatomic normal-mode calculation; Franck-Condon
  vibronic progressions; a first-order NMR multiplet simulator; and
  shared Gaussian/Lorentzian/Voigt lineshape utilities.
- `chemistrykit.structure` — a lightweight `Molecule` container; VSEPR
  geometry prediction with real 3D coordinate generation; point-group
  determination from 3D coordinates (via direct geometric symmetry-
  operation testing) and character tables; bond order from the Pauling
  length correlation and Huckel-theory MO coefficients; and
  formal-charge/oxidation-state assignment from a Lewis structure.
- `chemistrykit.electrochem` — the Nernst equation for standard and
  concentration cells with Debye-Huckel activity corrections; a curated
  standard-reduction-potential table with redox-couple balancing;
  Butler-Volmer electrode kinetics and Tafel-plot linearization; galvanic
  vs. electrolytic cells and Faraday's laws of electrolysis; and a
  simplified constant-current battery discharge model with Peukert's-law
  rate dependence.
- `chemistrykit.photochem` — Jablonski-diagram excited-state kinetics
  built on `chemistrykit.kinetics`'s reaction-network engine;
  fluorescence/phosphorescence quantum yields and the photochemical
  quantum yield via `chemistrykit.spectro`'s Beer-Lambert law;
  Stern-Volmer quenching with a static-vs-dynamic diagnostic; and
  photostationary-state kinetics for a two-state photoswitch.
- `chemistrykit.surface` — Langmuir, Freundlich, and BET adsorption
  isotherms sharing an `AdsorptionIsotherm` interface, each with its
  standard linearization for fitting parameters from data;
  Langmuir-Hinshelwood single- and dual-site surface-reaction kinetics;
  and a turnover-frequency/rate-enhancement catalysis model built on
  `chemistrykit.kinetics`'s Arrhenius equation.
- `chemistrykit.polymer` — ideal random-walk chain statistics and the
  Flory exponent for real chains under theta/good/poor solvent
  conditions, sharing a `PolymerChainModel` interface;
  molecular-weight-distribution statistics and the closed-form
  Flory-Schulz distribution; step-growth kinetics via the Carothers
  equation; and chain-growth/free-radical polymerization kinetics built
  on `chemistrykit.kinetics`'s reaction-network engine.
- `chemistrykit.crystal` — the 7 crystal systems and the general
  lattice-parameter unit-cell-volume formula; hard-sphere packing
  (packing fraction, coordination number, atoms per cell) for
  SC/BCC/FCC/ideal-HCP lattices sharing a `LatticePacking` interface;
  ionic-crystal lattice energy via the Born-Lande and Kapustinskii
  equations (`LatticeEnergyModel` interface), backed by a genuinely
  converging Evjen-method lattice summation for the NaCl Madelung
  constant (verified against the literature value ~1.747565); Bragg's
  law and powder-XRD peak positions with structure factors computed from
  atomic positions, including systematic absences for BCC/FCC verified
  combinatorially; and Schottky/Frenkel point-defect equilibrium. This
  completes all 14 domains from `chemistrykit-spec.md`.
- `chemistrykit.analytical` — redox and complexometric (EDTA)
  titration-curve simulation, extending (not duplicating) the acid-base
  titrations in `chemistrykit.solutions`; chromatographic plate theory
  and the van Deemter equation (optimum flow velocity, resolution,
  selectivity), reusing `chemistrykit.spectro`'s Gaussian lineshape for
  simulated chromatograms; linear-regression calibration curves with
  IUPAC-convention limits of detection/quantitation; and
  propagation-of-uncertainty formulas plus Dixon's Q-test outlier
  rejection with the Rorabacher (1991) critical-value table.

### Changed

- Dropped Python 3.9 support: `requires-python` is now `>=3.10` (matching
  physicskit), and the CI matrix, classifiers, ruff `target-version`, and
  mypy `python_version` follow. Current mypy and NumPy's own type stubs
  no longer support 3.9.
- `chemistrykit.crystal.powder_xrd_peaks` returns one peak per family of
  symmetry-equivalent planes, labelled by its canonical `h >= k >= l`
  member, instead of one per ordering (BCC's first line was three
  coincident `(1, 1, 0)`/`(1, 0, 1)`/`(0, 1, 1)` peaks). `XRDPeak` gained
  a `multiplicity` field (e.g. 12 for {110}).
- `VSEPRGeometry` now requires at least one bonded ligand
  (`lone_pairs < steric_number`) and has shape names for every remaining
  combination (e.g. HF as AX1E3, "linear"; AX3E3, "T-shaped").

### Fixed

- `chemistrykit.electrochem`: `exchange_current_density` applied the
  charge-transfer coefficient to the wrong species
  (`C_ox^(1-alpha) C_red^alpha`) relative to `butler_volmer_current_density`'s
  anodic `alpha`; it is now `C_ox^alpha C_red^(1-alpha)`, and every
  function in the module documents `alpha` as the anodic coefficient.
- `chemistrykit.polymer.free_radical_network` folded the initiator
  efficiency `f` into the decomposition rate constant (`f*kd`), so the
  initiator was consumed `1/f` times too slowly. It now decomposes at
  `kd` and yields `2f` radicals per event, keeping `R_i = 2 f kd [I]`.
- `chemistrykit.spectro.rotational_spectrum` weighted line intensities by
  `(J+1)(2J+1)`, double-counting the level degeneracy; the
  degeneracy-averaged transition dipole `(J+1)/(2J+1)` cancels it,
  giving `(J+1) exp(-E_J/kT)`.
- `chemistrykit.spectro.apparent_absorbance_with_stray_light` omitted
  the stray light reaching the detector during the blank reading, so a
  blank read a nonzero absorbance; it now uses
  `log10((P0 + Ps) / (P + Ps))`.
- `chemistrykit.analytical.EDTATitration` computed free `[M]` as
  `C_M - [MY]`, which cancels catastrophically at and past equivalence;
  it now solves the mass-action quadratic for `[M]` directly.
- The weak-acid/weak-base titrations in `chemistrykit.solutions` searched
  a fixed pH 0-14 bracket that concentrated solutions step outside of,
  and `find_positive_root`'s absolute tolerance (2e-12) exceeded a basic
  solution's `[H+]`. Both now use concentration-derived brackets and
  relative tolerances.
- `chemistrykit.constants.R` is computed as `NA * K_B` (exact by SI
  definition) rather than taken from `scipy.constants.R`, which older
  SciPy releases store truncated.
- The shared `chemistrykit.integrators` kernels that take an njit
  callback are no longer `@njit(cache=True)`: their signature includes the
  callback's per-instance dispatcher type, so the on-disk cache only
  accumulated dead entries and could raise
  `ReferenceError: underlying object has vanished` when re-saved.
- `chemistrykit.integrators.dopri5_integrate` (and so
  `ReactionNetwork.integrate(method="dopri5")`) now emits a
  `RuntimeWarning` when `max_steps` runs out before `t_end`, instead of
  silently returning a truncated trajectory.
- `ZeroOrder.concentration` no longer goes negative past the exhaustion
  time `C0/k` (it stays at 0), and `ZeroOrder.rate` drops to 0 there.
- `VSEPRGeometry.shape_name` raised a bare `KeyError` for valid inputs
  such as `(steric_number=4, lone_pairs=3)`.

[Unreleased]: https://github.com/cpoli/chemistrykit/compare/main...HEAD
