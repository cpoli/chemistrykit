# chemistrykit

[![CI](https://github.com/chemistrykit/chemistrykit/actions/workflows/ci.yml/badge.svg)](https://github.com/chemistrykit/chemistrykit/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/chemistrykit/badge/?version=latest)](https://chemistrykit.readthedocs.io)
[![PyPI](https://img.shields.io/pypi/v/chemistrykit.svg)](https://pypi.org/project/chemistrykit/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/chemistrykit.svg)](pyproject.toml)

A unified numerical toolkit for computational chemistry: pure-numerical
models across reaction kinetics, thermodynamics, quantum chemistry,
electrochemistry, spectroscopy, and more -- no cheminformatics
dependencies (no RDKit/ASE/PySCF/OpenMM), sharing common ODE integrators,
chemical constants, and a consistent NumPy-based API. Conventionally
imported as `ck`. chemistrykit clones the architecture and engineering
conventions of the sibling project
[physicskit](https://github.com/physicskit/physicskit) (`pk`).

All 14 domains from `chemistrykit-spec.md`'s build plan are implemented
-- see [Subpackages](#subpackages) for the full list, or browse the docs
at <https://chemistrykit.readthedocs.io>.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Quick start

```python
import chemistrykit as ck
import numpy as np

# Closed-form first-order decay and its half-life
rate_law = ck.kinetics.FirstOrder(k=0.1, C0=1.0)
print(rate_law.half_life())  # ln(2) / k

# A general stoichiometric reaction network: A -> B -> C
network = ck.kinetics.StoichiometricNetwork.consecutive(k1=1.0, k2=0.3, A0=1.0)
result = network.integrate((0.0, 10.0), dt=1e-3, method="rk4")
print(result.concentration("B")[-1])
```

## Subpackages

Domain subpackages, each with runnable examples linked below:

- [`chemistrykit.kinetics`](https://chemistrykit.readthedocs.io/en/latest/examples/kinetics.html) -- reaction kinetics: integrated rate laws and half-lives, the Arrhenius equation and activation-energy fitting, Michaelis-Menten enzyme kinetics with Lineweaver-Burk linearization and inhibition, a general stoichiometry-matrix reaction-network engine (parallel/consecutive/reversible/steady-state-approximation chains), and the Brusselator oscillating reaction network.
- [`chemistrykit.thermo`](https://chemistrykit.readthedocs.io/en/latest/examples/thermo.html) -- chemical thermodynamics: equations of state (ideal gas, van der Waals, Redlich-Kwong), Clausius-Clapeyron phase boundaries, reaction equilibrium (Kp/Kc, van't Hoff, and a Gibbs-energy-minimization equilibrium-composition solver), and Raoult's/Henry's law mixtures with colligative properties.
- [`chemistrykit.solutions`](https://chemistrykit.readthedocs.io/en/latest/examples/solutions.html) -- solution chemistry: pH/pOH and weak acid/base equilibria with Henderson-Hasselbalch buffers, acid-base titration curves, Ksp solubility equilibria and the common-ion effect, and Debye-Huckel activity coefficients.
- [`chemistrykit.md`](https://chemistrykit.readthedocs.io/en/latest/examples/md.html) -- molecular dynamics and force fields: the Lennard-Jones fluid in reduced units (periodic boundary conditions, a Verlet neighbor list, pressure, g(r)), Morse/Buckingham/harmonic bonded potentials, and velocity-rescaling/Nose-Hoover thermostats.
- [`chemistrykit.statmech`](https://chemistrykit.readthedocs.io/en/latest/examples/statmech.html) -- statistical mechanics of molecules: translational/rotational/vibrational partition functions and their thermodynamic functions, the Maxwell-Boltzmann speed distribution, and a canonical-ensemble lattice-gas adsorption model.
- [`chemistrykit.quantum`](https://chemistrykit.readthedocs.io/en/latest/examples/quantum.html) -- quantum chemistry: particle-in-a-box models (with the free-electron model of conjugated-dye color); the quantum harmonic oscillator vs. the exact Morse potential; the rigid rotor; hydrogen-like orbitals; Huckel molecular-orbital theory and its 4n+2 aromaticity rule; a minimal variational treatment of H2+; and Rayleigh-Schrodinger perturbation theory for the anharmonic oscillator.
- [`chemistrykit.spectro`](https://chemistrykit.readthedocs.io/en/latest/examples/spectro.html) -- spectroscopy: the Beer-Lambert absorbance law and its stray-light deviation from linearity; rigid-rotor rotational spectra with isotope shifts; harmonic vs. Morse vibrational band positions plus a Wilson GF-matrix triatomic normal-mode calculation; Franck-Condon vibronic progressions; and a first-order NMR multiplet simulator.
- [`chemistrykit.structure`](https://chemistrykit.readthedocs.io/en/latest/examples/structure.html) -- molecular structure and bonding: a lightweight `Molecule` container; VSEPR geometry prediction with real 3D coordinate generation; point-group determination from 3D coordinates and character tables; bond order from the Pauling length correlation and Huckel-theory MO coefficients; and formal-charge/oxidation-state assignment from a Lewis structure.
- [`chemistrykit.electrochem`](https://chemistrykit.readthedocs.io/en/latest/examples/electrochem.html) -- electrochemistry: the Nernst equation for standard and concentration cells with Debye-Huckel activity corrections; a curated standard-reduction-potential table with redox-couple balancing; Butler-Volmer electrode kinetics and Tafel-plot linearization; galvanic vs. electrolytic cells and Faraday's laws of electrolysis; and a simplified constant-current battery discharge model with Peukert's-law rate dependence.
- [`chemistrykit.photochem`](https://chemistrykit.readthedocs.io/en/latest/examples/photochem.html) -- photochemistry: Jablonski-diagram excited-state kinetics built on `chemistrykit.kinetics`'s reaction-network engine; fluorescence/phosphorescence quantum yields and the photochemical quantum yield via Beer-Lambert; Stern-Volmer quenching with a static-vs-dynamic diagnostic; and photostationary-state kinetics for a two-state photoswitch.
- [`chemistrykit.surface`](https://chemistrykit.readthedocs.io/en/latest/examples/surface.html) -- surface chemistry and catalysis: Langmuir, Freundlich, and BET adsorption isotherms with their standard linearizations for fitting parameters from data; Langmuir-Hinshelwood single- and dual-site surface-reaction kinetics; and a turnover-frequency/rate-enhancement catalysis model built on `chemistrykit.kinetics`'s Arrhenius equation.
- [`chemistrykit.polymer`](https://chemistrykit.readthedocs.io/en/latest/examples/polymer.html) -- polymer chemistry: ideal random-walk chain statistics and the Flory exponent for real chains under theta/good/poor solvent conditions; molecular-weight-distribution statistics and the closed-form Flory-Schulz distribution; step-growth kinetics via the Carothers equation; and chain-growth/free-radical polymerization kinetics built on `chemistrykit.kinetics`'s reaction-network engine.
- [`chemistrykit.crystal`](https://chemistrykit.readthedocs.io/en/latest/examples/crystal.html) -- crystallography and solid-state chemistry: the 7 crystal systems and general unit-cell volume; hard-sphere packing (packing fraction, coordination number) for SC/BCC/FCC/HCP lattices; ionic-crystal lattice energy via the Born-Lande and Kapustinskii equations, backed by a genuinely converging (Evjen-method) numerical Madelung constant; Bragg's law and powder-XRD peak positions with structure factors and systematic absences; and Schottky/Frenkel point-defect equilibrium.
- [`chemistrykit.analytical`](https://chemistrykit.readthedocs.io/en/latest/examples/analytical.html) -- analytical chemistry: redox and complexometric (EDTA) titration-curve simulation with equivalence-point detection, alongside `chemistrykit.solutions`'s acid-base titrations; chromatographic plate theory and the van Deemter equation (resolution, selectivity); linear-regression calibration curves with IUPAC-convention limits of detection/quantitation; and propagation-of-uncertainty formulas plus Dixon's Q-test for outlier rejection.

Shared infrastructure, used across the subpackages above rather than
standalone toolkits:

- `chemistrykit.constants` -- chemical constants (R, NA, k_B, Faraday's
  constant, ...) from `scipy.constants`, plus a small built-in
  periodic-table data table and a few well-defined unit conversions.
- `chemistrykit.integrators` -- shared numerical ODE integrators (RK4,
  leapfrog, Yoshida4, adaptive Dormand-Prince) used across the other
  subpackages.

## Test

Tests live alongside each subpackage, at `chemistrykit/<name>/tests/`.

```bash
pytest                                              # everything
pytest chemistrykit/kinetics/tests                  # a single subpackage

# docstring examples, across every subpackage:
MPLBACKEND=Agg pytest --doctest-modules chemistrykit \
    --ignore-glob="*/tests/*"
```

Both commands, plus `ruff check`/`ruff format --check`, run in CI on
every PR (`.github/workflows/ci.yml`) across Python 3.9-3.12 on Linux and
macOS. See [CONTRIBUTING.md](CONTRIBUTING.md) before opening a PR.

## Docs

Built docs are hosted at <https://chemistrykit.readthedocs.io>. To build
locally:

```bash
pip install -e ".[docs]"
cd docs && make html
```

## Citation

If you use chemistrykit in your research, please cite it — see
[CITATION.cff](CITATION.cff).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Please note that this project
follows the [Contributor Covenant](CODE_OF_CONDUCT.md).

## License

MIT -- see [LICENSE](LICENSE).
