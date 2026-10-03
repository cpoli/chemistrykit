# chemistrykit

| | |
|:--|:-:|
| Package | [![PyPI version](https://img.shields.io/pypi/v/chemistrykit)](https://pypi.org/project/chemistrykit/) [![Python versions](https://img.shields.io/pypi/pyversions/chemistrykit)](https://pypi.org/project/chemistrykit/) [![DOI](https://zenodo.org/badge/1377015934.svg)](https://zenodo.org/badge/latestdoi/1377015934) |
| Quality | [![License](https://img.shields.io/github/license/cpoli/chemistrykit)](https://github.com/cpoli/chemistrykit/blob/main/LICENSE) [![CI](https://github.com/cpoli/chemistrykit/actions/workflows/ci.yml/badge.svg)](https://github.com/cpoli/chemistrykit/actions/workflows/ci.yml) [![Coverage](https://img.shields.io/codecov/c/github/cpoli/chemistrykit)](https://codecov.io/gh/cpoli/chemistrykit) [![Coverage (manual)](https://img.shields.io/badge/coverage-96%25-brightgreen)](#coverage) |
| Documentation | [![Docs](https://img.shields.io/badge/docs-cpoli.github.io%2Fchemistrykit-blue)](https://cpoli.github.io/chemistrykit/) |
| Code style | [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) |
| Downloads | [![Downloads](https://static.pepy.tech/badge/chemistrykit)](https://pepy.tech/project/chemistrykit) [![Downloads/Month](https://static.pepy.tech/badge/chemistrykit/month)](https://pepy.tech/project/chemistrykit) |
| Community | [![GitHub Stars](https://img.shields.io/github/stars/cpoli/chemistrykit?style=social)](https://github.com/cpoli/chemistrykit) [![GitHub Forks](https://img.shields.io/github/forks/cpoli/chemistrykit?style=social)](https://github.com/cpoli/chemistrykit) [![Contributors](https://img.shields.io/github/contributors/cpoli/chemistrykit)](https://github.com/cpoli/chemistrykit/graphs/contributors) [![Last Commit](https://img.shields.io/github/last-commit/cpoli/chemistrykit)](https://github.com/cpoli/chemistrykit/commits/main) |

**See the equations of chemistry at work.**
chemistrykit is a Python toolkit for learning and teaching computational
chemistry, from Arrhenius kinetics and oscillating reactions to Hückel
aromaticity, Hartree-Fock, NMR multiplets and powder XRD. Every model is
built from its first-principles formula, every simulation returns an
inspectable result dataclass, every domain ships plotting helpers, and each domain's docs
walk through the field's breakthroughs in historical order, each one
linked to runnable code that reproduces it.

![Belousov-Zhabotinsky oscillations, ethanol's 1H NMR spectrum, and van der Waals isotherms of CO2, all drawn with chemistrykit](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_hero.png)

- **For students:** watch the Belousov-Zhabotinsky reaction oscillate,
  see where ethanol's triplet and quartet come from, find the
  ethanol-water azeotrope, all in a few lines each.
- **For instructors:** 14 domains, one consistent API, and over 200
  gallery examples, each downloadable as a Python script or Jupyter
  notebook, ready to hand out as course material.
- **Pure-numerical, built on numpy/scipy:** no cheminformatics
  dependencies (no RDKit/ASE/PySCF/OpenMM); library routines under the
  hood, with algorithms hand-rolled only where the steps themselves are
  what you're learning (see [Design](#design)).

chemistrykit is part of a family of packages --
[physicskit](https://github.com/cpoli/physicskit),
[mathematicskit](https://github.com/cpoli/mathematicskit) and
**chemistrykit** -- that share the same architecture, API conventions,
and history-driven documentation.

## Install

```bash
pip install chemistrykit
```

Then `import chemistrykit as ck`. Runtime dependencies are numpy, scipy,
matplotlib, [numba](https://numba.pydata.org/) (which compiles the inner
loops of `kinetics` reaction networks, `md` pair forces, and the shared
integrators), plotly, sympy, and tqdm.

For development: `pip install -e ".[dev]"` (see [CONTRIBUTING.md](CONTRIBUTING.md)).

## Quick start

```python
import chemistrykit as ck
from chemistrykit.kinetics.visualizers.kinetics_plots import plot_concentration_vs_time

# Closed-form first-order decay and its half-life
print(ck.kinetics.FirstOrder(k=0.1, C0=1.0).half_life())  # ln(2) / k = 6.93

# A general stoichiometric reaction network: A -> B -> C
network = ck.kinetics.StoichiometricNetwork.consecutive(k1=1.0, k2=0.3, A0=1.0)
result = network.integrate((0.0, 15.0), dt=1e-3, method="rk4")
B = result.concentration("B")
print(f"[B] peaks at t = {result.t[B.argmax()]:.2f}")  # ln(k1/k2) / (k1 - k2) = 1.72

plot_concentration_vs_time(result)
```

## Subpackages

Domain subpackages, each with runnable examples linked below:

- [`chemistrykit.kinetics`](https://cpoli.github.io/chemistrykit/api/gallery/kinetics/) -- reaction kinetics: integrated rate laws and half-lives, the Arrhenius equation and activation-energy fitting, Michaelis-Menten enzyme kinetics with Lineweaver-Burk linearization and inhibition, a general stoichiometry-matrix reaction-network engine (parallel/consecutive/reversible/steady-state-approximation chains) integrated by numba-compiled right-hand sides on `chemistrykit.integrators`, the Brusselator and Oregonator (Belousov-Zhabotinsky) oscillating reaction networks, Gillespie stochastic simulation (hand-rolled), Lindemann-Hinshelwood unimolecular falloff, and Semenov chain-branching explosion limits for H2/O2.

  ![Brusselator limit cycle, Oregonator oscillations, and Gillespie stochastic paths](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_kinetics.png)

- [`chemistrykit.thermo`](https://cpoli.github.io/chemistrykit/api/gallery/thermo/) -- chemical thermodynamics: equations of state (ideal gas, van der Waals, Redlich-Kwong), Clausius-Clapeyron and Antoine vapor-pressure curves, reaction equilibrium (Kp/Kc, van't Hoff, and a Gibbs-energy-minimization equilibrium-composition solver via `scipy.optimize.minimize`), Raoult's/Henry's law mixtures with colligative properties, and Margules/Wilson/NRTL/UNIQUAC activity models for non-ideal vapor-liquid equilibrium and azeotropes (roots via `scipy.optimize.brentq`).

  ![van der Waals isotherms of CO2, the ethanol-water azeotrope, and water's phase diagram](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_thermo.png)

- [`chemistrykit.solutions`](https://cpoli.github.io/chemistrykit/api/gallery/solutions/) -- solution chemistry: pH/pOH and weak acid/base equilibria with Henderson-Hasselbalch buffers, acid-base titration curves, Ksp solubility equilibria and the common-ion effect, polyprotic and metal-ligand complexation speciation, and Debye-Huckel activity coefficients (equilibria solved with `scipy.optimize.brentq`).

  ![Strong and weak acid titration curves, phosphoric acid speciation, and activity coefficients](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_solutions.png)

- [`chemistrykit.md`](https://cpoli.github.io/chemistrykit/api/gallery/md/) -- molecular dynamics and force fields: the Lennard-Jones fluid in reduced units (periodic boundary conditions, a Verlet neighbor list, pressure, g(r)), Morse/Buckingham/harmonic bonded potentials, velocity-rescaling/Nose-Hoover thermostats, mean-squared displacement and Green-Kubo diffusion, and XYZ trajectory export for VMD/OVITO (numba-compiled pair forces; integrators and thermostats hand-rolled, since they are the subject).

  ![Radial distribution function of liquid argon, MD speeds vs. Maxwell-Boltzmann, and mean-squared displacement](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_md.png)

- [`chemistrykit.statmech`](https://cpoli.github.io/chemistrykit/api/gallery/statmech/) -- statistical mechanics of molecules: translational/rotational/vibrational partition functions and their thermodynamic functions, the Maxwell-Boltzmann speed distribution, Einstein and Debye heat capacities (Debye functions via `scipy.integrate.quad`), Onsager's exact 2D Ising solution (elliptic integrals via `scipy.special`), and a canonical-ensemble lattice-gas adsorption model.

  ![Maxwell-Boltzmann speeds, Onsager's Ising magnetization, and Debye vs. Einstein heat capacity](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_statmech.png)

- [`chemistrykit.quantum`](https://cpoli.github.io/chemistrykit/api/gallery/quantum/) -- quantum chemistry: particle-in-a-box models (with the free-electron model of conjugated-dye color); the quantum harmonic oscillator vs. the exact Morse potential; the rigid rotor; hydrogen-like orbitals; Huckel molecular-orbital theory and its 4n+2 aromaticity rule; a minimal variational treatment of H2+; restricted Hartree-Fock SCF in an STO-3G basis for H2 and HeH+; and Rayleigh-Schrodinger perturbation theory for the anharmonic oscillator (eigenproblems via `scipy.linalg.eigh`, with the SCF loop hand-rolled).

  ![Hydrogen radial distributions, H2+ bonding and antibonding orbitals, and Frost circles](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_quantum.png)

- [`chemistrykit.spectro`](https://cpoli.github.io/chemistrykit/api/gallery/spectro/) -- spectroscopy: the Beer-Lambert absorbance law and its stray-light deviation from linearity; rigid-rotor rotational spectra with isotope shifts; harmonic vs. Morse vibrational band positions plus a Wilson GF-matrix triatomic normal-mode calculation; Franck-Condon vibronic progressions; Gaussian/Lorentzian/Voigt lineshapes (`scipy.special.wofz`); and first-order NMR multiplets, exact second-order (AB, ABX) spectra, and Fourier-transform NMR (`numpy.fft`).

  ![Ethanol's 1H NMR spectrum, a Franck-Condon progression, and the HCl rotational spectrum](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_spectro.png)

- [`chemistrykit.structure`](https://cpoli.github.io/chemistrykit/api/gallery/structure/) -- molecular structure and bonding: a lightweight `Molecule` container; VSEPR geometry prediction with real 3D coordinate generation; point-group determination from 3D coordinates and character tables (hand-rolled); bond order from the Pauling length correlation and Huckel-theory MO coefficients; and formal-charge/oxidation-state assignment from a Lewis structure.

  ![SF4's seesaw geometry in 3D, crystal-field splitting, and Pauling's bond-order/length correlation](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_structure.png)

- [`chemistrykit.electrochem`](https://cpoli.github.io/chemistrykit/api/gallery/electrochem/) -- electrochemistry: the Nernst equation for standard and concentration cells with Debye-Huckel activity corrections; a curated standard-reduction-potential table with redox-couple balancing; Butler-Volmer electrode kinetics and Tafel-plot linearization; galvanic vs. electrolytic cells and Faraday's laws of electrolysis; Cottrell, Ilkovič and Randles-Ševčík voltammetry; and a simplified constant-current battery discharge model with Peukert's-law rate dependence.

  ![Butler-Volmer partial currents, a two-ion polarogram, and a Tafel plot](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_electrochem.png)

- [`chemistrykit.photochem`](https://cpoli.github.io/chemistrykit/api/gallery/photochem/) -- photochemistry: Jablonski-diagram excited-state kinetics built on `chemistrykit.kinetics`'s reaction-network engine; fluorescence/phosphorescence quantum yields and the photochemical quantum yield via Beer-Lambert; Stern-Volmer quenching with a static-vs-dynamic diagnostic; Förster and Dexter energy transfer; and photostationary-state kinetics for a two-state photoswitch.

  ![Jablonski state populations, a Stern-Volmer plot, and FRET efficiency vs. distance](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_photochem.png)

- [`chemistrykit.surface`](https://cpoli.github.io/chemistrykit/api/gallery/surface/) -- surface chemistry and catalysis: Langmuir, Freundlich, BET, Temkin, and Dubinin-Radushkevich adsorption isotherms with their standard linearizations for fitting parameters from data; Langmuir-Hinshelwood single- and dual-site and Eley-Rideal surface-reaction kinetics; temperature-programmed desorption with Redhead analysis (`scipy.integrate.cumulative_trapezoid`); and a turnover-frequency/rate-enhancement catalysis model built on `chemistrykit.kinetics`'s Arrhenius equation.

  ![BET vs. Langmuir isotherms, Balandin's volcano curve, and temperature-programmed desorption](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_surface.png)

- [`chemistrykit.polymer`](https://cpoli.github.io/chemistrykit/api/gallery/polymer/) -- polymer chemistry: ideal random-walk chain statistics and the Flory exponent for real chains under theta/good/poor solvent conditions; molecular-weight-distribution statistics and the closed-form Flory-Schulz distribution; Flory-Huggins mixing thermodynamics and spinodals; step-growth kinetics via the Carothers equation; and chain-growth/free-radical polymerization kinetics built on `chemistrykit.kinetics`'s reaction-network engine.

  ![A freely jointed chain, the Flory-Schulz distribution, and Flory-Huggins spinodals](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_polymer.png)

- [`chemistrykit.crystal`](https://cpoli.github.io/chemistrykit/api/gallery/crystal/) -- crystallography and solid-state chemistry: the 7 crystal systems and general unit-cell volume; hard-sphere packing (packing fraction, coordination number) for SC/BCC/FCC/HCP lattices; ionic-crystal lattice energy via the Born-Lande and Kapustinskii equations, backed by a genuinely converging (Evjen-method, hand-rolled) numerical Madelung constant; Bragg's law and powder-XRD peak positions with structure factors and systematic absences; and Schottky/Frenkel point-defect equilibrium.

  ![Powder XRD patterns of Fe and Cu, Evjen's Madelung sum, and packing fractions](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_crystal.png)

- [`chemistrykit.analytical`](https://cpoli.github.io/chemistrykit/api/gallery/analytical/) -- analytical chemistry: redox and complexometric (EDTA) titration-curve simulation with equivalence-point detection, alongside `chemistrykit.solutions`'s acid-base titrations; chromatographic plate theory and the van Deemter equation (resolution, selectivity); linear-regression calibration curves with IUPAC-convention limits of detection/quantitation; Savitzky-Golay smoothing (`scipy.signal.savgol_filter`); and propagation-of-uncertainty formulas plus Dixon's Q-test and Grubbs' test for outlier rejection (critical values via `scipy.stats`).

  ![Tsvet's chromatogram of leaf pigments, the van Deemter curve, and EDTA titrations](https://raw.githubusercontent.com/cpoli/chemistrykit/main/docs/source/_static/images/readme_analytical.png)

Shared infrastructure, used across the subpackages above rather than
standalone toolkits:

- `chemistrykit.constants` -- chemical constants (R, NA, k_B, Faraday's
  constant, ...) from `scipy.constants`, plus a built-in periodic table
  (all 118 elements, with a chemical-formula parser, e.g.
  `molar_mass("Ca(OH)2")`) and a few well-defined unit conversions.
- `chemistrykit.stoichiometry` -- general equation balancing (including
  ionic and redox equations), limiting reagents, theoretical and percent
  yield, and empirical formulas from percent composition.
- `chemistrykit.integrators` -- shared numerical ODE integrators (RK4,
  leapfrog, Yoshida4, adaptive Dormand-Prince) used across the other
  subpackages.

## Units

Every function takes and returns plain floats and NumPy arrays; there is
no unit-tracking library and nothing checks that units are consistent.
Physical constants in `chemistrykit.constants` are SI values from
`scipy.constants`, and each docstring states the unit of every argument
and return value. Most models work in SI (J, m, Pa, K, mol/m^3). The
exceptions are the conventional ones, documented per function: solution
concentrations in mol/L, NMR shifts in ppm and couplings in Hz,
spectroscopic line positions in cm^-1, the Lennard-Jones fluid in
reduced units (sigma, epsilon, particle mass), and empirical correlations
(Antoine constants, activity-model and rate-constant parameters) in
whatever units their constants were tabulated in. Converting inputs to
the documented units is the caller's job.

## Design

chemistrykit calls `numpy`/`scipy` directly for anything they already
implement (root finding and minimization, quadrature, eigensolvers,
special functions, FFTs, Savitzky-Golay filtering, statistical
distributions), and hand-rolls an algorithm only where no
`numpy`/`scipy` equivalent exists (e.g. the stoichiometry-matrix
reaction-network engine, point-group determination, Madelung sums,
equation balancing) or where the algorithm's own steps are the
pedagogical subject (e.g. the Hartree-Fock SCF loop, MD integrators and
thermostats, Gillespie's stochastic simulation). No dependency on RDKit,
ASE, PySCF, OpenMM, or SMILES/PDB parsing. Performance-critical inner
loops (reaction-network right-hand sides, MD pair forces) are
numba-compiled, and the ODE integrators are shared across domains.

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
every PR (`.github/workflows/ci.yml`) across Python 3.10-3.15 on Linux and
macOS. See [CONTRIBUTING.md](CONTRIBUTING.md) before opening a PR.

### Coverage

```bash
MPLBACKEND=Agg pytest -q --cov=chemistrykit --cov-report=term
```

774 tests, 96% line coverage overall. Per-subpackage coverage:

| Subpackage | Coverage | | Subpackage | Coverage |
|:--|--:|---|:--|--:|
| `analytical` | 99% | | `solutions` | 99% |
| `crystal` | 99% | | `spectro` | 99% |
| `electrochem` | 99% | | `statmech` | 99% |
| `kinetics` | 93% | | `structure` | 98% |
| `md` | 89% | | `surface` | 99% |
| `photochem` | 95% | | `thermo` | 99% |
| `polymer` | 99% | | `integrators` | 47% |
| `quantum` | 98% | | `constants` | 84% |

`visualizers/` modules are smoke-tested only (correct return type/shape,
or that `anim.save()` succeeds) rather than covered line-by-line, per the
testing convention in [CLAUDE.md](CLAUDE.md). `integrators` sits lower
because several of its fixed-step/adaptive methods aren't exercised
directly by its own tests, only indirectly through the subpackages
(`kinetics`, `md`) that call into it; `constants` includes a few
rarely-used unit-conversion helpers not hit by any test.

## Docs

Built docs are hosted at <https://cpoli.github.io/chemistrykit/>, served
from the `gh-pages` branch. To build locally:

```bash
pip install -e ".[docs]"
cd docs && make html
```

See `docs/source/history/` for a chronology of each domain's foundational
breakthroughs, linked to the corresponding implementation at each step.
The README figures are regenerated with `python docs/make_readme_figure.py`
and `python docs/make_readme_subpackage_figures.py`.

## Citation

If you use chemistrykit in your research, please cite it — see
[CITATION.cff](CITATION.cff).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Please note that this project
follows the [Contributor Covenant](CODE_OF_CONDUCT.md).

## License

MIT -- see [LICENSE](LICENSE).
