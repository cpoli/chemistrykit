# Roadmap

What chemistrykit does not do yet, in priority order. chemistrykit is a
teaching toolkit. So the measure is not feature parity with research codes
(PySCF, Cantera, OpenMM). The measure is whether each domain covers the
canonical topics of a physical-chemistry course in that field, tells their
history, and checks that every example still shows the right chemistry.

Coverage is broad: 14 domains, 231 breakthroughs in `docs/source/history/`,
240 gallery examples, 984 tests. The models that were spot-checked give the
right numbers: RHF/STO-3G gives -1.1167 Eh for H2 and -2.8607 Eh for HeH+,
Evjen's sum gives 1.74757 for NaCl, and the Debye-Hückel limiting law
gives γ± = 0.889 at I = 0.01. What is missing is mostly not numbers that
are wrong. Three things are missing:

- **Nothing checks the gallery.** No workflow runs it (item 1), and 210 of
  the 240 examples assert nothing (item 2).
- **The release path has never run** (item 3).
- **Canonical topics are missing.** The biggest gaps are a stiff solver
  for kinetics, Monte Carlo, Hartree-Fock for any molecule beyond H2 and
  HeH+, the rovibrational spectrum, selection rules from group theory,
  equilibrium constants from partition functions, and simulated cyclic
  voltammetry. These are items 4 to 21.

Each item gives the reason, the current state, an outline of the
approach, and the condition for done. You can take each item in its own
prompt.

Some items edit the same files, so do not run them in parallel:

- Items 1 and 3 both edit `.github/workflows/`.
- Items 7 and 8 both edit `quantum/systems/hartree_fock.py`.
- Items 6 and 11 both edit `structure/`.
- Items 14 and 16 both edit `electrochem/`.

Items that add a breakthrough also edit a shared
`docs/source/history/*_breakthroughs.rst` file. If two sessions touch the
same domain, re-read that file before you commit. Items marked **Needs
the maintainer** have a step in a PyPI or GitHub account that only the
maintainer can do.

The project conventions apply to every item:

- NumPy-style docstrings with doctests.
- A cited source formula.
- Closed-form tests.
- `ArrayLike` annotations.
- For each new breakthrough, a history entry with its own gallery
  example(s). The example is titled after the breakthrough, and no other
  entry links to it.

## 1. Run the gallery in CI

**Why.** The gallery is the only place where the examples run, and the
examples are what the docs promise works. Only `ci.yml` (lint, tests,
doctests, mypy) and `release.yml` exist. No workflow runs `make html`, so
nothing catches a broken example on a PR. gh-pages is pushed by hand. The
last push (2026-10-16, "with sitemap and OpenGraph metadata") was built
from changes to `docs/source/conf.py`, `index.rst`, `pyproject.toml` and
`README.md` that are still uncommitted on `main`. So the live site cannot
be rebuilt from the repository. CONTRIBUTING.md asks contributors to
build the docs locally, and nothing enforces it.

**Now.** `docs/source/sg_execution_times.rst` records 1 s for 240 files.
That means the last build skipped every example, because they were
unchanged. No one has measured a cold build.

**Approach.** Commit the pending sitemap/OpenGraph change first, as its
own commit. Then add `.github/workflows/docs.yml`, based on physicskit's.
Learn from physicskit's 75-minute builds:

- Cache `NUMBA_CACHE_DIR` and the sphinx-gallery output between runs.
- On PRs, run only the changed examples, plus the examples that import a
  changed module. Do this through an env-driven `filename_pattern` in
  `conf.py`.
- Run the full gallery on pushes to `main`, and deploy gh-pages from
  that run.
- Make sure a failing `assert` in an example fails the job.

Update CONTRIBUTING.md and the CLAUDE.md "docs" note to match.

**Done when.**

- `docs.yml` runs on PRs and on `main`.
- A PR that breaks one example fails in under 15 minutes.
- gh-pages is deployed only by the workflow, from a committed tree.

## 2. Give every example a `Check` cell

**Why.** physicskit ends each example with a `# Check` cell that asserts
the result against theory. Without one, an example can plot nonsense and
the docs build still passes. 210 of the 240 chemistrykit examples assert
nothing:

| domain | examples with an assert |
|:--|--:|
| analytical, electrochem, kinetics, md, photochem, quantum, solutions, spectro, statmech | 0 of 155 |
| polymer | 1 of 16 |
| surface | 1 of 15 |
| thermo | 1 of 21 |
| crystal | 14 of 17 |
| structure | 13 of 16 |

**Approach.** Do one domain per prompt. In each example, add a final cell
in physicskit's form (`# %%`, `# Check`, `# -----`, one comment line that
states the identity, then `assert` / `np.testing.assert_allclose`). It
should check what the example claims to show, against a closed form or a
textbook value. Some examples:

- The Brusselator oscillates only above `b = 1 + a^2`.
- [B] peaks at `ln(k1/k2)/(k1-k2)`.
- The MD speed histogram matches Maxwell-Boltzmann.
- The H2 SCF energy is -1.1167 Eh.

Use tolerances that pass on CI and still fail if the physics breaks. If
an example makes no checkable claim, rewrite it until it does.

**Done when.** Every `examples/**/plot_*.py` ends with a `Check` cell
that has at least one assertion on the chemistry. The full gallery build
passes (item 1).

## 3. Make the release workflow publish  — **Needs the maintainer**

**Why.** `release.yml` was added after v0.2.1, and no tag has triggered
it since. So the PyPI trusted publisher has never been tested. physicskit
added the same workflow, and its first run failed with
`invalid-publisher`. The version history also has a gap: tag `v0.2.0`
exists, but PyPI only has 0.1.0 and 0.2.1. conda-forge has 0.2.1.

**Now.** The GitHub environment `pypi` exists. Whether PyPI has a
matching trusted publisher (`cpoli/chemistrykit`, `release.yml`,
environment `pypi`) is unknown. `release.yml` does not rebuild gh-pages.

**Approach.**

1. The maintainer confirms or adds the trusted publisher on pypi.org.
2. Add a `workflow_dispatch` dry run that publishes to TestPyPI, or a
   pre-release tag.
3. After item 1 lands, add a job that rebuilds gh-pages after a
   successful publish.
4. Add a note to CHANGELOG.md that 0.2.0 was never published to PyPI.

**Done when.** A pushed tag publishes to PyPI and creates the GitHub
Release, and Zenodo mints a version DOI, all without manual steps.

## 4. A stiff solver for reaction networks

**Why.** Real chemical kinetics is stiff: radical chains, combustion,
atmospheric mechanisms, the Oregonator. Every integrator in
`chemistrykit/integrators/` is explicit (`rk4`, `dopri5`, leapfrog,
Yoshida). Three docstrings recommend `dopri5` for "stiff-ish" networks:
`integrators/__init__.py:18`, `integrators/adaptive.py:14` and
`kinetics/core/base_system.py:224`. That advice is wrong, because an
explicit method's step size is set by stability, not by accuracy.
`examples/kinetics/oscillators/plot_03_oregonator_bz_reaction.py` needs
`max_steps=2_000_000` to finish.

**Approach.**

- Add `method="bdf"` and `method="radau"` to `ReactionNetwork.integrate`.
  Use `scipy.integrate.solve_ivp` (the Design rule says to use scipy
  where scipy has the algorithm).
- Pass an analytic Jacobian, built from the stoichiometry matrix and
  mass-action rate laws.
- Correct the three docstrings.
- Switch the Oregonator example to the stiff solver.
- Add Robertson's (1966) problem as a gallery example and a kinetics
  breakthrough. Gear (1971) is the matching numerics breakthrough.

**Done when.**

- Robertson's problem at t = 40 gives y1 = 0.7158, y2 = 9.185e-6 and
  y3 = 0.2842 (Hairer & Wanner, *Solving ODEs II*).
- The Jacobian is tested against finite differences.
- The Oregonator example runs in under 1 s.
- No docstring recommends an explicit method for stiff problems.

## 5. Metropolis Monte Carlo

**Why.** Monte Carlo is half of molecular simulation, and the package has
none of it. The md history entry "1953 -- Metropolis et al." credits the
paper only for periodic boundaries and the minimum-image convention. It
links an example of those, and leaves out the algorithm the paper is
known for. statmech has Onsager's exact 2D Ising solution but no way to
simulate the model and compare with it.

**Approach.**

- Add a hand-rolled Metropolis 2D Ising model in `statmech`, numba-compiled
  in the same way as `md`, with energy, magnetization and heat capacity.
- Add NVT Metropolis MC of the Lennard-Jones fluid in `md`. Reuse the
  existing pair potential and minimum image. Tune the acceptance ratio,
  and add Widom insertion for the chemical potential.
- Rewrite the 1953 entry around the algorithm, with its own example.
  Keep the minimum-image example under a separate entry or in the
  entry's text, following the gallery rule.

**Done when.**

- The 2D Ising magnetization matches Onsager's to within statistical
  error away from Tc.
- LJ MC and LJ MD at the same state point give the same pressure and
  energy, and the same g(r), within error.
- Detailed balance is tested on a two-state system.

## 6. Crystal-field theory in `structure`

**Why.** The README's structure figure has a panel titled "Crystal-field
splitting in an Oh field". `docs/make_readme_subpackage_figures.py:331`
draws it by hand with `hlines`, and no library code computes it. So the
README shows a feature that does not exist. Crystal-field and
ligand-field theory is also the core of an inorganic course, and nothing
covers it.

**Approach.** Add a `structure/systems/crystal_field.py` with:

- d-orbital splitting for octahedral, tetrahedral and square-planar
  fields;
- high-spin and low-spin filling and the CFSE (including the pairing
  energy);
- spin-only magnetic moments;
- the spectrochemical series as data.

Redraw the README panel from this module. Add history entries for Bethe
(1929) and Van Vleck (1932).

**Done when.**

- CFSE for d0 to d10, in both spin states, matches the standard table
  (for example, low-spin d6 is -2.4 Δo + 2P).
- μ = sqrt(n(n+2)) μB.
- The README figure calls the library.

## 7. Hartree-Fock for real molecules

**Why.** `RestrictedHartreeFock` supports only contracted s-type
Gaussians, so it covers H2 and HeH+ only. Students cannot run H2O, NH3,
CH4 or any first-row molecule. Yet the quantum history already tells the
Boys (1950) and Roothaan-Hall (1951) story.

**Approach.**

- Add p-type (and optionally d-type) Cartesian Gaussian integrals:
  overlap, kinetic, nuclear attraction and ERI. Use McMurchie-Davidson
  or Obara-Saika, with the general Boys function from
  `scipy.special.hyp1f1` or `gammainc`.
- Add STO-3G parameters for H to Ne.
- Accept a `structure.Molecule` as input.
- Add DIIS, Mulliken populations and the dipole moment.
- Keep the SCF loop hand-rolled, because it is the subject being taught.

**Done when.**

- H2O at the Crawford Programming Project #3 geometry gives
  E(RHF/STO-3G) = -74.942080 Eh.
- The H2 and HeH+ results do not change.
- The integrals are tested against the published Crawford integral
  files or closed forms.

## 8. Electron correlation: H2 dissociation

**Why.** The most-taught failure of Hartree-Fock is that the RHF H2
curve dissociates to the wrong limit. Configuration interaction and
unrestricted HF both fix it. The package has no correlated method and no
UHF, so it cannot show the idea of correlation energy at all.

**Approach.** Use the existing s-only STO-3G H2 integrals (this item does
not depend on item 7). Add:

- UHF with symmetry breaking past the Coulson-Fischer point;
- the 2×2 full CI (ground and doubly excited determinants);
- MP2.

Plot RHF, UHF, MP2 and FCI potential curves. Add history entries for
Coulson-Fischer (1949) and Møller-Plesset (1934).

**Done when.**

- FCI at R = 1.4 bohr gives -1.1373 Eh, a correlation energy of
  -0.0206 Eh (Szabo & Ostlund Ch. 4).
- FCI and UHF go to 2 E(H) = -0.9332 Eh (STO-3G) at large R, and RHF
  does not.

## 9. A numerical 1D Schrödinger solver

**Why.** Every quantum model in the package is a closed form: box,
oscillator, Morse, rotor, hydrogen. There is no way to solve an arbitrary
potential. So the package cannot show tunneling, double wells (ammonia
inversion) or the effect of anharmonicity on a fitted PES. "Tunnel"
does not appear anywhere in the code.

**Approach.**

- Add a finite-difference (or Numerov) eigensolver in
  `quantum/systems/`, using `scipy.linalg.eigh_tridiagonal`.
- Add a transmission coefficient for a 1D barrier, based on a transfer
  matrix.
- Add the Wigner tunneling correction to the rate constant in
  `kinetics.rate_theory`.

**Done when.**

- The harmonic and Morse eigenvalues match the closed forms to 1e-6.
- Rectangular-barrier transmission matches the analytic formula.
- A symmetric double well shows the tunneling splitting, and the
  splitting shrinks exponentially with barrier height.

## 10. The rovibrational (IR) spectrum

**Why.** The HCl vibration-rotation spectrum is a standard
physical-chemistry lab. `spectro` has pure rotational lines and pure
vibrational band positions. It has no P and R branches, no
vibration-rotation coupling αe, no band gap at the origin and no
combination differences.

**Approach.**

- Add P- and R-branch line positions with B_v = B_e - α_e(v + ½),
  centrifugal distortion, and Boltzmann intensities with (2J + 1)
  degeneracy.
- Add H35Cl/H37Cl isotope doubling, reusing `isotope_shift_ratio`.
- Add a combination-differences fit that recovers B″ and B′.
- Add a Raman (ΔJ = 0, ±2) counterpart in the same module.

**Done when.**

- The HCl band reproduces the literature constants: ν0 ≈ 2886 cm⁻¹,
  B_e ≈ 10.59 cm⁻¹ and α_e ≈ 0.307 cm⁻¹ (Herzberg).
- The fit recovers the input constants from synthetic lines.
- The intensity maximum sits at J_max = sqrt(kT/2hcB) - ½.

## 11. Selection rules from group theory

**Why.** `structure` finds point groups and stores character tables.
`spectro` computes normal modes. Nothing connects them, although that
link is the reason chemists learn group theory.

**Approach.**

- Reduce a reducible representation into irreps, using the reduction
  formula.
- Build Γ_3N from a `Molecule` (the unshifted-atoms character), then
  subtract translations and rotations to get Γ_vib.
- Label each irrep IR-active (it transforms like x, y or z) and/or
  Raman-active (it transforms like a quadratic). The character tables
  may need their basis-function columns.

**Done when.**

- H2O gives 2A1 + B2, with all three modes IR- and Raman-active.
- CO2 gives Σg⁺ (Raman only) + Πu + Σu⁺ (IR only), which is the mutual
  exclusion rule.
- CH4 gives A1 + E + 2T2, with only T2 IR-active.

## 12. Equilibrium constants from partition functions

**Why.** Computing K from molecular data is the climax of molecular
statistical mechanics. `statmech` has translational,
linear-rotor and harmonic-vibrational partition functions. It has no
electronic partition function, no nonlinear rotor and no link to
`thermo.equilibrium`.

**Approach.** Add:

- an electronic partition function (degeneracies and term energies);
- a nonlinear-rotor partition function with a symmetry number;
- the zero-point energy D0 reference;
- `equilibrium_constant(reactants, products, T)` from the molecular q's.

**Done when.**

- K for H2 + D2 ⇌ 2HD tends to 4 at high T, and agrees with McQuarrie's
  value near 298 K.
- K for Na2 ⇌ 2Na or I2 ⇌ 2I matches a cited textbook value.
- ln K is consistent with the existing Gibbs-energy route.

## 13. Temperature-dependent heat capacities

**Why.** Everywhere in `thermo`, ΔH and ΔS are constant in T. van't Hoff
assumes constant ΔH. Nothing implements Kirchhoff's law. Real-gas
thermochemistry (ΔrH at 700 K for ammonia synthesis, or an adiabatic
flame temperature) cannot be done.

**Approach.**

- Add Shomate or NASA-7 heat-capacity polynomials, with a small cited
  data set from the NIST WebBook (N2, H2, NH3, O2, H2O, CO2, CH4).
- Add closed-form H(T), S(T) and G(T).
- Add Kirchhoff's law, and a van't Hoff equation with ΔCp.
- Add adiabatic flame temperature as a root solve.

**Done when.**

- The closed-form integrals match `quad`.
- ΔrH and K for NH3 synthesis at 700 K are within 1% of JANAF.
- The adiabatic flame temperature of CH4 in air is about 2225 K.

## 14. Simulated cyclic voltammetry

**Why.** Cyclic voltammetry is the central technique of modern
electrochemistry. `electrochem.voltammetry` has only the Randles-Ševčík
peak-current formula. It has no voltammogram, so it cannot show the
duck-shaped curve, peak separation, or the move from reversible to
irreversible behavior.

**Approach.** Add an explicit finite-difference simulation of the
diffusion problem (Bard & Faulkner, App. B; Britz, *Digital
Simulation*). Couple it to the existing Butler-Volmer kinetics at the
electrode and use a triangular potential waveform. The output is a
result dataclass with E, i and t.

**Done when.**

- In the reversible limit, the peak current matches
  `randles_sevcik_peak_current` to 1%, and ΔEp = 57/n mV at 25 °C.
- Nicholson's ψ(ΔEp) relation is reproduced as k0 decreases.

## 15. Pourbaix diagrams

**Why.** Potential-pH diagrams join Nernst and acid-base equilibria. They
are standard in inorganic chemistry, corrosion and geochemistry, and
neither `electrochem` nor `solutions` has them.

**Approach.** Build E-pH predominance boundaries from a list of half-reactions
(including H⁺) and acid-base equilibria. Reuse `nernst_potential` and the
standard-potential table. Add a plotter that also shows the water
stability lines.

**Done when.**

- The water lines are E = 1.229 - 0.0592 pH and E = -0.0592 pH.
- The Fe-H2O diagram's boundaries match a cited diagram (Pourbaix's
  *Atlas*) at the stated activity.

## 16. The electrode interface: double layer and impedance

**Why.** There is no electrical double layer (Gouy-Chapman-Stern) and no
electrochemical impedance spectroscopy. Both are standard
electrochemistry, and the Randles circuit is the model of an electrode
that students use most.

**Approach.**

- Add the Gouy-Chapman potential profile, the Grahame equation, the
  Debye length and the Stern series capacitance.
- Add the Randles-circuit impedance (Rs, Rct, Cdl and the Warburg
  element), with Nyquist and Bode plotters and a complex nonlinear
  least-squares fit (`scipy.optimize.least_squares`).

**Done when.**

- The Debye length of a 1:1 electrolyte at 25 °C is 0.304/√I nm.
- The Nyquist semicircle has diameter Rct and peaks at ω = 1/(Rct·Cdl).
- The fit recovers the input parameters from synthetic noisy data.

## 17. Ewald summation

**Why.** `md` can only simulate neutral Lennard-Jones particles. Ionic
liquids, molten salts and water all need long-range Coulomb sums.
`crystal` already computes the Madelung constant by Evjen's method, so
Ewald summation has a built-in check.

**Approach.** Add a hand-rolled Ewald sum (real space, reciprocal space
and self term), with numba applied to the real-space loop as for LJ
forces. Use it in a charged-particle `md` system.

**Done when.**

- The Ewald energy of the NaCl lattice gives M = 1.747565, the same as
  `crystal.madelung`.
- The result does not depend on the splitting parameter α to 1e-8.
- Forces match finite differences of the energy.

## 18. Marcus theory of electron transfer

**Why.** Marcus theory (Nobel 1992) appears only as a comment in
`photochem/tests/test_electron_transfer.py:20`, which says "no Marcus
inverted region". The inverted region is one of the great predictions of
physical chemistry.

**Approach.** In `photochem.electron_transfer` (or in `kinetics`), add:

- the classical Marcus rate, k = A·exp(-(ΔG° + λ)²/4λkT);
- the four-sphere-cross relation;
- Marcus-Hush-Chidsey for electrodes, as an optional link to
  `electrochem`.

Add a history entry (Marcus 1956) and Miller-Calcaterra-Closs (1984), the
experimental confirmation.

**Done when.**

- The rate is maximal at -ΔG° = λ and symmetric in the
  parabola's variable.
- The cross relation reproduces a cited k12.

## 19. Cooperative ligand binding

**Why.** There is no binding theory beyond Langmuir/Michaelis-Menten:
no Hill equation, no Scatchard plot, no Adair or MWC model.
Hemoglobin's sigmoidal oxygen binding is a standard case in
biophysical chemistry.

**Approach.** Add Hill, Adair and MWC (Monod-Wyman-Changeux) binding
models, with Scatchard and Hill-plot linearizations and a Hill fit. Put
them next to `kinetics.enzyme` or in `solutions`. Add history entries
for Hill (1910) and MWC (1965).

**Done when.**

- n = 1 Hill and MWC with L = 0 reduce to the hyperbolic isotherm.
- The MWC Hill slope at half-saturation matches the closed form.
- A hemoglobin example reproduces P50 ≈ 26 mmHg and n ≈ 2.8.

## 20. EPR hyperfine spectra

**Why.** NMR is well covered. EPR is not covered at all: no g-factor
resonance and no hyperfine multiplets. It is the standard way radicals
are identified.

**Approach.** Generate hyperfine stick spectra from sets of equivalent
nuclei with spin I. Generalize `pascals_triangle_intensities` to
I > ½. Add the resonance field from g, and first-derivative lineshapes
built from `spectro.lineshapes`.

**Done when.**

- The methyl radical gives a 1:3:3:1 quartet.
- The naphthalene radical anion gives 25 lines.
- One I = 1 nucleus gives a 1:1:1 triplet.
- B_res = hν/(g μB).

## 21. Mass-spectrometry isotope patterns

**Why.** `analytical` has no mass spectrometry. Isotope patterns are the
part of MS that is first-principles: a convolution of isotope
abundances. The periodic table already parses formulas.

**Approach.** Add IUPAC isotope masses and abundances as plain data in
`constants`/`periodic_table`. Compute the isotopic fine structure by
polynomial convolution, then bin it to nominal and monoisotopic masses.
Add a stick-spectrum plotter.

**Done when.**

- Cl2 gives 9:6:1, Br2 gives about 1:2:1 and CH2Cl2 gives the textbook
  M:M+2:M+4 ratios.
- The monoisotopic mass of caffeine is 194.0804.
- Abundances sum to 1.

## 22. Tutorials

**Why.** `chemistrykit-spec.md` promises `docs/source/tutorials/*.rst`,
which are longer cross-domain walkthroughs. None exist. physicskit has 20.
The gallery examples are each about one formula. Nothing shows how the
domains combine, and that is the package's main argument for teaching.

**Approach.** Write one tutorial per prompt. Each tutorial imports from
three or more subpackages and has an executable check. Candidates:

- From mechanism to rate law (kinetics: networks, steady state, stiff
  solver).
- A molecule from scratch: VSEPR geometry, point group, normal modes,
  IR-active bands (needs items 10 and 11).
- Liquid argon: LJ potential, MD, g(r), diffusion coefficient, compared
  with experiment.
- The ethanol-water azeotrope, from activity model to distillation.
- From titration curve to buffer design (solutions + analytical).

**Done when.** At least five tutorials are linked from the docs landing
page, and each tutorial executes in the docs build.

## 23. Repository hygiene

**Why.** Small inaccuracies in the docs undermine trust in the parts that
are correct.

- CLAUDE.md says that CI tests Python 3.10-3.12. `ci.yml` tests
  3.10-3.15.
- The README says "774 tests", but there are now 984 test functions.
- `hypothesis` is a test dependency that no test imports.
- `mypy` runs with `check_untyped_defs = false`.

**Approach.** Fix the two docs. Either use `hypothesis` (property tests
for mass conservation in reaction networks, and for equation balancing
in `stoichiometry`) or remove it. Turn on `check_untyped_defs` one
subpackage at a time, or open a tracking note for it.

**Done when.** The docs match CI, and every declared dependency is
used.

---

### Considered and not listed

Each of these is reasonable, but none is canonical enough yet to come
before the items above:

- reaction-diffusion and Turing patterns;
- the Hammett equation and kinetic isotope effects;
- DFT;
- NPT barostats;
- RRKM;
- radiochemistry (much of it is already in `kinetics` and physicskit);
- chemometrics (PCA/PLS);
- standard addition;
- DLVO;
- CSTR/PFR reactors (engineeringkit covers these).
