chemistrykit
============

.. include:: /_generated/vars.rst

**chemistrykit** is a unified numerical toolkit for computational
chemistry, spanning |num_subpackages| domains -- from reaction kinetics to
crystallography -- sharing common ODE integrators and chemical constants
under one NumPy-based API. All 14 domains from ``chemistrykit-spec.md``
(in the repository root) are implemented, listed below.

Every subpackage is grounded in the chemistry it implements, not just
coded against it: public functions carry runnable, CI-checked examples,
and each subpackage's :doc:`history </history/index>` page traces the
breakthroughs behind it -- from Wilhelmy's 1850 first rate law to the
Brunauer-Emmett-Teller isotherm, whose co-discoverer Paul Emmett went on
to develop uranium-separation barrier materials for the Manhattan
Project -- each one linked to the code that reproduces it.

- :mod:`chemistrykit.analytical` -- analytical chemistry: redox and
  complexometric (EDTA) titration-curve simulation with
  equivalence-point detection, alongside :mod:`chemistrykit.solutions`'s
  acid-base titrations; chromatographic plate theory and the van Deemter
  equation; linear-regression calibration curves with IUPAC-convention
  limits of detection/quantitation; and propagation-of-uncertainty
  formulas plus Dixon's Q-test for outlier rejection.
- :mod:`chemistrykit.crystal` -- crystallography and solid-state
  chemistry: the 7 crystal systems and general unit-cell volume;
  hard-sphere packing (packing fraction, coordination number) for
  SC/BCC/FCC/HCP lattices; ionic-crystal lattice energy via the
  Born-Lande and Kapustinskii equations, backed by a genuinely
  converging (Evjen-method) numerical Madelung constant; Bragg's law and
  powder-XRD peak positions with structure factors and systematic
  absences; and Schottky/Frenkel point-defect equilibrium.
- :mod:`chemistrykit.electrochem` -- electrochemistry: the Nernst
  equation for standard and concentration cells with Debye-Huckel
  activity corrections; a curated standard-reduction-potential table
  with redox-couple balancing; Butler-Volmer electrode kinetics and
  Tafel-plot linearization; galvanic vs. electrolytic cells and
  Faraday's laws of electrolysis; and a simplified constant-current
  battery discharge model with Peukert's-law rate dependence.
- :mod:`chemistrykit.kinetics` -- reaction kinetics: integrated rate laws,
  the Arrhenius equation, Michaelis-Menten enzyme kinetics, a general
  stoichiometric reaction-network engine, and the Brusselator oscillator.
- :mod:`chemistrykit.md` -- molecular dynamics and force fields: the
  Lennard-Jones fluid in reduced units (periodic boundary conditions,
  a Verlet neighbor list, pressure, g(r)), Morse/Buckingham/harmonic
  bonded potentials, and velocity-rescaling/Nose-Hoover thermostats.
- :mod:`chemistrykit.photochem` -- photochemistry: Jablonski-diagram
  excited-state kinetics built on :mod:`chemistrykit.kinetics`'s
  reaction-network engine; fluorescence/phosphorescence quantum yields
  and the photochemical quantum yield via Beer-Lambert; Stern-Volmer
  quenching with a static-vs-dynamic diagnostic; and photostationary-
  state kinetics for a two-state photoswitch.
- :mod:`chemistrykit.polymer` -- polymer chemistry: ideal random-walk
  chain statistics and the Flory exponent for real chains under
  theta/good/poor solvent conditions; molecular-weight-distribution
  statistics and the closed-form Flory-Schulz distribution;
  step-growth kinetics via the Carothers equation; and chain-growth/
  free-radical polymerization kinetics built on
  :mod:`chemistrykit.kinetics`'s reaction-network engine.
- :mod:`chemistrykit.quantum` -- quantum chemistry: particle-in-a-box
  models (with the free-electron model of conjugated-dye color); the
  quantum harmonic oscillator vs. the exact Morse potential; the rigid
  rotor; hydrogen-like orbitals; Huckel molecular-orbital theory and its
  4n+2 aromaticity rule; a minimal variational treatment of H2+; and
  Rayleigh-Schrodinger perturbation theory for the anharmonic oscillator.
- :mod:`chemistrykit.solutions` -- solution chemistry: pH/pOH and weak
  acid/base equilibria with Henderson-Hasselbalch buffers, titration
  curves, Ksp solubility equilibria, and Debye-Huckel activity
  coefficients.
- :mod:`chemistrykit.spectro` -- spectroscopy: the Beer-Lambert
  absorbance law and its stray-light deviation from linearity;
  rigid-rotor rotational spectra with isotope shifts; harmonic vs. Morse
  vibrational band positions plus a Wilson GF-matrix triatomic
  normal-mode calculation; Franck-Condon vibronic progressions; and a
  first-order NMR multiplet simulator.
- :mod:`chemistrykit.statmech` -- statistical mechanics of molecules:
  translational/rotational/vibrational partition functions and their
  thermodynamic functions, the Maxwell-Boltzmann speed distribution, and
  a canonical-ensemble lattice-gas adsorption model.
- :mod:`chemistrykit.structure` -- molecular structure and bonding: a
  lightweight ``Molecule`` container; VSEPR geometry prediction with real
  3D coordinate generation; point-group determination from 3D
  coordinates and character tables; bond order from the Pauling length
  correlation and Huckel-theory MO coefficients; and formal-charge/
  oxidation-state assignment from a Lewis structure.
- :mod:`chemistrykit.surface` -- surface chemistry and catalysis:
  Langmuir, Freundlich, and BET adsorption isotherms with their standard
  linearizations for fitting parameters from data; Langmuir-Hinshelwood
  single- and dual-site surface-reaction kinetics; and a
  turnover-frequency/rate-enhancement catalysis model built on
  :mod:`chemistrykit.kinetics`'s Arrhenius equation.
- :mod:`chemistrykit.thermo` -- chemical thermodynamics: equations of
  state (ideal gas, van der Waals, Redlich-Kwong), Clausius-Clapeyron
  phase boundaries, reaction equilibrium (Kp/Kc, van't Hoff, and a
  Gibbs-energy-minimization equilibrium-composition solver), and
  Raoult's/Henry's law mixtures with colligative properties.

Conventionally imported as ``ck``:

.. code-block:: python

   import chemistrykit as ck

   network = ck.kinetics.StoichiometricNetwork.consecutive(k1=1.0, k2=0.3)
   result = network.integrate((0.0, 10.0), dt=1e-3, method="rk4")

   acid = ck.solutions.WeakAcid(Ca=0.1, Ka=1.8e-5)
   print(acid.pH())

   fluid = ck.md.LJFluid.from_lattice(n_per_side=6, density=0.7, temperature=1.0)
   result = fluid.run(dt=0.002, n_steps=2000)

.. toctree::
   :maxdepth: 2
   :caption: Subpackages
   :hidden:

   subpackages/index

.. toctree::
   :maxdepth: 2
   :caption: History
   :hidden:

   history/index

.. toctree::
   :maxdepth: 2
   :caption: Examples
   :hidden:

   examples/index

.. toctree::
   :maxdepth: 2
   :caption: API
   :hidden:

   api/index
