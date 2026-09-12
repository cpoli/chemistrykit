chemistrykit
============

**chemistrykit** is a unified numerical toolkit for computational
chemistry, sharing common ODE integrators and chemical constants across
domain subpackages. This is an early, in-progress build -- see
``chemistrykit-spec.md`` in the repository root for the full 14-domain
plan; only the domains below exist so far.

- :mod:`chemistrykit.kinetics` -- reaction kinetics: integrated rate laws,
  the Arrhenius equation, Michaelis-Menten enzyme kinetics, a general
  stoichiometric reaction-network engine, and the Brusselator oscillator.
- :mod:`chemistrykit.thermo` -- chemical thermodynamics: equations of
  state (ideal gas, van der Waals, Redlich-Kwong), Clausius-Clapeyron
  phase boundaries, reaction equilibrium (Kp/Kc, van't Hoff, and a
  Gibbs-energy-minimization equilibrium-composition solver), and
  Raoult's/Henry's law mixtures with colligative properties.
- :mod:`chemistrykit.solutions` -- solution chemistry: pH/pOH and weak
  acid/base equilibria with Henderson-Hasselbalch buffers, titration
  curves, Ksp solubility equilibria, and Debye-Huckel activity
  coefficients.
- :mod:`chemistrykit.md` -- molecular dynamics and force fields: the
  Lennard-Jones fluid in reduced units (periodic boundary conditions,
  a Verlet neighbor list, pressure, g(r)), Morse/Buckingham/harmonic
  bonded potentials, and velocity-rescaling/Nose-Hoover thermostats.
- :mod:`chemistrykit.statmech` -- statistical mechanics of molecules:
  translational/rotational/vibrational partition functions and their
  thermodynamic functions, the Maxwell-Boltzmann speed distribution, and
  a canonical-ensemble lattice-gas adsorption model.
- :mod:`chemistrykit.quantum` -- quantum chemistry: particle-in-a-box
  models (with the free-electron model of conjugated-dye color); the
  quantum harmonic oscillator vs. the exact Morse potential; the rigid
  rotor; hydrogen-like orbitals; Huckel molecular-orbital theory and its
  4n+2 aromaticity rule; a minimal variational treatment of H2+; and
  Rayleigh-Schrodinger perturbation theory for the anharmonic oscillator.
- :mod:`chemistrykit.spectro` -- spectroscopy: the Beer-Lambert
  absorbance law and its stray-light deviation from linearity;
  rigid-rotor rotational spectra with isotope shifts; harmonic vs. Morse
  vibrational band positions plus a Wilson GF-matrix triatomic
  normal-mode calculation; Franck-Condon vibronic progressions; and a
  first-order NMR multiplet simulator.
- :mod:`chemistrykit.structure` -- molecular structure and bonding: a
  lightweight ``Molecule`` container; VSEPR geometry prediction with real
  3D coordinate generation; point-group determination from 3D
  coordinates and character tables; bond order from the Pauling length
  correlation and Huckel-theory MO coefficients; and formal-charge/
  oxidation-state assignment from a Lewis structure.

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
   :caption: API reference

   api/kinetics
   api/thermo
   api/solutions
   api/md
   api/statmech
   api/quantum
   api/spectro
   api/structure

.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples/kinetics
   examples/thermo
   examples/solutions
   examples/md
   examples/statmech
   examples/quantum
   examples/spectro
   examples/structure
