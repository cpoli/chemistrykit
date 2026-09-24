Breakthroughs in Molecular Dynamics
======================================


.. include:: /_generated/nav/md.rst

.. epigraph::

   "Everything that living things do can be understood in terms of the
   jigglings and wigglings of atoms." -- Richard Feynman, *The Feynman
   Lectures on Physics*, Vol. I, Ch. 1 (1963)

Molecular dynamics is the project of taking that jiggling and wiggling
literally: writing down Newton's second law for every particle in a
system, and simply integrating it forward in time on a computer rather
than solving it in closed form. :mod:`chemistrykit.md` retraces that
project from its first idealized computer experiments through the
specific potentials, integrators, and temperature-control schemes that
turned it into the standard tool of computational chemistry it is today.
This chronology traces that thread, with a pointer to the corresponding
implementation in this package at each stop.

.. contents:: Timeline
   :local:
   :depth: 1

1860 -- 1877 -- Maxwell, Boltzmann, and the Kinetic Theory of Gases
------------------------------------------------------------------------

Decades before anyone could simulate a single molecule's trajectory,
James Clerk Maxwell (1860) and Ludwig Boltzmann (working most
systematically in his 1877 memoir connecting entropy to probability)
worked out what the velocities of a huge population of colliding gas
molecules must look like at thermal equilibrium, purely by statistical
reasoning about an ensemble: the distribution of speeds follows a
specific, universal curve fixed entirely by the temperature and the
molecular mass. In a sense this was the first molecular dynamics result
ever obtained -- without a single simulated trajectory in sight. Boltzmann
went further still, arguing in what is now called the ergodic hypothesis
that a *single* system's time-averaged behavior, followed for long
enough, should reproduce the same statistics as an ensemble average over
many independent copies -- the entire theoretical justification for why
running one simulated trajectory forward in time and averaging over it,
rather than simulating a large ensemble of independent copies, is a
legitimate way to extract equilibrium thermodynamic quantities from
molecular dynamics at all.

*Connection:*
:meth:`~chemistrykit.md.LJFluid.from_lattice` initializes
particle velocities by drawing each Cartesian component independently
from a Gaussian of variance :math:`k_BT/m` -- exactly Maxwell's velocity
distribution -- and
``chemistrykit.md.visualizers.md_plots.plot_speed_distribution()``
checks the distribution of speeds that actually emerges from a *running*
simulated trajectory against the analytic Maxwell-Boltzmann prediction
with no fitting involved: a direct demonstration of the ergodic
equivalence Boltzmann argued for.

*References:* J. C. Maxwell, "Illustrations of the Dynamical Theory of
Gases," Philosophical Magazine 19, 19-32 (1860); L. Boltzmann, "Ueber die
Beziehung zwischen dem zweiten Hauptsatze der mechanischen Waermetheorie
und der Wahrscheinlichkeitsrechnung," Wiener Berichte 76, 373-435 (1877).

.. minigallery:: ../../examples/md/lj_fluid/plot_02_maxwell_boltzmann_check.py

1905 -- Einstein: Diffusion as a Random Walk
------------------------------------------------

Albert Einstein, in the same year as his papers on special relativity
and the photoelectric effect, gave the first quantitative molecular
theory of Brownian motion -- the ceaseless jittering of pollen-sized
particles suspended in a liquid. He argued that each particle is kicked
at random by the surrounding molecules, so that its displacement is a
random walk: it does not travel a distance proportional to the elapsed
time, as a particle in free flight would, but wanders off so that its
*mean-squared* displacement grows linearly with time, with the diffusion
coefficient :math:`D` as the constant of proportionality,

.. math::

   \langle |\vec r(t) - \vec r(0)|^2 \rangle = 2dDt

in :math:`d` dimensions. Jean Perrin's measurements of this spreading
(1908-1909) gave one of the decisive determinations of Avogadro's number
and helped settle the reality of atoms. The same relation became, half
a century later, the standard way of extracting a diffusion coefficient
from a molecular dynamics trajectory: every atom in a simulated liquid
performs exactly the random walk Einstein described, and :math:`D` is
simply read off the long-time slope of its mean-squared displacement.

*Implementation:*
:func:`~chemistrykit.md.systems.transport.unwrap_trajectory` undoes the
periodic wrapping of a sampled trajectory,
:func:`~chemistrykit.md.systems.transport.mean_squared_displacement`
computes the time-origin-averaged mean-squared displacement, and
:func:`~chemistrykit.md.systems.transport.einstein_diffusion_coefficient`
fits Einstein's straight line to it and returns :math:`D`.

*References:* A. Einstein, "Über die von der molekularkinetischen
Theorie der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten
suspendierten Teilchen," Annalen der Physik 17, 549-560 (1905).

.. minigallery:: ../../examples/md/lj_fluid/plot_07_einstein_mean_squared_displacement.py

1924 -- Lennard-Jones and the Intermolecular Potential
------------------------------------------------------------

John Lennard-Jones, studying the equation of state of real (non-ideal)
gases, proposed in 1924 a simple two-parameter empirical form for the
interaction energy between a pair of neutral atoms: a steeply repulsive
term at short range, standing in for the Pauli exclusion principle's
resistance to electron-cloud overlap, and an attractive :math:`r^{-6}`
term at longer range representing the genuine London dispersion force.
The specific repulsive exponent of 12 is, by Lennard-Jones's own later
account, chosen mainly for its convenient algebraic relationship to the
:math:`r^{-6}` attractive term rather than derived from first-principles
quantum theory -- the true short-range repulsion is closer to an
exponential, as R. A. Buckingham's later (1938, below) alternative
potential reflects -- but the resulting form fits noble-gas equations of state well
enough, and is cheap enough to evaluate, that it has remained the default
reference interaction of computational statistical mechanics for a full
century since.

.. math::

   U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]

*Implementation:*
:class:`chemistrykit.md.systems.lj_fluid.LennardJones` implements exactly
this two-term potential, with
:attr:`~chemistrykit.md.LennardJones.r_min` giving the
exact potential-minimum separation; :class:`~chemistrykit.md.systems.lj_fluid.LJFluid`
builds an entire many-body fluid out of nothing but pairwise sums of this
one interaction, evaluated over a periodically rebuilt neighbor list (see
1967, below).

*References:* J. E. Lennard-Jones, "On the Determination of Molecular
Fields," Proc. R. Soc. Lond. A 106, 463-477 (1924); J. E. Lennard-Jones,
"On the Determination of Molecular Fields. II," Proc. R. Soc. Lond. A
106, 709-718 (1924).

.. minigallery:: ../../examples/md/pair_potentials/plot_04_lennard_jones_potential.py

1929 -- Morse and the Anharmonic Diatomic Potential
----------------------------------------------------------

Philip McCord Morse, working on the emerging quantum mechanics of
molecular vibration, proposed in 1929 an exactly solvable potential-energy
curve for a diatomic bond stretch that -- unlike the harmonic oscillator
-- correctly flattens out to a finite dissociation energy at large
separation, and correctly captures the softening (anharmonicity) of real
vibrational spectra as a bond approaches dissociation. Morse's own
motivation was spectroscopic: the exact quantum energy levels of his
potential reproduce the anharmonicity term chemists were already fitting
empirically to observed vibration-rotation spectra, giving that empirical
parameter a genuine physical model behind it for the first time.
Classically -- the sense used throughout this package -- the same
potential is simply the natural anharmonic generalization of a harmonic
bond, agreeing with it exactly at the equilibrium separation's curvature
and diverging increasingly far from it away from equilibrium.

.. math::

   U(r) = D_e\left[1-e^{-a(r-r_e)}\right]^2 - D_e

*Implementation:* :class:`chemistrykit.md.systems.pair_potentials.Morse`
implements exactly this potential, and its
:meth:`~chemistrykit.md.Morse.from_force_constant`
constructor builds one whose curvature at equilibrium matches a given
:class:`~chemistrykit.md.systems.pair_potentials.HarmonicBond` exactly, so
the two can be compared directly as
:class:`~chemistrykit.md.systems.pair_potentials.DiatomicOscillator`
drives increasingly large-amplitude vibrations that pull them apart.

*References:* P. M. Morse, "Diatomic Molecules According to the Wave
Mechanics. II. Vibrational Levels," Phys. Rev. 34, 57-64 (1929).

.. minigallery:: ../../examples/md/pair_potentials/plot_01_morse_vs_harmonic.py

1938 -- Buckingham's Exponential-Repulsion Potential
------------------------------------------------------------

Richard Buckingham revisited the short-range repulsive wall Lennard-Jones
had modeled, for algebraic convenience, as :math:`r^{-12}` and replaced
it with a form closer to what quantum mechanics actually predicts for two
electron clouds forced to overlap: an exponential, since the overlap of
two exponentially decaying atomic wavefunctions falls off exponentially
with separation, not as an inverse power of distance at all.

.. math::

   U(r) = Ae^{-Br} - \frac{C}{r^6}

Pairing this exponential repulsion with the same :math:`r^{-6}`
dispersion attraction Lennard-Jones used keeps the potential's
long-range physics identical while fixing its short-range physics to be
more defensible -- at the cost of a real practical flaw the pure :math:`r^{-12}`
form does not share: because the exponential term eventually loses to
the unbounded :math:`-C/r^6` attraction as :math:`r\to0`, the Buckingham
potential turns over and plunges to :math:`U\to-\infty` at zero
separation instead of diverging to :math:`+\infty` like Lennard-Jones's.
That unphysical inner turnover means a Buckingham-potential simulation
can, in principle, let two particles collapse into each other if they
ever stray inside it, a failure mode Lennard-Jones's purely repulsive
core rules out by construction -- part of why Lennard-Jones's cruder but
better-behaved form, rather than Buckingham's more physically motivated
one, became computational chemistry's default. The exponential form
nonetheless remains the standard choice wherever that inner-turnover risk
is manageable and the extra physical realism is worth the added
parameter and computational cost, most prominently for rare-gas
potentials and the short-range repulsion between ions in a crystal
lattice.

*Implementation:* :class:`chemistrykit.md.systems.pair_potentials.Buckingham`
implements exactly this exponential-plus-dispersion form, as an
alternative to
:class:`~chemistrykit.md.systems.lj_fluid.LennardJones`
sharing the same :math:`r^{-6}` attractive tail but replacing its
:math:`r^{-12}` repulsive term with Buckingham's exponential one.

*References:* R. A. Buckingham, "The Classical Equation of State of
Gaseous Helium, Neon and Argon," Proc. R. Soc. Lond. A 168, 264-283
(1938).

.. minigallery:: ../../examples/md/pair_potentials/plot_03_buckingham_vs_lj.py

1946 -- Hill, Westheimer, and Mayer: The Founding of Molecular Mechanics
-------------------------------------------------------------------------

Two papers appearing side by side in the same 1946 issue of the *Journal
of Chemical Physics* -- one by Terrell Hill, the other by Frank
Westheimer and Joseph Mayer -- independently proposed the same
essentially mechanical picture for computing a molecule's strain energy:
model each bond as a stretched or compressed spring and each bond angle
as a bent hinge, each obeying a simple harmonic restoring law with an
empirically or spectroscopically determined force constant, and sum the
contributions from every bond and angle in the molecule relative to a
strain-free reference geometry. Westheimer and Mayer applied the idea
specifically to explain why certain hindered biphenyl derivatives resist
twisting into their otherwise-preferred planar geometry, an early success
that helped establish the method's credibility. This "spring-and-hinge"
picture -- crude next to a full quantum-mechanical treatment, but vastly
cheaper to compute and, near equilibrium geometries, often good enough --
is the direct ancestor of every classical molecular-mechanics force field
used in computational chemistry since, up to and including the elaborate
multi-term force fields used for protein and nucleic-acid simulation
today.

*Implementation:*
:class:`chemistrykit.md.systems.pair_potentials.HarmonicBond` and
:class:`~chemistrykit.md.systems.pair_potentials.HarmonicAngle` are
exactly Hill, Westheimer, and Mayer's spring and hinge terms, and
:class:`~chemistrykit.md.systems.pair_potentials.HarmonicMolecule`
assembles an arbitrary small molecule's bonds and angles into precisely
the kind of additive strain-energy force field their papers introduced --
built here, as in their derivation, from nothing but Hooke's law applied
to bond lengths and bond angles.

*References:* T. L. Hill, "On Steric Effects," J. Chem. Phys. 14, 465
(1946); F. H. Westheimer and J. E. Mayer, "The Theory of the
Racemization of Optically Active Derivatives of Biphenyl," J. Chem. Phys.
14, 733-738 (1946).

.. minigallery:: ../../examples/md/pair_potentials/plot_02_harmonic_molecule.py

1953 -- Metropolis et al. and the Minimum-Image Convention
------------------------------------------------------------------

Nicholas Metropolis, Arianna and Marshall Rosenbluth, and Augusta and
Edward Teller, working at Los Alamos on one of the first Monte Carlo
calculations ever run on a digital computer (the MANIAC), needed to
simulate a fluid of interacting hard disks without the impossibly
impractical alternative of enclosing it in real, fixed, molecule-scattering
walls that would dominate the results for any computationally affordable
number of particles. Their solution -- adopted since as the default
choice in essentially every particle simulation, molecular dynamics
included -- was to place the particles in a box with periodic boundary
conditions: a particle exiting one face reappears at the opposite face,
and every pairwise interaction is evaluated against whichever periodic
image of a neighbor happens to be closest (the minimum-image convention),
so a modestly sized simulated box behaves, for the interactions that
matter, like a small patch cut out of an effectively infinite bulk fluid
rather than a finite droplet dominated by its own surface.

*Implementation:*
:func:`~chemistrykit.md.wrap_positions` and
:func:`~chemistrykit.md.minimum_image_displacement` implement
exactly these two pieces of the scheme, and every periodic
:class:`chemistrykit.md.systems.lj_fluid.LJFluid` simulation in this
package depends on both: positions are wrapped back into the primary cell
after every step, and every pairwise force and radial-distribution-function
calculation is evaluated under the minimum-image convention via
:func:`~chemistrykit.md.build_neighbor_list`.

*References:* N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H.
Teller, and E. Teller, "Equation of State Calculations by Fast Computing
Machines," J. Chem. Phys. 21, 1087-1092 (1953).

.. minigallery:: ../../examples/md/periodic_boundaries/plot_01_minimum_image_convention.py

1954 -- 1957 -- Green, Kubo, and Transport from Time Correlations
----------------------------------------------------------------------

Einstein's relation extracts a diffusion coefficient from how far
particles wander; Melville Green (1954) and Ryogo Kubo (1957) found a
second, far more general route. Working independently, they showed that
every linear transport coefficient -- diffusion, viscosity, thermal
conductivity, electrical conductivity -- is the time integral of an
*equilibrium* time-correlation function of the corresponding
microscopic flux. The fluctuations a system undergoes spontaneously at
rest already encode how it will respond when pushed out of equilibrium.
For self-diffusion the flux is simply a particle's velocity, and the
Green-Kubo relation reads

.. math::

   D = \frac{1}{d}\int_0^\infty \langle \vec v(0)\cdot\vec v(t)\rangle\,dt ,

exactly equivalent to Einstein's mean-squared-displacement formula but
built from how long a particle "remembers" its velocity. Molecular
dynamics made these correlation functions directly computable: Rahman's
1964 argon simulation (below) found that the velocity autocorrelation
of a dense liquid dips below zero, as an atom rebounds off the cage of
its neighbors, a feature no simple kinetic theory had predicted.

*Implementation:*
:func:`~chemistrykit.md.systems.transport.velocity_autocorrelation`
computes the time-origin-averaged velocity autocorrelation function of
a sampled trajectory, and
:func:`~chemistrykit.md.systems.transport.green_kubo_diffusion_coefficient`
integrates it to give :math:`D`.

*References:* M. S. Green, "Markoff Random Processes and the Statistical
Mechanics of Time-Dependent Phenomena. II. Irreversible Processes in
Fluids," J. Chem. Phys. 22, 398-413 (1954); R. Kubo,
"Statistical-Mechanical Theory of Irreversible Processes. I. General
Theory and Simple Applications to Magnetic and Conduction Problems,"
J. Phys. Soc. Jpn. 12, 570-586 (1957).

.. minigallery:: ../../examples/md/lj_fluid/plot_08_green_kubo_velocity_autocorrelation.py

1957 -- Alder and Wainwright: The First Molecular Dynamics Simulation
-------------------------------------------------------------------------

Berni Alder and Thomas Wainwright, at Lawrence Livermore, used one of the
earliest large scientific computers to do something previously only
attempted with pencil, paper, and drastic approximation: numerically
follow the exact trajectories of several hundred idealized hard-sphere
particles -- elastic billiard balls with no interaction at all except an
instantaneous, momentum-conserving collision whenever their surfaces
touch -- colliding deterministically according to Newton's laws, event by
event, with no continuous force ever evaluated (a hard-sphere "potential"
is either zero or infinite, so there is no smooth force to integrate; the
correct dynamics is instead a sequence of exact collision events). This
is now generally regarded as the first true molecular dynamics
simulation, predating Rahman's use of a continuous, physically realistic
force field (below) by seven years. Alder and Wainwright's most striking,
and initially controversial, result was that this simulated hard-sphere
fluid -- with literally no attractive force anywhere in the model --
nonetheless underwent a sharp, first-order-like transition from a
disordered fluid to an ordered, crystal-like solid as the density
increased, driven by nothing but the particles' mutual exclusion of one
another's volume: an entropic, packing-driven ordering effect, not an
energetic one. It took further analytical and numerical work through the
1960s and beyond, on larger systems than were affordable in 1957, to
confirm that this hard-sphere freezing transition was genuine and not a
finite-size artifact.

*Connection:* this package has no event-driven hard-sphere collision
engine of its own --
:class:`chemistrykit.md.core.base_system.MolecularDynamicsSystem`
integrates smooth, continuous forces via velocity-Verlet, not
instantaneous collision events -- so rather than reproducing Alder and
Wainwright's exact algorithm, the example below approximates a
hard-sphere fluid with a purely repulsive, continuous stand-in built from
:class:`~chemistrykit.md.systems.lj_fluid.LJFluid` (see the
Weeks-Chandler-Andersen entry below), and shows the same qualitative,
packing-driven ordering effect their simulations first revealed, at
low and high density.

*References:* B. J. Alder and T. E. Wainwright, "Phase Transition for a
Hard Sphere System," J. Chem. Phys. 27, 1208-1209 (1957).

.. minigallery:: ../../examples/md/lj_fluid/plot_03_hard_sphere_freezing.py

1964 -- Rahman: The First Realistic Molecular Dynamics Simulation
-----------------------------------------------------------------------

Aneesur Rahman, at Argonne National Laboratory, ran what is generally
regarded as the first "realistic" molecular dynamics simulation: 864
atoms interacting through an actual physically motivated pair potential
(a Lennard-Jones-type form fitted to real argon data) rather than Alder
and Wainwright's idealized hard spheres, integrated forward in time for
long enough to compute genuinely dynamical -- not merely structural --
quantities: self-diffusion coefficients, velocity autocorrelation
functions, and the radial distribution function, compared directly
against real experimental neutron-scattering and diffusion data for
liquid argon. The agreement, close though not perfect, established for
the first time that a simulation built from a simple, generic
interatomic potential could stand in for a genuinely difficult
experiment on a real liquid, launching molecular dynamics as a
predictive tool for real materials rather than a purely idealized-model
exercise.

*Implementation:*
:meth:`~chemistrykit.md.LJFluid.from_lattice`, run at
liquid-argon-like reduced density and temperature, together with
:meth:`~chemistrykit.md.LJFluid.radial_distribution_function`,
reproduce exactly the kind of calculation -- a Lennard-Jones model
liquid's structure, computed from an integrated trajectory rather than
assumed -- that Rahman's paper pioneered; the example below runs it at
Rahman's own argon state point and converts g(r) back to ångströms.

*References:* A. Rahman, "Correlations in the Motion of Atoms in Liquid
Argon," Phys. Rev. 136, A405-A411 (1964).

.. minigallery:: ../../examples/md/lj_fluid/plot_01_rahman_liquid_argon.py

1967 -- Verlet: The Integration Algorithm, the Neighbor List, and g(r)
-------------------------------------------------------------------------

Loup Verlet, extending Rahman's approach to a full computation of the
Lennard-Jones fluid's structure and phase diagram, introduced in 1967 the
numerical integration scheme that bears his name: advance the position
alone using its two most recent values and the currently computed force,
from which the velocity -- needed only for kinetic energy and
temperature, not for advancing the trajectory itself -- can be recovered
afterward as a centered finite difference. Verlet's scheme is
algebraically nothing more than a rearranged Taylor expansion, but it has
a property few comparably simple schemes share: it is symplectic,
exactly conserving (up to a small bounded oscillation, never a systematic
drift) a discretized version of the system's total energy however long
the simulation runs -- decisive for molecular dynamics, where runs of
millions of steps are routine and any per-step energy leak would
otherwise accumulate into nonsense. The same 1967 paper also introduced
the neighbor list that bears Verlet's name: build the list of
interacting pairs out to a slightly enlarged cutoff radius (a "skin"),
and reuse that same list for several subsequent steps, since a pair can
only enter the true interaction range after moving by at least the skin
distance -- turning a force evaluation that would otherwise need
rebuilding from scratch every single step into a far cheaper amortized
calculation.

The specific finite-difference bookkeeping Verlet described in 1967
(positions only) is mathematically equivalent to, but not identical in
implementation to, the more commonly used "velocity Verlet" form --
tracking position and velocity together in a symmetric kick-drift-kick
sequence -- introduced later, and given that name, by Swope, Andersen,
Berens, and Wilson in 1982; both are the same underlying symplectic
algorithm, differing only in which auxiliary quantities are stored
between steps.

*Implementation:*
:meth:`~chemistrykit.md.MolecularDynamicsSystem.step`
advances every system in this subpackage using
``chemistrykit.integrators.velocity_verlet_step()`` -- an alias, as
``chemistrykit.integrators``'s own module docstring notes, for the
same kick-drift-kick leapfrog scheme, in the velocity-Verlet bookkeeping
of Swope et al. (1982) -- and
:class:`chemistrykit.md.utils.neighbor_list.VerletNeighborList` implements
Verlet's own skin-list optimization directly, rebuilding only every
``rebuild_every`` calls rather than at every step.
:meth:`~chemistrykit.md.LJFluid.radial_distribution_function`
computes the same g(r) diagnostic Verlet used to characterize the
Lennard-Jones fluid's structure.

*References:* L. Verlet, "Computer 'Experiments' on Classical Fluids. I.
Thermodynamical Properties of Lennard-Jones Molecules," Phys. Rev. 159,
98-103 (1967); W. C. Swope, H. C. Andersen, P. H. Berens, and K. R.
Wilson, "A computer simulation method for the calculation of equilibrium
constants for the formation of physical clusters of molecules,"
J. Chem. Phys. 76, 637-649 (1982).

.. minigallery:: ../../examples/md/lj_fluid/plot_05_verlet_integrator.py
   ../../examples/md/lj_fluid/plot_06_verlet_neighbor_list.py

1971 -- Weeks, Chandler, and Andersen: The Repulsive-Force Decomposition
--------------------------------------------------------------------------

John Weeks, David Chandler, and Hans Andersen asked a pointed question
about why the Lennard-Jones fluid -- and, by extension, simple liquids
generally -- has the structure it does: is a liquid's characteristic
short-range order set mainly by the attractive part of the potential,
pulling neighboring atoms together, or by the repulsive part, simply
preventing them from overlapping? Their answer, confirmed by comparing
simulations of the full Lennard-Jones potential against simulations of
only its purely repulsive part (the potential truncated and shifted to
be continuous exactly at its own minimum, leaving nothing but the
steeply rising repulsive wall), was decisive: for a dense liquid, the
structure is overwhelmingly set by repulsion and simple packing, with the
attractive tail acting as little more than a smooth, structurally
unimportant background. This "WCA" (Weeks-Chandler-Andersen)
decomposition -- splitting a potential into a purely repulsive reference
system plus a weak perturbation -- became the basis of an entire branch
of liquid-state perturbation theory, and gave Alder and Wainwright's
earlier, seemingly artificial hard-sphere idealization a rigorous
justification: a dense real liquid's structure really is, to good
approximation, a hard-sphere-like packing problem.

*Implementation:* setting
:class:`chemistrykit.md.systems.lj_fluid.LJFluid`'s ``cutoff`` parameter
to exactly
:attr:`~chemistrykit.md.LennardJones.r_min` (:math:`=
2^{1/6}\sigma`) reproduces precisely the WCA purely repulsive reference
potential -- the Lennard-Jones force is repulsive at every separation up
to its own minimum, so truncating there removes the attractive tail
entirely while leaving the potential continuous, courtesy of the
truncated-and-shifted construction already built into
:class:`~chemistrykit.md.systems.lj_fluid.LJFluid`. The example below
simulates this purely repulsive fluid alongside the full Lennard-Jones
fluid at the same dense liquid state point and compares their radial
distribution functions, the comparison at the heart of WCA's argument.

*References:* J. D. Weeks, D. Chandler, and H. C. Andersen, "Role of
Repulsive Forces in Determining the Equilibrium Structure of Simple
Liquids," J. Chem. Phys. 54, 5237-5247 (1971).

.. minigallery:: ../../examples/md/lj_fluid/plot_04_wca_repulsive_reference.py

1977 -- Ryckaert, Ciccotti, and Berendsen: SHAKE and Rigid Bonds
-----------------------------------------------------------------------

The time step of a molecular dynamics run is limited by the fastest
motion in the system, and in a molecule that is almost always the
stretching of bonds to hydrogen atoms -- vibrations roughly ten
femtoseconds long, far faster than the torsions, bends, and diffusive
motions a simulation is usually run to study. Jean-Paul Ryckaert,
Giovanni Ciccotti, and Herman Berendsen removed that bottleneck by
treating such bonds as rigid: each bond length :math:`d_{ij}` becomes a
holonomic constraint :math:`|\vec r_i - \vec r_j| = d_{ij}`, held
exactly by constraint forces directed along the bonds. Their SHAKE
algorithm makes an ordinary, unconstrained Verlet position update and
then corrects it, cycling through the constraints and moving each
bonded pair along its previous bond vector, in inverse proportion to the
atoms' masses, until every bond length is restored to within a
tolerance. Applied first to liquid *n*-alkanes, it let the time step be
set by the slower motions that remain, and constrained bonds (via SHAKE
and its velocity-level companion RATTLE, introduced by Hans Andersen in
1983) became standard in biomolecular simulation.

*Implementation:* :func:`chemistrykit.md.systems.constraints.shake`
implements the SHAKE iteration, and
:class:`~chemistrykit.md.systems.constraints.ShakeMolecule` extends
:class:`~chemistrykit.md.systems.pair_potentials.HarmonicMolecule` with
rigid bonds, applying SHAKE within each velocity-Verlet step followed by
RATTLE's velocity correction.

*References:* J.-P. Ryckaert, G. Ciccotti, and H. J. C. Berendsen,
"Numerical integration of the cartesian equations of motion of a system
with constraints: molecular dynamics of n-alkanes," J. Comput. Phys. 23,
327-341 (1977); H. C. Andersen, "Rattle: A 'velocity' version of the
Shake algorithm for molecular dynamics calculations," J. Comput. Phys.
52, 24-34 (1983).

.. minigallery:: ../../examples/md/constraints/plot_01_shake_rigid_bonds.py

1980 -- 1985 -- Andersen, Nose, and Hoover: Extended-Lagrangian Thermostats
-------------------------------------------------------------------------------

An ordinary velocity-Verlet molecular dynamics run conserves total energy
and therefore samples the microcanonical (NVE) ensemble -- appropriate
for an isolated system, but not for the far more common experimental
situation of a system held at fixed temperature by contact with a much
larger heat bath (the canonical, NVT, ensemble). Hans Andersen, in 1980,
introduced the idea of extending a system's equations of motion with an
additional, fictitious dynamical variable representing coupling to an
external reservoir -- in his original paper, a "piston" variable for
constant-pressure simulation, alongside a separate stochastic scheme for
constant-temperature simulation. Shuichi Nose generalized Andersen's
extended-system idea specifically to temperature control in 1984,
introducing a single additional coordinate and its conjugate momentum
whose equations of motion, derived from an extended Lagrangian, drive the
physical system's time-averaged kinetic temperature to a target value
while still generating genuine canonical-ensemble fluctuations, rather
than simply clamping the temperature the way a cruder rescaling scheme
does. William Hoover reformulated Nose's rather awkward, time-scaled
original equations, within a year, into the far more commonly implemented
"friction coefficient" form -- an extra term :math:`-\xi \vec v_i` added
directly to each particle's equation of motion, with the friction
coefficient :math:`\xi` itself evolving according to how far the
instantaneous kinetic energy sits from its target -- that is now simply
called the Nose-Hoover thermostat.

.. math::

   \dot{\vec v}_i = \frac{\vec F_i}{m_i} - \xi \vec v_i, \qquad
   Q\dot\xi = \left(\sum_i m_i v_i^2\right) - \text{dof}\cdot k_BT

*Implementation:*
:class:`chemistrykit.md.systems.thermostats.NoseHooverThermostat`
implements exactly Hoover's friction-coefficient reformulation, evolving
:math:`\xi` and rescaling velocities by :math:`e^{-\xi\,dt}` once after
every velocity-Verlet step; this simplified, single-step-coupled scheme
is, as its own docstring notes, an approximation to the reference
operator-split ("Trotter factorization") integrator of Martyna,
Tuckerman, Tobias, and Klein (1996) rather than a literal implementation
of it, though it still relaxes the time-averaged kinetic temperature
toward the target exactly as the Nose-Hoover equations require.

*References:* H. C. Andersen, "Molecular dynamics simulations at
constant pressure and/or temperature," J. Chem. Phys. 72, 2384-2393
(1980); S. Nose, "A unified formulation of the constant temperature
molecular dynamics methods," J. Chem. Phys. 81, 511-519 (1984); W. G.
Hoover, "Canonical dynamics: Equilibrium phase-space distributions,"
Phys. Rev. A 31, 1695-1697 (1985).

.. minigallery:: ../../examples/md/thermostats/plot_01_nose_hoover_thermostat.py

1984 -- Berendsen's Weak-Coupling Thermostat
------------------------------------------------

Herman Berendsen and coworkers, the same year as Nose's paper but working
from a very different and far simpler starting intuition, proposed
coupling a simulated system's temperature to a target value the way a
physical thermostat couples a room to a heater: not through an extra
dynamical variable in an extended Lagrangian, but by a direct, first-order
relaxation of the kinetic energy toward its target with a chosen time
constant :math:`\tau`, applied as a velocity rescaling at every step. The
Berendsen ("weak-coupling") thermostat is enormously popular for exactly
this simplicity and its typically fast, robust equilibration behavior,
but it carries a well-documented flaw: because it actively damps
kinetic-energy *fluctuations* rather than merely nudging their average, a
system thermostatted this way does not sample genuine canonical-ensemble
fluctuations at all -- an issue not fully appreciated (and given the name
"flying ice cube" for one of its more dramatic failure modes, spuriously
funneling a system's kinetic energy into overall translation) until well
over a decade after the method's original introduction.

*Implementation:*
:class:`chemistrykit.md.systems.thermostats.BerendsenThermostat`
implements Berendsen's scheme directly, rescaling every velocity after
each step by :math:`\lambda=\sqrt{1+(\Delta t/\tau)(T_0/T-1)}` so that
the temperature relaxes toward :math:`T_0` with time constant
:math:`\tau`. Its :math:`\tau=\Delta t` limit is the immediate, complete
rescaling of
:class:`~chemistrykit.md.systems.thermostats.VelocityRescalingThermostat`,
the crudest member of the same family.

*References:* H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
A. DiNola, and J. R. Haak, "Molecular dynamics with coupling to an
external bath," J. Chem. Phys. 81, 3684-3690 (1984).

.. minigallery:: ../../examples/md/thermostats/plot_02_berendsen_weak_coupling.py

2007 -- Bussi, Donadio, and Parrinello: Stochastic Velocity Rescaling
--------------------------------------------------------------------------

Giovanni Bussi, Davide Donadio, and Michele Parrinello repaired the
flaw in Berendsen's thermostat while keeping what made it popular. In a
canonical ensemble the total kinetic energy :math:`K` of :math:`N_f`
degrees of freedom is not fixed but fluctuates, following a gamma
distribution with mean :math:`\bar K = N_f k_BT_0/2` and relative width
:math:`\sqrt{2/N_f}`. Berendsen's rule relaxes :math:`K` toward its
mean and in doing so squeezes those fluctuations away. Bussi and
coworkers instead let :math:`K` itself obey a stochastic differential
equation: Berendsen's deterministic relaxation plus a random term whose
size is fixed by requiring that the canonical distribution be stationary,

.. math::

   dK = (\bar K - K)\frac{dt}{\tau} + 2\sqrt{\frac{K\bar K}{N_f\tau}}\,dW ,

and then rescale all velocities by a common factor to reach the new
:math:`K`. The equation can be integrated exactly over a time step, the
coupling time :math:`\tau` keeps its intuitive meaning, and since only
a single global factor is applied the dynamics are disturbed very
little. The result samples the canonical ensemble correctly, and it has
since become a standard thermostat in production codes (it is GROMACS's
``v-rescale``).

*Implementation:*
:class:`chemistrykit.md.systems.thermostats.StochasticVelocityRescalingThermostat`
applies the exact finite-step update of Bussi, Donadio, and Parrinello's
appendix, drawing one Gaussian and one :math:`\chi^2_{N_f-1}` random
number per step.

*References:* G. Bussi, D. Donadio, and M. Parrinello, "Canonical
sampling through velocity rescaling," J. Chem. Phys. 126, 014101 (2007).

.. minigallery:: ../../examples/md/thermostats/plot_03_stochastic_velocity_rescaling.py

See Also
--------

- :doc:`/api/md`
- :doc:`/history/solutions_breakthroughs`
