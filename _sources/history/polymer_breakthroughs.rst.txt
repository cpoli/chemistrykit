Breakthroughs in Polymer Chemistry
====================================


.. include:: /_generated/nav/polymer.rst

.. epigraph::

   "The most important result of this work... is the discovery that
   substances of colloidal properties may possess a definite molecular
   structure, and that the size of the molecule, in the true chemical
   sense, may reach or exceed that of the largest colloidal particle." --
   Hermann Staudinger, Nobel Lecture, 1953

Polymer chemistry is the study of molecules built by linking small,
repeating units into chains that can run to hundreds of thousands of
atoms -- long enough that a single molecule's own statistics, not just
its chemical bonds, start to determine the material's bulk properties.
For most of the nineteenth century and the first two decades of the
twentieth, rubber, cellulose, and proteins were widely believed to be
colloidal aggregates of small molecules, held together by some
ill-defined association force rather than ordinary covalent bonds; the
recognition that they are instead genuine, if enormous, single molecules
-- and the working out of exactly how such a chain's size, its
statistical distribution of lengths, and its reaction kinetics behave --
is one of twentieth-century chemistry's most contested and consequential
stories. The systems in :mod:`chemistrykit.polymer` retrace that story,
from the first proposal that macromolecules are real, through the two
distinct polymerization mechanisms (step-growth and chain-growth) and
the chain statistics that describe a polymer's size and shape in
solution. This chronology traces the major conceptual breakthroughs
behind the package, with a pointer to the corresponding implementation
at each stop.

.. contents:: Timeline
   :local:
   :depth: 1

1920 -- Staudinger's Macromolecular Hypothesis
--------------------------------------------------

In a 1920 paper, Hermann Staudinger proposed that rubber, cellulose,
starch, and proteins are not colloidal aggregates of small molecules held
together by some weak, ill-defined association force -- the prevailing
view among organic chemists of the day -- but are instead genuine
single molecules, thousands of atoms long, joined end to end by ordinary
covalent bonds like any other molecule, merely very large ones
("macromolecules," a term Staudinger himself coined). The hypothesis was
met with open hostility from much of the chemistry establishment; senior
colleagues reportedly urged Staudinger to abandon the idea before it
damaged his career, and one is said to have compared it to claiming an
elephant could be found wandering the streets of Zurich. Staudinger spent
the following decade and a half accumulating evidence -- viscosity
measurements, chemical-degradation studies, and eventually X-ray
diffraction of stretched cellulose and rubber showing regular,
molecular-scale repeat spacing -- before the macromolecular picture was
generally accepted. He was awarded the 1953 Nobel Prize in Chemistry
"for his discoveries in the field of macromolecular chemistry."

*Implementation:*
:func:`~chemistrykit.polymer.staudinger_specific_viscosity` implements
Staudinger's viscosity rule, :math:`\eta_\text{sp}/c=K_mM`, the
chain-length dependence behind his viscosity evidence. More broadly,
every model in :mod:`chemistrykit.polymer` presupposes
Staudinger's hypothesis as settled fact rather than a live controversy --
a polymer chain's end-to-end distance
(``chemistrykit.polymer.systems.chain_statistics``), the closed-form
distribution of chain lengths produced by a polymerization
(``molecular_weight_distribution``), and
the degree of polymerization reached by a given reaction mechanism
(``step_growth``,
``chain_growth``) are all quantities
that only make sense once a polymer is understood to be one long, real,
covalently bonded molecule of a *definite* length, not a colloidal
aggregate with no molecular identity of its own.

*References:* H. Staudinger, "Ueber Polymerisation," *Ber. Dtsch. Chem.
Ges.* 53 (1920), 1073-1085; "The Nobel Prize in Chemistry 1953,"
NobelPrize.org; H. Staudinger and W. Heuer, *Ber. Dtsch. Chem.
Ges.* 63 (1930), 222-234 (the viscosity rule illustrated in the gallery
example).

.. minigallery:: ../../examples/polymer/solution_properties/plot_01_staudinger_viscosity_law.py

1929 -- 1936 -- Carothers, Nylon, and the Carothers Equation
-----------------------------------------------------------------

Wallace Carothers, leading a small fundamental-research group at DuPont
that his employer hoped would generate patentable discoveries almost as
an afterthought, set out from 1928 onward to test Staudinger's
macromolecular hypothesis directly, by deliberately *synthesizing* long
chains through simple, well-understood condensation reactions (an
alcohol and an acid forming an ester and releasing water, repeated at
both ends of a growing chain) rather than merely characterizing chains
that occurred naturally. His group's systematic study of these
step-growth ("condensation") polymerizations, beginning with a
foundational 1929 paper laying out the general theory of how such
reactions build up chain length, led directly to the synthesis of the
first nylon fiber (nylon 6,6) in February 1935 -- a fully synthetic
replacement for silk that DuPont announced to the public in 1938 and
brought to market as nylon stockings in 1940, arguably the first
synthetic polymer to become an immediate and famous consumer product.

Carothers's systematic study of these reactions produced a foundational
quantitative result now called the Carothers equation: since any two
molecules with compatible functional groups -- monomers, dimers, or
chains of any length -- can react with each other in an ideal step-growth
polymerization, the number-average degree of polymerization depends only
on the extent of reaction :math:`p` (the fraction of functional groups
that have reacted), not on how long the reaction has run in absolute
time:

.. math::

    \bar{X}_n = \frac{1}{1-p}

This relationship diverges only as :math:`p\to1`, meaning a useful
molecular weight requires driving the reaction to very high (often
>99%) conversion -- a demanding practical requirement that shaped how
step-growth polymers are manufactured industrially ever since, and a
sharp qualitative contrast with the free-radical chain-growth mechanism
below, where high molecular weight appears essentially immediately, even
at low overall monomer conversion.

*Implementation:*
:func:`~chemistrykit.polymer.degree_of_polymerization`
implements exactly the Carothers equation, with
:func:`~chemistrykit.polymer.extent_of_reaction_for_DP`
inverting it to find the conversion needed for a target chain length, and
:func:`~chemistrykit.polymer.degree_of_polymerization_stoichiometric_imbalance`
generalizing it to the case of a stoichiometric imbalance between the two
functional groups (or a deliberately added monofunctional "chain
stopper"), which caps the attainable degree of polymerization even at
complete conversion of the limiting group.

*References:* W. H. Carothers, "Studies on Polymerization and Ring
Formation. I. An Introduction to the General Theory of Condensation
Polymers," *J. Am. Chem. Soc.* 51 (1929), 2548-2559; W. H. Carothers,
"Polymers and Polyfunctionality," *Trans. Faraday Soc.* 32 (1936),
39-49 (the paper the Carothers equation itself is usually cited to).

.. minigallery:: ../../examples/polymer/step_growth/plot_01_carothers.py

1934 -- Kuhn's Random-Walk Model of the Polymer Chain
-----------------------------------------------------------

Werner Kuhn showed that a flexible polymer chain, free to rotate about
each of its backbone bonds and ignoring for the moment the fact that it
cannot pass through itself, behaves statistically exactly like a random
walk of :math:`n` freely jointed segments ("Kuhn segments") of length
:math:`b` -- a chain's real, chemically detailed local structure can be
coarse-grained into an *effective* segment length and count, chosen so
the coarse-grained random walk reproduces the real chain's overall size.
Applying the central limit theorem to this random walk gives, exactly (no
approximation beyond the random-walk idealization itself),

.. math::

    \langle R^2\rangle = nb^2

for the mean-square end-to-end distance -- linear in the number of
segments, the same :math:`\sqrt{n}` scaling of end-to-end distance with
chain length that a random walk in any other physical context obeys.
This ideal-chain result, and the associated random-walk machinery Kuhn
introduced to derive it, became the reference point every more realistic
chain model (accounting for excluded volume, solvent quality, or chain
stiffness) is compared against, exactly the role the free particle plays
in quantum mechanics or the ideal gas in thermodynamics.

*Implementation:*
:class:`chemistrykit.polymer.systems.chain_statistics.IdealChain`
implements exactly this exact random-walk result,
:math:`\langle R^2\rangle = nb^2` and
:math:`\langle R_g^2\rangle = nb^2/6` for the mean-square radius of
gyration, via the shared
:class:`~chemistrykit.polymer.core.base_system.PolymerChainModel`
interface that
:class:`~chemistrykit.polymer.systems.chain_statistics.RealChain` (see
1953, below) also implements, so the two can be compared side by side at
the same chain length and segment size;
:func:`~chemistrykit.polymer.freely_jointed_chain` samples explicit
freely jointed conformations whose averaged :math:`R^2` converges to
:math:`nb^2`.

*References:* W. Kuhn, "Ueber die Gestalt fadenfoermiger Molekuele in
Loesungen," *Kolloid-Zeitschrift* 68 (1934), 2-15 (exact page range as
commonly cited in secondary literature; not independently re-verified
against the original volume).

.. minigallery:: ../../examples/polymer/chain_statistics/plot_03_kuhn_freely_jointed_chain.py

1936 -- Flory and the Molecular-Weight Distribution of Step-Growth Polymers
-------------------------------------------------------------------------------

Paul Flory recognized that Carothers's equation, though correct, gives
only the *average* chain length reached at a given extent of reaction --
a real step-growth polymerization produces a whole distribution of chain
lengths, some far shorter and some far longer than the average, and
knowing that distribution's shape matters for essentially every physical
property of the resulting material. Using elementary probability (each
additional repeat unit in a chain is one more condensation step that did,
or did not, occur), Flory derived the exact closed-form chain-length
distribution for an ideal linear step-growth polymerization: the
probability that a randomly chosen chain has exactly :math:`x` repeat
units follows a geometric distribution, :math:`N_x=(1-p)p^{x-1}`, often
called the "most probable" distribution and, because the closely related
distribution G. V. Schulz derived shortly afterward for chain-growth
systems shares the same functional form, frequently referred to jointly
as the Flory-Schulz distribution. Its number- and weight-average degrees
of polymerization,

.. math::

    \bar{X}_n = \frac{1}{1-p}, \qquad \bar{X}_w = \frac{1+p}{1-p}

reproduce the Carothers equation for :math:`\bar X_n` and give a
polydispersity index :math:`\bar X_w/\bar X_n=1+p`, which approaches --
but, however far the reaction is driven, never exceeds -- exactly 2 as
:math:`p\to1`, one of the most quoted results in all of polymer
chemistry.

*Implementation:*
:func:`~chemistrykit.polymer.flory_schulz_number_fraction`
and
:func:`~chemistrykit.polymer.flory_schulz_weight_fraction`
implement exactly this distribution;
:func:`~chemistrykit.polymer.flory_schulz_number_average_DP`,
:func:`~chemistrykit.polymer.flory_schulz_weight_average_DP`,
and
:func:`~chemistrykit.polymer.flory_schulz_pdi`
give its closed-form moments, cross-checked in
``chemistrykit.polymer.tests.test_molecular_weight_distribution``
against direct numerical summation of the distribution via the shared
``chemistrykit.polymer.utils.moments.number_average`` /
``weight_average`` helpers, which also compute :func:`~chemistrykit.polymer.number_average_molar_mass`
and :func:`~chemistrykit.polymer.weight_average_molar_mass`
directly from any measured (count, molar mass) distribution, Flory-Schulz
or otherwise.

*References:* P. J. Flory, "Molecular Size Distribution in Linear
Condensation Polymers," *J. Am. Chem. Soc.* 58 (1936), 1877-1885.

.. minigallery:: ../../examples/polymer/molecular_weight_distribution/plot_01_flory_schulz.py

1937 -- Flory and the Steady-State Kinetics of Free-Radical Polymerization
-------------------------------------------------------------------------------

Free-radical (chain-growth) polymerization proceeds through three
distinct elementary steps -- initiation of a radical, propagation as
that radical adds one monomer after another, and termination when two
radical chains meet -- and Paul Flory showed how to extract simple,
testable closed-form rate laws from this three-step mechanism using the
*steady-state approximation*: because termination is enormously faster
than initiator decomposition, the total radical concentration relaxes to
a slowly drifting quasi-equilibrium almost immediately, at which point
setting the rate of radical production equal to the rate of radical
consumption gives

.. math::

    [M^\bullet]_{ss} = \sqrt{\frac{fk_d[I]}{k_t}}

Flory's analysis produced the classic, experimentally distinctive
signature of free-radical chain polymerization: the overall rate of
polymerization scales as the *square root* of initiator concentration
(rather than linearly, as a naive mass-action guess might suggest), and
he introduced the *kinetic chain length* :math:`\nu` -- the number of
monomer units added, on average, per radical chain generated -- directly
relating the kinetics to the resulting polymer's degree of
polymerization.

*Implementation:*
:func:`~chemistrykit.polymer.free_radical_network`
builds the full initiation/propagation/termination mechanism as a
:class:`chemistrykit.kinetics.systems.networks.StoichiometricNetwork` and
integrates it numerically (reusing the kinetics domain's general
mass-action reaction-network engine rather than hand-rolling new ODE
machinery);
:func:`~chemistrykit.polymer.steady_state_radical_concentration`
implements exactly Flory's closed-form :math:`[M^\bullet]_{ss}`, used as
a direct cross-check against the numerically integrated network in
``chemistrykit.polymer.tests.test_chain_growth``;
:func:`~chemistrykit.polymer.steady_state_rate_of_polymerization`
and
:func:`~chemistrykit.polymer.kinetic_chain_length`
implement the resulting :math:`R_p\propto\sqrt{[I]}` rate law and the
kinetic chain length itself.

*References:* P. J. Flory, "The Mechanism of Vinyl Polymerization," *J.
Am. Chem. Soc.* 59 (1937), 241-253.

.. minigallery:: ../../examples/polymer/chain_growth/plot_01_free_radical_kinetics.py

1938 -- 1940 -- Mark, Houwink, and the Intrinsic-Viscosity Law
------------------------------------------------------------------

Staudinger's viscosity rule, :math:`\eta_\text{sp}/c=K_mM`, proved too
simple: careful measurements on fractionated samples showed intrinsic
viscosity growing more slowly than linearly with molar mass. Herman Mark
(1938) and Roelof Houwink (1940), and independently Ichiro Sakurada in
Japan, proposed the power law

.. math::

    [\eta] = KM^a

with constants :math:`K` and :math:`a` characteristic of a given polymer,
solvent, and temperature, and :math:`a` typically between 0.5 and 0.8
for flexible coils. Calibrated once against absolute molar masses, the
Mark-Houwink (or Mark-Houwink-Sakurada) equation turned a cheap
capillary-viscometer measurement into a routine molar-mass determination,
still used today (including as the "universal calibration" of
size-exclusion chromatography). Flory and Fox later explained the
exponent: :math:`[\eta]\propto R^3/M` with :math:`R\propto M^\nu` gives
:math:`a=3\nu-1`, so :math:`a=1/2` in a theta solvent and :math:`a=4/5`
with Flory's good-solvent :math:`\nu=3/5`.

*Implementation:*
:func:`~chemistrykit.polymer.mark_houwink_intrinsic_viscosity` implements
the power law, :func:`~chemistrykit.polymer.fit_mark_houwink` recovers
:math:`K` and :math:`a` from a log-log fit of measurements, and
:func:`~chemistrykit.polymer.mark_houwink_exponent_from_flory` gives the
Flory-Fox exponent :math:`a=3\nu-1`.

*References:* H. Mark, in *Der feste Koerper*, ed. R. Saenger (Leipzig:
Hirzel, 1938); R. Houwink, "Zusammenhang zwischen viscosimetrisch und
osmotisch bestimmten Polymerisationsgraden bei Hochpolymeren," *J.
Prakt. Chem.* 157 (1940), 15-18; P. J. Flory and T. G Fox, "Treatment of
Intrinsic Viscosities," *J. Am. Chem. Soc.* 73 (1951), 1904-1908.

.. minigallery:: ../../examples/polymer/solution_properties/plot_02_mark_houwink.py

1941 -- 1942 -- Flory-Huggins Lattice Theory of Polymer Solutions
-----------------------------------------------------------------------

Paul Flory and Maurice Huggins, working independently, applied a simple
lattice model -- imagining the solvent and polymer segments as occupying
sites on a regular lattice, with the polymer's segments constrained to
occupy a connected sequence of neighboring sites -- to derive the first
successful statistical-thermodynamic theory of polymer solutions, giving
a closed-form free energy of mixing in terms of the polymer's degree of
polymerization and a single dimensionless interaction parameter,
conventionally written :math:`\chi`, that captures the net energetic
preference between polymer-solvent and polymer-polymer (or
solvent-solvent) contacts. This single parameter gives, for the first
time, a *quantitative* meaning to what "solvent quality" means physically:
a small or negative :math:`\chi` (polymer-solvent contacts favored, or at
least not disfavored) corresponds to a good solvent that swells the
chain; a large positive :math:`\chi` (polymer-polymer contacts favored)
corresponds to a poor solvent that collapses it; and the special
crossover value at which the two effects exactly balance defines the
theta condition -- the solvent quality at which a real chain's excluded
volume is effectively cancelled out and it behaves, to leading order,
exactly like Kuhn's ideal random walk above.

In lattice units the free energy of mixing per site is

.. math::

    \frac{\Delta F_\text{mix}}{k_BT} = \frac{\phi}{N}\ln\phi
    + (1-\phi)\ln(1-\phi) + \chi\phi(1-\phi)

for polymer volume fraction :math:`\phi` and :math:`N` segments per
chain, with a critical point at :math:`\phi_c=1/(1+\sqrt N)`,
:math:`\chi_c=\tfrac12(1+1/\sqrt N)^2`, which tends to the theta value
:math:`\chi=1/2` for long chains.

*Implementation:*
:func:`~chemistrykit.polymer.flory_huggins_free_energy`,
:func:`~chemistrykit.polymer.flory_huggins_spinodal_chi`, and
:func:`~chemistrykit.polymer.flory_huggins_critical_point` implement the
lattice free energy, its spinodal, and the critical point;
``chemistrykit.polymer.systems.chain_statistics`` uses the same
theta/good/poor solvent-quality vocabulary
(``FLORY_EXPONENTS``, :func:`~chemistrykit.polymer.flory_exponent`) at the
level of the resulting chain-size exponent :math:`\nu`.

*References:* P. J. Flory, "Thermodynamics of High Polymer Solutions,"
*J. Chem. Phys.* 9 (1941), 660; M. L. Huggins, "Solutions of Long Chain
Compounds," *J. Chem. Phys.* 9 (1941), 440; M. L. Huggins, "Theory of
Solutions of High Polymers," *J. Am. Chem. Soc.* 64 (1942), 1712-1719.

.. minigallery:: ../../examples/polymer/solution_properties/plot_03_flory_huggins.py

1941 -- 1943 -- Flory, Stockmayer, and the Theory of Gelation
-----------------------------------------------------------------

When some monomers carry three or more reactive groups, a step-growth
polymerization makes branched molecules, and at a well-defined extent of
reaction the mixture abruptly sets into a gel -- a single molecule
spanning the whole sample, as in the glyptal (glycerol/phthalic
anhydride) resins studied by Roy Kienle. Carothers estimated
the gel point by setting :math:`\bar X_n` to infinity,
:math:`p_c=2/f_\text{avg}`. Paul Flory (1941) instead treated the
branched molecules as trees built by independent reaction events and
found that the network appears when the expected number of further
branches reached from each branch point exceeds one; for
self-condensation of an :math:`f`-functional monomer this gives

.. math::

    p_c = \frac{1}{f-1}, \qquad
    \bar X_w = \frac{1+p}{1-(f-1)p}, \qquad
    \bar X_n = \frac{1}{1-fp/2}

so the *weight*-average size diverges at :math:`p_c` while the number
average stays small (:math:`\bar X_n=4` at the gel point for
:math:`f=3`). Walter Stockmayer (1943) derived the full size
distribution of the branched molecules. The Flory-Stockmayer theory was
an early example of what is now called percolation, and its gel points
bracket experiment from below, with Carothers's from above.

*Implementation:* :func:`~chemistrykit.polymer.flory_stockmayer_gel_point`
and :func:`~chemistrykit.polymer.carothers_gel_point` give the two gel
point estimates;
:func:`~chemistrykit.polymer.branching_number_average_DP` and
:func:`~chemistrykit.polymer.branching_weight_average_DP` give the
pre-gel averages, reducing to the Carothers and Flory-Schulz results for
:math:`f=2`.

*References:* P. J. Flory, "Molecular Size Distribution in Three
Dimensional Polymers. I. Gelation," *J. Am. Chem. Soc.* 63 (1941),
3083-3090; W. H. Stockmayer, "Theory of Molecular Size Distribution and
Gel Formation in Branched-Chain Polymers," *J. Chem. Phys.* 11 (1943),
45-55.

.. minigallery:: ../../examples/polymer/step_growth/plot_02_flory_stockmayer_gelation.py

1944 -- Debye's Light-Scattering Measurement of Weight-Average Molar Mass
-------------------------------------------------------------------------------

Peter Debye showed that the intensity of light scattered by a dilute
polymer solution, extrapolated to zero scattering angle and zero
concentration, gives a direct absolute measurement of the polymer's
molar mass -- and, crucially, gives specifically the *weight*-average
molar mass :math:`M_w`, since larger molecules scatter light in
proportion to the square of their mass while contributing to a
concentration-based average only in proportion to the mass itself. This
mattered because the classical methods already in use (osmotic pressure,
end-group titration) instead measure the *number*-average molar mass
:math:`M_n`, since they count molecules rather than weighing scattered
light -- so the two families of technique are not interchangeable
measurements of "the" molecular weight, but genuinely different moments
of the same underlying distribution, and comparing them (via the
polydispersity index :math:`M_w/M_n`) became, for the first time, a
practical experimental measurement rather than a purely theoretical
distinction.

*Implementation:*
:func:`~chemistrykit.polymer.weight_average_molar_mass`
implements exactly the moment light scattering measures,
:math:`M_w=\sum N_iM_i^2/\sum N_iM_i`, via the shared
``chemistrykit.polymer.utils.moments.weight_average()`` helper;
:func:`~chemistrykit.polymer.number_average_molar_mass`
implements the complementary :math:`M_n` an osmometry- or
end-group-based measurement would instead give, and
:func:`~chemistrykit.polymer.polydispersity_index`
computes their ratio;
:func:`~chemistrykit.polymer.rayleigh_ratio_dilute_mixture`,
:func:`~chemistrykit.polymer.debye_Kc_over_R`, and
:func:`~chemistrykit.polymer.osmotic_pressure_dilute_mixture` model the
two measurements themselves.

*References:* P. Debye, "Light Scattering in Solutions," *J. Appl.
Phys.* 15 (1944), 338-342.

.. minigallery:: ../../examples/polymer/molecular_weight_distribution/plot_02_debye_light_scattering.py

1944 -- Mayo and Lewis: The Copolymer Composition Equation
-------------------------------------------------------------

When two monomers polymerize together, the copolymer that forms usually
has a different composition from the monomer feed, because each growing
radical prefers one monomer over the other. Frank Mayo and Frederick
Lewis (and, in the same year, Turner Alfrey and George Goldfinger)
assumed that a radical's reactivity depends only on its terminal unit,
so there are four propagation rate constants :math:`k_{11}, k_{12},
k_{21}, k_{22}`. Applying a steady state to the two radical types gives
the instantaneous copolymer mole fraction :math:`F_1` in terms of the
feed mole fraction :math:`f_1`:

.. math::

    F_1 = \frac{r_1f_1^2+f_1f_2}{r_1f_1^2+2f_1f_2+r_2f_2^2}, \qquad
    r_1=\frac{k_{11}}{k_{12}},\; r_2=\frac{k_{22}}{k_{21}}

Just two reactivity ratios capture whether a pair copolymerizes
randomly (:math:`r_1r_2=1`), tends to alternate (:math:`r_1,r_2\to0`), or
drifts in composition, and when both ratios are below one the curve
crosses :math:`F_1=f_1` at an azeotropic feed
:math:`f_1^*=(1-r_2)/(2-r_1-r_2)`. Mayo and Lewis tested it on
styrene/methyl methacrylate, and tabulated reactivity ratios became the
standard way to design copolymers.

*Implementation:*
:func:`~chemistrykit.polymer.mayo_lewis_copolymer_composition` implements
the copolymer equation and
:func:`~chemistrykit.polymer.azeotropic_feed_composition` the azeotropic
feed.

*References:* F. R. Mayo and F. M. Lewis, "Copolymerization. I. A Basis
for Comparing the Behavior of Monomers in Copolymerization; The
Copolymerization of Styrene and Methyl Methacrylate," *J. Am. Chem. Soc.*
66 (1944), 1594-1601; T. Alfrey Jr. and G. Goldfinger, "The Mechanism of
Copolymerization," *J. Chem. Phys.* 12 (1944), 205-209.

.. minigallery:: ../../examples/polymer/chain_growth/plot_03_mayo_lewis_copolymerization.py

1949 -- Kratky and Porod's Worm-Like Chain
----------------------------------------------

Kuhn's freely jointed chain suits very flexible polymers, but many
chains -- cellulose derivatives, and later DNA -- bend only gradually.
Otto Kratky and Guenther Porod, analyzing small-angle X-ray scattering
from dissolved chain molecules, modeled such a polymer as a continuous
filament whose direction loses memory of itself exponentially along the
contour, with a characteristic persistence length :math:`P`. For contour
length :math:`L` this gives

.. math::

    \langle R^2\rangle = 2PL\left[1-\frac{P}{L}\left(1-e^{-L/P}\right)\right]

which behaves as a rigid rod (:math:`R^2\approx L^2`) when
:math:`L\ll P` and as an ideal Kuhn coil with Kuhn length :math:`b=2P`
(:math:`R^2\approx2PL`) when :math:`L\gg P`. The worm-like chain is now
the standard model for semiflexible biopolymers; the double-stranded DNA
persistence length of about 50 nm is quoted in these terms.

*Implementation:*
:func:`~chemistrykit.polymer.worm_like_chain_mean_square_end_to_end`
implements the Kratky-Porod formula, with its rod and coil limits checked
in the tests against :math:`L^2` and
:class:`~chemistrykit.polymer.systems.chain_statistics.IdealChain`'s
:math:`nb^2` with :math:`b=2P`.

*References:* O. Kratky and G. Porod, "Roentgenuntersuchung geloester
Fadenmolekuele," *Recl. Trav. Chim. Pays-Bas* 68 (1949), 1106-1122.

.. minigallery:: ../../examples/polymer/chain_statistics/plot_04_kratky_porod_worm_like_chain.py

1953 -- 1963 -- Ziegler, Natta, and Coordination (Insertion) Polymerization
-----------------------------------------------------------------------------

Karl Ziegler, investigating why an aluminum-alkyl-catalyzed ethylene
oligomerization kept unexpectedly stalling at short chain lengths in the
presence of trace nickel contamination, traced the effect to specific
transition-metal impurities and, by 1953, had turned that accidental
observation into something far more consequential: a
titanium-tetrachloride/aluminum-alkyl catalyst system that polymerizes
ethylene at ordinary pressure and temperature into high-density, largely
linear polyethylene -- a sharp departure from the high-pressure,
free-radical process (branched, lower-density polyethylene) that was
until then the only way to make the polymer industrially. Giulio Natta,
learning of Ziegler's catalyst within the year, extended it to propylene
and discovered something Ziegler's own ethylene chemistry had no way to
reveal, since ethylene has no stereocenter to control: the catalyst could
place each propylene monomer into the growing chain with the same
spatial orientation every time, producing "isotactic" polypropylene --
every methyl side-group on the same side of the extended chain -- a
crystalline, mechanically useful stereoregular material, in sharp
contrast to the irregular, low-melting atactic polypropylene a
free-radical mechanism produces.

Mechanistically, both results share a common thread entirely different
from the free-radical chain-growth mechanism above: the growing chain end
stays coordinated to the transition-metal catalyst throughout, and each
new monomer first coordinates to the metal alongside it before
*inserting* into the metal-carbon bond -- "coordination" or "insertion"
polymerization -- so the catalyst's own geometry, not random radical
collision, controls both the resulting polymer's regularity (linear
versus branched) and, for a prochiral monomer like propylene, its
stereochemistry. Ziegler and Natta shared the 1963 Nobel Prize in
Chemistry "for their discoveries in the field of the chemistry and
technology of high polymers."

*Implementation:* the catalyst chemistry itself is outside the scope of
:mod:`chemistrykit.polymer`, but its most visible outcome -- the
stereoregularity of the chain -- is modeled statistically:
:func:`~chemistrykit.polymer.bernoullian_triad_fractions` gives the
meso/racemo triad fractions :math:`[mm]=P_m^2`,
:math:`[mr]=2P_m(1-P_m)`, :math:`[rr]=(1-P_m)^2` for a meso-placement
probability :math:`P_m` (close to 1 for an isotactic Ziegler-Natta
polypropylene, near 1/2 for an atactic free-radical one),
:func:`~chemistrykit.polymer.mean_isotactic_run_length` the average
isotactic run :math:`1/(1-P_m)`, and
:func:`~chemistrykit.polymer.sample_dyad_sequence` explicit sequences.
(Bernoullian tacticity statistics were formalized for NMR analysis by
Bovey and Tiers in 1960.)

*References:* K. Ziegler, E. Holzkamp, H. Breil, and H. Martin, "Das
Mulheimer Normaldruck-Polyaethylen-Verfahren," *Angew. Chem.* 67 (1955),
541-547; G. Natta, P. Pino, P. Corradini, F. Danusso, E. Mantica, G. Mazzanti,
and G. Moraglio, "Crystalline High Polymers of alpha-Olefins," *J. Am.
Chem. Soc.* 77 (1955), 1708-1710; F. A. Bovey and G. V. D. Tiers,
"Polymer NSR Spectroscopy. II. The High Resolution Spectra of Methyl
Methacrylate Polymers Prepared with Free Radical and Anionic
Initiators," *J. Polym. Sci.* 44 (1960), 173-182; "The Nobel Prize in
Chemistry 1963," NobelPrize.org.

.. minigallery:: ../../examples/polymer/stereochemistry/plot_01_ziegler_natta_tacticity.py

1953 -- 1974 -- Flory's Statistical Thermodynamics of Chain Conformations
-------------------------------------------------------------------------------

Paul Flory's 1953 textbook *Principles of Polymer Chemistry* synthesized
two decades of his own and others' work -- step-growth kinetics, the
Flory-Schulz distribution, Flory-Huggins solution theory, and chain
statistics -- into the field's foundational reference, and included his
own mean-field argument for how a real chain's size scales with the
number of segments once excluded volume (the fact that two segments
cannot occupy the same space) is accounted for: rather than the ideal
chain's exact :math:`R\sim n^{1/2}`, a real chain in a good solvent
swells to :math:`R\sim bn^\nu` with a solvent-quality-dependent exponent
:math:`\nu`, which Flory's mean-field balance of the entropic cost of
stretching a chain against the energetic cost of excluded-volume
overlaps estimated at exactly :math:`\nu=3/5` in three dimensions -- a
poor solvent instead collapses the chain to a dense globule of
essentially constant density, :math:`\nu=1/3`, and the theta solvent in
between exactly recovers the ideal chain's :math:`\nu=1/2`. Flory's
approximate mean-field derivation turned out to be remarkably close to
the exact answer (see 1979, below), a fact appreciated only once
renormalization-group methods became available decades later to check
it. Flory's cumulative body of work on the statistical thermodynamics of
macromolecules earned him the 1974 Nobel Prize in Chemistry, "for his
fundamental achievements, both theoretical and experimental, in the
physical chemistry of the macromolecules."

*Implementation:*
:class:`chemistrykit.polymer.systems.chain_statistics.RealChain`
implements exactly this scaling law, :math:`R=bn^\nu`, with
:meth:`~chemistrykit.polymer.RealChain.theta_solvent`,
:meth:`~chemistrykit.polymer.RealChain.good_solvent`,
and
:meth:`~chemistrykit.polymer.RealChain.poor_solvent`
building it with exactly Flory's :math:`\nu=1/2,\,3/5,\,1/3`
(``chemistrykit.polymer.systems.chain_statistics.FLORY_EXPONENTS``);
at the theta point it reduces exactly to
:class:`~chemistrykit.polymer.systems.chain_statistics.IdealChain`'s
end-to-end distance, though (as the module's own docstring flags
explicitly) its radius-of-gyration prefactor is carried over from the
exact ideal-chain result as a scaling-law approximation rather than an
exact result for a real, self-avoiding chain.

*References:* P. J. Flory, *Principles of Polymer Chemistry* (Ithaca,
NY: Cornell University Press, 1953); "The Nobel Prize in Chemistry
1974," NobelPrize.org.

.. minigallery:: ../../examples/polymer/chain_statistics/plot_01_chain_scaling.py

1954 -- Bevington, Melville, and Taylor: Combination vs. Disproportionation
--------------------------------------------------------------------------------

Free-radical chain termination can happen in either of two distinct
ways: two growing radical chains can fuse directly into a single dead
chain ("combination"), or one radical can abstract a hydrogen atom from
the other, leaving two separate dead chains, one with a saturated and
one with an unsaturated chain end ("disproportionation") -- and which
mechanism dominates depends on the specific monomer and temperature.
John Bevington, Harry Melville, and Reginald Taylor worked out how to
distinguish the two experimentally, principally through end-group
analysis (using radioactively labeled initiator to count how many
initiator fragments end up per dead chain -- one per chain for
disproportionation, but only one for every *two* chains for
combination) and through the resulting degree of polymerization's
relationship to the kinetic chain length :math:`\nu` Flory had defined
above: :math:`\bar X_n=2\nu` if termination is by combination (each dead
chain carries two radicals' worth of added monomer), but
:math:`\bar X_n=\nu` if by disproportionation (each dead chain carries
only one radical's worth).

*Implementation:*
:func:`~chemistrykit.polymer.free_radical_network`'s
`mode` parameter selects between exactly these two termination
stoichiometries -- one dead chain per termination event for
``mode="combination"``, two for ``mode="disproportionation"`` -- while
leaving the initiator, monomer, and radical dynamics themselves
identical between the two (only the bookkeeping of how many dead chains
result from a given number of termination events differs), directly
reproducing the factor-of-two relationship between :math:`\bar X_n` and
:math:`\nu` Bevington, Melville, and Taylor's analysis predicts.

*References:* J. C. Bevington, H. W. Melville, and R. P. Taylor, "The
Termination Reaction in Radical Polymerizations," *J. Polym. Sci.* 12
(1954), 449-459, and its sequel in the same volume (exact page ranges
as commonly cited in secondary literature; not independently
re-verified against the original volume).

.. minigallery:: ../../examples/polymer/chain_growth/plot_02_combination_vs_disproportionation.py

1956 -- Szwarc and Living Polymerization
--------------------------------------------

Michael Szwarc found that styrene polymerized in tetrahydrofuran by
sodium naphthalenide, an electron-transfer initiator, has no termination
step: when the monomer runs out, the chain ends stay active -- the
solution keeps its red color -- and adding more monomer, even a different
one, restarts growth. He called these "living" polymers. With all chains
started at once and growing until the monomer is used up, each chain
receives monomers as independent random events, so the chain lengths
follow a Poisson distribution, as Flory had worked out in 1940 for
ethylene oxide polymerization. With :math:`\nu` monomers per initiator,

.. math::

    N_x = \frac{e^{-\nu}\nu^{x-1}}{(x-1)!}, \qquad
    \frac{\bar X_w}{\bar X_n} = 1+\frac{\nu}{(1+\nu)^2}

so the polydispersity index approaches 1 for long chains -- far narrower
than the value of about 2 of step-growth or conventional free-radical
polymerization. Living anionic polymerization made well-defined block
copolymers and near-monodisperse molar-mass standards possible, and it
inspired the later "controlled" radical methods.

*Implementation:*
:func:`~chemistrykit.polymer.poisson_number_fraction`,
:func:`~chemistrykit.polymer.poisson_number_average_DP`,
:func:`~chemistrykit.polymer.poisson_weight_average_DP`, and
:func:`~chemistrykit.polymer.poisson_pdi` implement the distribution and
its moments, and
:func:`~chemistrykit.polymer.simulate_living_polymerization` grows living
chains by random monomer addition to check them.

*References:* M. Szwarc, "'Living' Polymers," *Nature* 178 (1956),
1168-1169; M. Szwarc, M. Levy, and R. Milkovich, "Polymerization
Initiated by Electron Transfer to Monomer. A New Method of Formation of
Block Polymers," *J. Am. Chem. Soc.* 78 (1956), 2656-2657; P. J. Flory,
"Molecular Size Distribution in Ethylene Oxide Polymers," *J. Am. Chem.
Soc.* 62 (1940), 1561-1565.

.. minigallery:: ../../examples/polymer/chain_growth/plot_04_szwarc_living_polymerization.py

1979 -- de Gennes and the Renormalization-Group Refinement of the Flory Exponent
--------------------------------------------------------------------------------------

Pierre-Gilles de Gennes showed, beginning with a 1972 paper mapping the
self-avoiding-walk problem onto the :math:`n\to0` limit of an
:math:`n`-component magnetic spin model, that the full machinery of
critical-phenomena renormalization-group theory -- developed to compute
precise, universal exponents for phase transitions -- could be carried
over essentially unchanged to compute the Flory exponent :math:`\nu`
exactly, rather than by Flory's original mean-field estimate. The
resulting best modern value, :math:`\nu\approx0.588` in three dimensions
(refined through the 1970s by de Gennes and, independently, precise
renormalization-group calculations from Jean Zinn-Justin and
collaborators), is remarkably close to -- but not exactly equal to --
Flory's original :math:`\nu=3/5=0.6`, an accuracy his simple mean-field
argument had no obvious right to. De Gennes synthesized this and the
broader analogy between polymer statistics and critical phenomena in his
1979 book *Scaling Concepts in Polymer Physics*, and was awarded the
1991 Nobel Prize in Physics "for discovering that methods developed for
studying order phenomena in simple systems can be generalized to more
complex forms of matter, in particular to liquid crystals and polymers."

*Connection:*
:class:`~chemistrykit.polymer.systems.chain_statistics.RealChain`'s
default good-solvent exponent
(:meth:`~chemistrykit.polymer.RealChain.good_solvent`,
:math:`\nu=3/5`, from ``FLORY_EXPONENTS``) is Flory's original classic
value; the more precise renormalization-group estimate above is
available separately via
:meth:`~chemistrykit.polymer.RealChain.good_solvent_renormalization_group`
(``FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP``,
:math:`\nu\approx0.588`) for direct comparison. The two *exponents*
differ by only about 2% -- small enough that
``chemistrykit.polymer.systems.chain_statistics``'s qualitative point (a
good solvent swells a chain measurably beyond the ideal, theta-solvent
scaling) is unaffected by which value is used -- but because chain size
is a power of `n`, that 2% exponent gap compounds with chain length: the
two models' predicted sizes diverge further apart the longer the chain,
as the gallery example below shows directly.

*References:* P.-G. de Gennes, "Exponents for the Excluded Volume
Problem as Derived by the Wilson Method," *Phys. Lett. A* 38 (1972),
339-340; P.-G. de Gennes, *Scaling Concepts in Polymer Physics* (Ithaca,
NY: Cornell University Press, 1979); "The Nobel Prize in Physics 1991,"
NobelPrize.org.

.. minigallery:: ../../examples/polymer/chain_statistics/plot_02_flory_vs_de_gennes.py

See Also
--------

- :doc:`/api/polymer`
- :doc:`/history/surface_breakthroughs`
