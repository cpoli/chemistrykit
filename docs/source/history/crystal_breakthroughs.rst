Breakthroughs in Crystallography
=================================


.. include:: /_generated/nav/crystal.rst

.. epigraph::

   "There must be some definite cause why, whenever snow begins to fall,
   its initial formations invariably display the shape of a six-cornered
   starlet." -- Johannes Kepler, *Strena Seu de Nive Sexangula* (*A New
   Year's Gift, or On the Six-Cornered Snowflake*), 1611

Crystallography began as a question about *shape*: why does a mineral,
however it is broken, keep reappearing in the same family of angles and
faces? Two centuries of answers -- first geometric, then atomic, then wave-
mechanical -- turned that question into the modern description of a solid
as a periodic array of points decorated with atoms, probed by diffraction
and quantified by a handful of lattice sums and empirical equations. The
systems in :mod:`chemistrykit.crystal` retrace that arc: from the purely
geometric classification of a unit cell, through the packing and
diffraction consequences of periodicity, to the electrostatic and
statistical-mechanical models of an ionic crystal's energy and defects.
This chronology traces the major conceptual breakthroughs behind the
package, with a pointer to the corresponding implementation at each stop.

.. contents:: Timeline
   :local:
   :depth: 1

1611 -- Kepler's Snowflake and the Close-Packing Conjecture
------------------------------------------------------------

In a short New Year's pamphlet on why snowflakes are six-cornered,
Johannes Kepler asked how equal spheres -- small pellets, or the seeds of a
pomegranate -- can be stacked most tightly. He described the stacking in
which each layer of spheres nests in the hollows of the layer below, so
that every sphere touches twelve others, and asserted that this packing
would be the tightest possible. That arrangement fills a fraction

.. math::

   \frac{\pi}{\sqrt{18}} \approx 0.7405

of space, and it is the face-centered cubic packing (with the hexagonal
close packing sharing the same density by stacking the layers
differently). Kepler's pamphlet is one of the first attempts to explain
a crystal's outward shape by the regular stacking of small identical
units, and his density claim -- the Kepler conjecture -- resisted proof
for almost four centuries, until Thomas Hales's computer-assisted proof
(published 2005).

*Implementation:* :class:`~chemistrykit.crystal.systems.packing.FaceCenteredCubicPacking`
and :class:`~chemistrykit.crystal.systems.packing.HexagonalClosePacking`
both return Kepler's :math:`\pi/\sqrt{18}` from
:meth:`~chemistrykit.crystal.core.base_system.LatticePacking.packing_fraction`,
while :class:`~chemistrykit.crystal.systems.packing.SimpleCubicPacking`
(:math:`\pi/6\approx0.524`) and
:class:`~chemistrykit.crystal.systems.packing.BodyCenteredCubicPacking`
(:math:`\sqrt3\pi/8\approx0.680`) fall below the bound.

*References:* J. Kepler, *Strena Seu de Nive Sexangula* (Frankfurt:
Godfrey Tampach, 1611); T. C. Hales, "A Proof of the Kepler Conjecture,"
Ann. Math. 162 (2005), 1065-1185.

.. minigallery:: ../../examples/crystal/packing/plot_01_packing_fractions.py

1669 -- Steno and the Constancy of Interfacial Angles
------------------------------------------------------

Nicolas Steno (Niels Stensen), a Danish anatomist working in Florence,
compared cross sections of many quartz crystals and found that however
unequally their faces had grown -- some faces large, others shrunk almost
to nothing -- the angles between corresponding faces were always the
same. The external shape of a crystal can vary without limit, but the
angles cannot: they are a fingerprint of the substance. Steno reported
the observation in his 1669 *Prodromus*, and Jean-Baptiste Romé de
l'Isle, measuring hundreds of minerals with Arnould Carangeot's contact
goniometer, generalized it into a law of crystallography in 1783. The law
is the first hint that the faces of a crystal are fixed by an internal
structure rather than by the accidents of growth; Hauy's stacking
hypothesis, next, explained why.

*Implementation:* in modern terms a face is labelled by its Miller
indices :math:`(hkl)`, and the angle between two faces of a cubic crystal
depends on the indices alone,

.. math::

   \cos\phi = \frac{h_1h_2+k_1k_2+l_1l_2}
                  {\sqrt{h_1^2+k_1^2+l_1^2}\sqrt{h_2^2+k_2^2+l_2^2}},

which :func:`~chemistrykit.crystal.interplanar_angle_cubic` evaluates --
for example :math:`54.74°` between a cube face and an octahedron face,
whatever the crystal's size or shape.

*References:* N. Steno, *De solido intra solidum naturaliter contento
dissertationis prodromus* (Florence, 1669); J.-B. L. Romé de l'Isle,
*Cristallographie*, 2nd ed. (Paris, 1783).

.. minigallery:: ../../examples/crystal/crystal_systems/plot_02_steno_interfacial_angles.py

1784 -- 1801 -- Hauy and the Law of Rational Indices
-----------------------------------------------------

Rene Just Hauy, examining a piece of calcite (Iceland spar) that had
shattered when dropped, noticed something a casual glance would miss:
however the crystal broke, every fragment cleaved along the same family
of planes, meeting at the same characteristic angles regardless of the
original specimen's external shape. Hauy generalized the observation into
a structural hypothesis -- a crystal is built from identical, invisibly
small "integrant molecules" stacked in an orderly three-dimensional
array -- and, from that hypothesis, derived the *law of rational indices*:
every natural face of a crystal can be described by intercepts on the
three crystallographic axes that stand in a ratio of small whole numbers.
Hauy first laid out the stacking hypothesis in his 1784 *Essai d'une
théorie sur la structure des crystaux* and gave it its mature, systematic
treatment in the four-volume *Traité de Minéralogie* of 1801, for which he
is often called the founder of modern crystallography. The law of rational
indices is the founding empirical result of the whole field: it is the
first quantitative statement that a crystal's outward geometry is the
visible signature of an internal, discrete periodic order, decades before
anyone had direct evidence -- diffraction or otherwise -- of what such an
order might look like at the atomic scale.

*Implementation:* Hauy's stacking hypothesis, freed of its literal
"integrant molecule" mechanism, is exactly the modern idea that a crystal
reduces to a lattice classified by the equalities and inequalities among
six parameters -- what
:func:`~chemistrykit.crystal.classify_crystal_system`
computes directly, sorting a unit cell into one of the 7 crystal systems
purely from :math:`(a,b,c,\alpha,\beta,\gamma)`, with
:func:`~chemistrykit.crystal.unit_cell_volume`
giving the corresponding general cell-volume formula. Hauy's rational-
index law itself survives essentially unchanged in modern notation as the
integer Miller indices :math:`(h,k,l)` (formalized by Whewell and Miller a
few decades later) that index every reflection and atomic-plane spacing
throughout ``chemistrykit.crystal.systems.xrd`` -- the `hkl` argument
of :func:`~chemistrykit.crystal.structure_factor` and
:func:`~chemistrykit.crystal.d_spacing_cubic` is, in substance,
Hauy's own rational triple of intercepts, and
:func:`~chemistrykit.crystal.miller_indices_from_intercepts` performs the
conversion from a face's rational axial intercepts to that integer triple.

*References:* R. J. Hauy, *Essai d'une théorie sur la structure des
crystaux* (Paris: Gogué & Née de la Rochelle, 1784); R. J. Hauy, *Traité
de Minéralogie*, 4 vols. (Paris: Chez Louis, 1801).

.. minigallery:: ../../examples/crystal/crystal_systems/plot_03_hauy_rational_indices.py

1815 -- Weiss and the Crystal Systems
--------------------------------------

Christian Samuel Weiss, professor of mineralogy in Berlin, reorganized
crystal description around *axes* rather than around Hauy's stacked
molecules. He referred every crystal form to a set of crystallographic
axes and grouped forms into "systems" according to how many of those axes
are equivalent and how they are inclined to one another -- the regular
(cubic) system with three equal perpendicular axes, systems with one
distinct perpendicular axis or with three unequal perpendicular axes, and
the hexagonal system built on a sixfold axis. Friedrich Mohs and Carl
Friedrich Naumann added the oblique (monoclinic and triclinic) systems in
the 1820s, completing the scheme that survives as the seven crystal
systems. Weiss also introduced the axial intercepts ("Weiss parameters")
from which Miller's indices were later derived.

*Implementation:* :func:`~chemistrykit.crystal.classify_crystal_system`
sorts a unit cell into one of the seven systems purely from the
equalities among its axial lengths :math:`a,b,c` and interaxial angles
:math:`\alpha,\beta,\gamma`, exactly Weiss's criterion, and
:func:`~chemistrykit.crystal.unit_cell_volume` gives the cell volume for
any system from the same six parameters.

*References:* C. S. Weiss, "Übersichtliche Darstellung der verschiedenen
natürlichen Abtheilungen der Krystallisationssysteme," Abhandlungen der
Königlichen Akademie der Wissenschaften zu Berlin, physikalische Klasse,
1814-1815 (published 1818).

.. minigallery:: ../../examples/crystal/crystal_systems/plot_01_crystal_systems.py

1848 -- Bravais and the 14 Space Lattices
--------------------------------------------

Auguste Bravais asked a purely combinatorial question about Hauy's
stacking hypothesis: if a crystal is a periodic array of points, how many
geometrically distinct ways can such an array actually be arranged in
three dimensions? Bravais's 1848 memoir answered it exactly: there are 14
distinct lattice types (a correction of Moritz Frankenheim's 1842 count of
15, one of which Bravais showed was a duplicate already present among the
others), classified by the combination of the 7 crystal systems above with
the possible centerings (primitive, body-centered, face-centered,
base-centered) consistent with each system's symmetry. The result gave
crystallography its first rigorous, purely geometric answer to "what kinds
of periodic order are even possible" -- independent of what atoms, if any,
sit at the lattice points -- and remains the organizing classification
underneath every space group and every concrete lattice model in this
package.

*Implementation:* three of the four hard-sphere lattices in
``chemistrykit.crystal.systems.packing`` --
:class:`~chemistrykit.crystal.systems.packing.SimpleCubicPacking`,
:class:`~chemistrykit.crystal.systems.packing.BodyCenteredCubicPacking`,
and :class:`~chemistrykit.crystal.systems.packing.FaceCenteredCubicPacking`
-- are three of Bravais's 14 lattices in the cubic system (P, I, and F
centering respectively), sharing the common
:class:`~chemistrykit.crystal.core.base_system.LatticePacking` interface.
:class:`~chemistrykit.crystal.systems.packing.HexagonalClosePacking` is a
deliberate exception worth flagging precisely because it is easy to
misstate: HCP is *not* itself one of the 14 Bravais lattices -- it is a
2-atom basis (the ABAB stacking sequence) decorating the single primitive
hexagonal Bravais lattice, which is why its coordination number (12)
matches FCC's close packing even though the underlying lattice is
different in kind, not just in centering. The same three cubic Bravais
lattices reappear as the conventional-cell bases (`_CUBIC_BASES`) that
:func:`~chemistrykit.crystal.powder_xrd_peaks` sums over below.
:data:`~chemistrykit.crystal.BRAVAIS_LATTICES` lists all 14 lattices as
the centerings allowed in each crystal system, and
:func:`~chemistrykit.crystal.cubic_lattice_points` generates the P, I,
and F cubic lattices point by point.

*References:* A. Bravais, "Mémoire sur les systèmes formés par des points
distribués régulièrement sur un plan ou dans l'espace," J. Ecole
Polytechnique 19 (1850), 1-128 (the memoir is dated and generally credited
to 1848, the year of its presentation to the Académie des Sciences, though
formal publication followed in 1850).

.. minigallery:: ../../examples/crystal/crystal_systems/plot_04_bravais_cubic_lattices.py

1912 -- von Laue, Friedrich, and Knipping: X-ray Diffraction by Crystals
---------------------------------------------------------------------------

Max von Laue conjectured that if crystals really are periodic arrays of
discrete scattering centers, as Hauy and Bravais's geometry implied, and
if X-rays really are electromagnetic waves of a wavelength comparable to
the spacing between those centers (both open questions at the time), then
a crystal should diffract a beam of X-rays exactly as a ruled grating
diffracts visible light. Walter Friedrich and Paul Knipping tested the
idea in Munich in 1912, directing a beam of X-rays through a crystal and
recording the transmitted beam on a photographic plate placed behind it:
rather than a single spot, the plate showed a symmetric pattern of sharp,
discrete spots surrounding the direct beam -- the unmistakable signature of
interference between many coherently scattering, periodically arranged
centers. The experiment settled both open questions simultaneously: X-rays
are waves capable of interference, and crystals are periodic on a length
scale close to the X-ray wavelength -- turning Hauy and Bravais's purely
geometric lattice hypothesis into a directly, physically observable fact
for the first time. Von Laue received the 1914 Nobel Prize in Physics for
the discovery; Friedrich and Knipping, who ran the actual experiment, did
not share it, a point of enduring controversy in the historical record.

*Implementation:* Laue's own diffraction condition is usually stated in
terms of the scattered-wave momentum transfer matching a reciprocal-
lattice vector, a transmission-geometry formulation this package does not
implement directly; instead, following the Braggs' equivalent and more
widely used reflection-geometry recasting below,
:func:`~chemistrykit.crystal.structure_factor` computes exactly
the kinematic amplitude
:math:`F_{hkl}=\sum_j f_j\exp[2\pi i(hx_j+ky_j+lz_j)]` whose squared
modulus is what a real diffraction experiment -- Laue's transmission
photograph or a modern powder pattern alike -- actually measures, with
:func:`~chemistrykit.crystal.powder_xrd_peaks` assembling full
simulated patterns from it.

*References:* W. Friedrich, P. Knipping, and M. Laue, "Interferenz-
Erscheinungen bei Röntgenstrahlen," Sitzungsberichte der Königlich
Bayerischen Akademie der Wissenschaften, math.-phys. Klasse (1912),
303-322, with the quantitative follow-up "Eine quantitative Prüfung der
Theorie für die Interferenz-Erscheinungen bei Röntgenstrahlen," same
volume, 363-373.

.. minigallery:: ../../examples/crystal/xrd/plot_02_laue_interference.py

1912 -- 1913 -- The Braggs and Bragg's Law
---------------------------------------------

William Henry Bragg and his son William Lawrence Bragg -- 22 years old at
the time -- recast von Laue's diffraction condition into a far more
usable form: treat each family of parallel lattice planes, spaced `d`
apart, as a partially reflecting mirror, and demand that waves reflected
from successive planes arrive back in phase. The resulting condition,

.. math::

   n\lambda = 2d\sin\theta,

is simple enough to solve for `d` (and hence for atomic positions) from a
measured diffraction angle :math:`\theta`, and W. L. Bragg used it,
together with an X-ray spectrometer his father built, to solve the first
crystal structures ever determined this way -- sodium chloride, zinc
blende, and diamond among them -- launching X-ray crystallography as a
practical tool for structure determination rather than a novel
diffraction curiosity. The Braggs shared the 1915 Nobel Prize in Physics
for the work, the only father-son pair to share a Nobel Prize in the same
field to date, with W. L. Bragg -- 25 at the time of the award -- still
the youngest Nobel laureate in the sciences.

*Implementation:* :func:`~chemistrykit.crystal.bragg_angle`
solves exactly the equation above for :math:`\theta`, given a `d`-spacing
and wavelength; :func:`~chemistrykit.crystal.d_spacing_cubic`
supplies `d` for a cubic lattice's :math:`(hkl)` planes, and
:func:`~chemistrykit.crystal.powder_xrd_peaks` combines both
with :func:`~chemistrykit.crystal.structure_factor` -- dropping
reflections with zero structure factor (systematic absences, the modern
tool the Braggs lacked but which falls directly out of the same theory) --
to enumerate an entire simulated powder pattern exactly the way the
Braggs' spectrometer would trace one out peak by peak.

*References:* W. H. Bragg and W. L. Bragg, "The Reflection of X-rays by
Crystals," Proc. R. Soc. Lond. A 88 (1913), 428-438; W. L. Bragg, "The
Diffraction of Short Electromagnetic Waves by a Crystal," Proc. Camb.
Phil. Soc. 17 (1913), 43-57 (read to the Society on 11 November 1912).

.. minigallery:: ../../examples/crystal/xrd/plot_01_powder_xrd.py

1918 -- Scherrer and Crystallite Size from Line Broadening
------------------------------------------------------------

Peter Debye and Paul Scherrer had just shown (1916) that a powder of
randomly oriented crystallites diffracts X-rays into sharp cones, one per
family of lattice planes. Scherrer then asked what happens when the
crystallites are very small -- colloidal gold particles a few nanometres
across. A crystal only a few hundred planes thick cannot cancel slightly
off-Bragg waves completely, so each line is broadened, by an amount
inversely proportional to the crystal's thickness:

.. math::

   \beta = \frac{K\lambda}{\tau\cos\theta},

with :math:`\beta` the line width (in radians of :math:`2\theta`),
:math:`\tau` the crystallite size, and `K` a shape factor near 0.9. The
Scherrer equation remains the standard quick estimate of nanocrystal size
from a powder pattern.

*Implementation:* :func:`~chemistrykit.crystal.scherrer_crystallite_size`
inverts the equation for :math:`\tau` from a measured line width and
position.

*References:* P. Scherrer, "Bestimmung der Grösse und der inneren
Struktur von Kolloidteilchen mittels Röntgenstrahlen," Nachr. Ges. Wiss.
Göttingen, Math.-Phys. Kl. (1918), 98-100; A. L. Patterson, "The Scherrer
Formula for X-Ray Particle Size Determination," Phys. Rev. 56 (1939),
978-982.

.. minigallery:: ../../examples/crystal/xrd/plot_03_scherrer_crystallite_size.py

1918 -- Born and Lande's Lattice-Energy Equation
----------------------------------------------------

Max Born and Alfred Lande set out to explain a specific, measurable
number: the energy released when a mole of gaseous ions condenses into an
ionic crystal, extractable experimentally from a Born-Haber thermochemical
cycle. Their model balances two effects at the equilibrium ion spacing
:math:`r_0` -- the attractive electrostatic energy of the whole lattice of
point charges (bundled into the Madelung constant `M`, the subject of
Madelung's companion paper below) against a short-range repulsion between
overlapping electron clouds, modeled as an empirical inverse power law
:math:`\propto 1/r^n`:

.. math::

   U = -\frac{N_A M |z_+z_-|e^2}{4\pi\varepsilon_0 r_0}\left(1-\frac1n\right)

The Born exponent `n` was fixed independently, from each ion's measured
compressibility, and tabulated by isoelectronic noble-gas configuration --
a scheme still in use essentially unchanged today. It was, together with
Madelung's lattice sum below, the first quantitative theory to connect a
crystal's cohesive energy to nothing but its structure, its ions' charges,
and a single empirically fitted repulsion exponent.

*Implementation:* :class:`chemistrykit.crystal.systems.lattice_energy.BornLandeLatticeEnergy`
evaluates exactly this equation;
``chemistrykit.crystal.utils.reference_data.BORN_EXPONENTS`` reproduces
Born and Lande's own noble-gas-configuration exponent table, and
``average_born_exponent()``
implements the standard prescription (arithmetic mean of the two ions'
exponents) for a salt of two different ion types.

*References:* M. Born and A. Lande, "Über die Berechnung der
Kompressibilität regulärer Kristalle aus der Gittertheorie," Verh. Dtsch.
Phys. Ges. 20 (1918), 210-216.

.. minigallery:: ../../examples/crystal/lattice_energy/plot_01_born_lande.py

1918 -- Madelung's Electrostatic Lattice Sum
------------------------------------------------

Alongside Born and Lande's equation above, Erwin Madelung worked out the
purely geometric half of the same problem: for ions of unit charge sitting
at the points of a given lattice, what is the net electrostatic energy of
one ion in the field of every other ion in the (infinite) crystal? The
answer is a single dimensionless number per structure type -- the Madelung
constant `M` -- defined so that the electrostatic lattice energy per ion
pair is :math:`U=-Mz_+z_-e^2/(4\pi\varepsilon_0 r_0)`, exactly the
quantity Born and Lande's equation needs as an input. Madelung's own 1918
paper is a general treatment of the electric field produced by any regular
array of point charges, not confined to any one structure; specializing it
to the rock-salt (NaCl) structure gives the specific constant
:math:`M_{NaCl}\approx1.7476` used throughout this package. What Madelung's
paper does not resolve -- and what would remain a subtlety for over a
decade -- is that the defining lattice sum is only *conditionally*
convergent, a complication taken up by Evjen below.

*Implementation:* :func:`~chemistrykit.crystal.madelung_constant_nacl`
computes exactly this NaCl Madelung constant (via the genuinely convergent
Evjen-method summation described in the next entry, rather than the naive
sum Madelung's own paper does not distinguish from the correct one), and
``MADELUNG_CONSTANT_NACL_LITERATURE``
records the accepted literature value it reproduces to 5 significant
figures; :class:`chemistrykit.crystal.systems.lattice_energy.BornLandeLatticeEnergy`
consumes it directly as the `madelung_constant` feeding Born and Lande's
equation above.

*References:* E. Madelung, "Das elektrische Feld in Systemen von
regelmäßig angeordneten Punktladungen," Phys. Z. 19 (1918), 524-533 (this
page range is the one most consistently given in secondary crystallography
references; an occasional secondary source instead cites p. 32 of the same
volume, and the discrepancy has not been independently resolved against
the original issue).

.. minigallery:: ../../examples/crystal/madelung/plot_02_madelung_lattice_sum.py

1926 -- Frenkel and the Vacancy-Interstitial Defect
-------------------------------------------------------

Yakov Frenkel, studying the thermal motion of atoms in solids and liquids,
recognized that a perfectly ordered crystal is a `T=0` idealization: at any
finite temperature, thermal fluctuations will occasionally kick an ion out
of its regular lattice site into a normally unoccupied interstitial gap,
leaving a vacancy behind. Because forming such a pair costs a fixed
enthalpy but also raises the crystal's configurational entropy (there are
combinatorially many ways to place a handful of vacancy-interstitial pairs
among many sites), minimizing the resulting Gibbs energy gives a small but
strictly nonzero equilibrium population of these defects at any `T>0` --
not a flaw to be engineered away, but a necessary consequence of
thermodynamics. This vacancy-interstitial pairing, now called a Frenkel
defect, was the first of the two canonical point-defect mechanisms in
ionic solids to be identified; the second, the Schottky defect, followed
four years later.

*Implementation:* :func:`~chemistrykit.crystal.frenkel_defect_concentration`
implements exactly the resulting Boltzmann-factor equilibrium population,
:math:`n_F=\sqrt{NN_i}\exp(-\Delta H_F/2k_BT)`, with `N` the number of
normal lattice sites and :math:`N_i` the number of available interstitial
sites.

*References:* J. Frenkel, "Über die Wärmebewegung in festen und flüssigen
Körpern," Z. Phys. 35 (1926), 652-669.

.. minigallery:: ../../examples/crystal/defects/plot_01_frenkel_defects.py

1926 -- Goldschmidt's Tolerance Factor
---------------------------------------

Victor Moritz Goldschmidt, working in Oslo, set out to predict crystal
structures from the sizes of the ions alone. Using a new set of ionic
radii determined by his group, he formulated "laws of crystal chemistry":
the structure a compound adopts is governed by the numbers, sizes, and
polarizabilities of its ions. His best-known example is the perovskite
structure :math:`ABX_3`, in which a small B cation sits in an octahedron
of X anions and a large A cation fills the twelve-coordinate cavity
between octahedra. In the ideal cubic cell the A-X distance is
:math:`\sqrt2` times the B-X distance, so hard-sphere ions fit exactly
when the tolerance factor

.. math::

   t = \frac{r_A + r_X}{\sqrt{2}\,(r_B + r_X)}

equals 1. Values somewhat below 1 lead to tilted octahedra and lower
symmetry, and values above 1 to hexagonal or polar distortions. The
tolerance factor is still the first screen applied to candidate
perovskites, from ferroelectrics to solar-cell absorbers.

*Implementation:* :func:`~chemistrykit.crystal.goldschmidt_tolerance_factor`
evaluates :math:`t` from three ionic radii (for example Shannon's).

*References:* V. M. Goldschmidt, "Die Gesetze der Krystallochemie,"
Naturwissenschaften 14 (1926), 477-485.

.. minigallery:: ../../examples/crystal/crystal_chemistry/plot_02_goldschmidt_tolerance_factor.py

1929 -- Pauling's Rules and the Radius-Ratio Rule
--------------------------------------------------

Linus Pauling condensed the structural chemistry of ionic crystals into
five rules. The first says that each cation is surrounded by a
coordination polyhedron of anions, with the cation-anion distance set by
the sum of the radii and the coordination number set by the radius ratio
:math:`r_+/r_-`. The geometry is simple: a polyhedron is stable only
while its anions can touch the central cation without overlapping each
other, which gives exact lower limits for each coordination --
:math:`\sqrt{3/2}-1\approx0.225` for a tetrahedron,
:math:`\sqrt2-1\approx0.414` for an octahedron, and
:math:`\sqrt3-1\approx0.732` for a cube. The remaining rules (electrostatic
valence, sharing of edges and faces, and parsimony) explain how the
polyhedra link together. Pauling's rules remain the standard qualitative
guide to ionic and mineral structures, even though real ions are not
hard spheres and the radius-ratio rule fails for a sizeable minority of
salts.

*Implementation:* :func:`~chemistrykit.crystal.radius_ratio_coordination`
returns a :class:`~chemistrykit.crystal.RadiusRatioPrediction` (ratio,
coordination number, polyhedron) using the geometric limits tabulated in
:data:`~chemistrykit.crystal.RADIUS_RATIO_LIMITS`.

*References:* L. Pauling, "The Principles Determining the Structure of
Complex Ionic Crystals," J. Am. Chem. Soc. 51 (1929), 1010-1026.

.. minigallery:: ../../examples/crystal/crystal_chemistry/plot_01_pauling_radius_ratio.py

1930 -- Wagner and Schottky: The Paired-Vacancy Defect
------------------------------------------------------------

Carl Wagner and Walter Schottky, developing a general thermodynamic theory
of "ordered mixed phases" (crystals whose stoichiometry can deviate
slightly from an ideal integer ratio), identified the second canonical
point-defect mechanism: rather than displacing an ion to an interstitial
site, a crystal can instead remove a stoichiometric pair of ions --
one cation and one anion -- from the bulk entirely, depositing them at the
crystal's surface and leaving a paired cation-anion vacancy behind. Because
this mechanism removes ions in a fixed charge-neutral ratio rather than
merely relocating one, it dominates in ionic solids where an isolated
interstitial cation (or anion) would be too large or too costly to
accommodate, and it became the founding concept of what is now called
defect chemistry -- the systematic study of how nonstoichiometry and ionic
conductivity in real crystals are controlled by exactly this kind of
point-defect population, rather than being simply the failure of an
otherwise-perfect lattice.

*Implementation:* :func:`~chemistrykit.crystal.schottky_defect_concentration`
implements the corresponding Boltzmann-factor equilibrium population,
:math:`n_S=N\exp(-\Delta H_S/2k_BT)`, and
``chemistrykit.crystal.visualizers.crystal_plots.plot_defect_concentration_vs_temperature()``
plots it directly alongside the Frenkel population above, on the same
semi-log temperature axis, making the shared Boltzmann-factor functional
form -- and the very different magnitudes two different formation
enthalpies produce -- visually explicit.

*References:* W. Schottky and C. Wagner, "Theorie der geordneten
Mischphasen," Z. Phys. Chem. B 11 (1930), 163-210.

.. minigallery:: ../../examples/crystal/defects/plot_02_schottky_defects.py

1932 -- Evjen's Method for a Genuinely Convergent Madelung Sum
--------------------------------------------------------------------

The Madelung lattice sum that Madelung's 1918 paper defines has an
unpleasant property that went unaddressed for over a decade: it is only
*conditionally* convergent. Summing the alternating :math:`\pm1/r` series
over a naively growing cube or sphere of ions does not settle down to a
single value as the cutoff grows -- it oscillates, because the outermost
shell added at each step of the expansion is not itself electrically
neutral, so its contribution never shrinks to zero the way an
absolutely-convergent series' tail would. H. M. Evjen supplied the fix: give
any ion whose position lies exactly on the truncation region's boundary a
fractional weight (one-half for a face, one-quarter for an edge,
one-eighth for a corner of the truncating cube) chosen precisely so that
every finite partial sum is built from an electrically neutral region,
restoring genuine, rapidly converging convergence to the correct
conditionally-convergent limit. Evjen's method remains the standard
textbook technique for evaluating Madelung constants by direct summation.

*Implementation:* ``chemistrykit.crystal.utils.lattice_sums.evjen_lattice_sum_cubic_alternating()``
implements exactly this fractional-boundary-weighting scheme over a simple
cubic lattice of alternating unit charges -- the object NaCl's Madelung
sum reduces to once its two interpenetrating FCC sublattices are combined
onto one finer simple-cubic mesh -- and its module docstring includes a
runnable demonstration that a naive, unweighted truncated sum genuinely
does oscillate (standard deviation over `n_shells` = 6..10 well above
0.05) while the Evjen-weighted sum has already converged
(standard deviation below :math:`10^{-4}`) over the same range;
:func:`~chemistrykit.crystal.madelung_constant_nacl`
is the NaCl-specific wrapper that consumes it.

*References:* H. M. Evjen, "On the Stability of Certain Heteropolar
Crystals," Phys. Rev. 39 (1932), 675-687.

.. minigallery:: ../../examples/crystal/madelung/plot_01_madelung_convergence.py

1956 -- Kapustinskii's Structure-Independent Lattice-Energy Equation
---------------------------------------------------------------------------

Anatoli Kapustinskii noticed an empirical regularity in Born-Lande-type
lattice energies across many known salts: the Madelung constant of a given
structure type, divided by the number of ions per formula unit, is nearly
the same from one structure to the next. Exploiting that near-constancy,
he proposed a lattice-energy equation built entirely from ionic radii and
charges -- no Madelung constant, and hence no assumed crystal structure,
required at all:

.. math::

   U = -\frac{\kappa\,\nu\,|z_+z_-|}{r_++r_-}\left(1-\frac{d}{r_++r_-}\right)

with :math:`\kappa` and `d` empirical constants fitted so the equation
reproduces Born-Lande and experimental lattice energies to within about
5%. Because it needs no prior knowledge of a salt's actual crystal
structure -- only its ionic radii, which can be estimated for essentially
any ion pair -- the Kapustinskii equation became, and remains, the
standard quick estimate for the lattice energy of a hypothetical or
poorly characterized salt, exactly the situation the structure-dependent
Born-Lande equation above cannot handle at all.

*Implementation:* :class:`chemistrykit.crystal.systems.lattice_energy.KapustinskiiLatticeEnergy`
implements exactly this equation, sharing the
:class:`~chemistrykit.crystal.core.base_system.LatticeEnergyModel`
interface with :class:`~chemistrykit.crystal.systems.lattice_energy.BornLandeLatticeEnergy`
so the two independently-derived estimates can be compared directly on
the same salt.

*References:* A. F. Kapustinskii, "Lattice Energy of Ionic Crystals," Q.
Rev. Chem. Soc. 10 (1956), 283-294.

.. minigallery:: ../../examples/crystal/lattice_energy/plot_02_kapustinskii.py

1976 -- Shannon's Revised Effective Ionic Radii
-----------------------------------------------------

R. D. Shannon compiled and critically re-analyzed thousands of measured
interatomic distances in oxide and halide crystals to produce a
self-consistent table of "effective ionic radii" for essentially every
common ion, superseding the older, less systematic scales of Goldschmidt
and Pauling from the 1920s. Shannon's key methodological insight was to
treat the radius not as a single fixed number per ion but as a
value depending explicitly on coordination number and, for transition
metals, on spin state -- so that an ion's tabulated radius comes with the
coordination environment it was fitted for, rather than being applied
uniformly regardless of context. The resulting table, still in routine
use exactly as published, is the standard reference source for ionic
radii throughout structural and solid-state chemistry, radius-ratio
rules, and, as below, the Kapustinskii equation's radius-only lattice-
energy estimates.

*Implementation:* ``chemistrykit.crystal.utils.reference_data.SHANNON_IONIC_RADII_PM``
tabulates Shannon's 6-coordinate (octahedral) effective ionic radii for
the common main-group ions used in this domain's Kapustinskii-equation
examples and tests -- consumed directly as the `r_cation_pm`/`r_anion_pm`
arguments of :class:`chemistrykit.crystal.systems.lattice_energy.KapustinskiiLatticeEnergy`.

*References:* R. D. Shannon, "Revised Effective Ionic Radii and
Systematic Studies of Interatomic Distances in Halides and
Chalcogenides," Acta Cryst. A32 (1976), 751-767.

.. minigallery:: ../../examples/crystal/lattice_energy/plot_03_shannon_ionic_radii.py

See Also
--------

- :doc:`/api/crystal`
- :doc:`/history/analytical_breakthroughs`
