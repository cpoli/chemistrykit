Breakthroughs in Electrochemistry
====================================


.. include:: /_generated/nav/electrochem.rst

.. epigraph::

   "The apparatus of which I speak, and which will doubtless astonish
   you, is only the assemblage of a number of good conductors of
   different kinds arranged in a certain manner." -- Alessandro Volta,
   letter to Sir Joseph Banks announcing the voltaic pile, 20 March 1800

Electrochemistry began as an argument between two Italians over whether
electricity lived in living tissue or in the contact of two metals, and
resolved itself, within a generation, into a quantitative science of
charge, potential, and rate -- one exact enough that a battery's
discharge curve, an electroplating bath's yield, and a corroding pipe's
current can all be predicted from the same handful of equations. The
systems in :mod:`chemistrykit.electrochem` retrace that arc: from the
first cell and the first electrolysis, through Faraday's laws and
nomenclature, the Nernst and Butler-Volmer equations that connect
electrode potential and current to concentration and kinetics, to the
international convention that finally fixed what an electrode potential
even means. This chronology traces the major breakthroughs behind the
package, with a pointer to the corresponding implementation at each stop.

.. contents:: Timeline
   :local:
   :depth: 1

1780s -- 1800 -- Galvani, Volta, and the Animal-Electricity Dispute
--------------------------------------------------------------------

In the mid-1780s, Luigi Galvani, a professor of anatomy at Bologna,
found that a dissected frog's leg twitched when its nerve and muscle
were touched simultaneously by two different metals joined in a
circuit -- and, in another version of the experiment, when a scalpel
touched an exposed nerve near a working electrostatic machine. Galvani
concluded he had discovered "animal electricity," an intrinsic
electrical fluid secreted by the nerves and stored in the muscle,
with the metals merely serving as passive conductors that discharged it.
Alessandro Volta, a physicist at Pavia, at first accepted Galvani's
interpretation but, on repeating and extending the experiments, grew
convinced the frog's leg was not the source of the electricity at all but
merely a sensitive detector of it: what actually mattered, Volta argued,
was the contact of the two *different* metals themselves, with any moist
tissue (or, as he would soon show, brine-soaked cardboard) serving only
as an electrolytic go-between. The two men's dispute over this point ran
for the rest of Galvani's life. Volta settled the argument the way
physics arguments are best settled -- by building a device that worked
without any animal tissue whatsoever: stacking pairs of zinc and silver
(or copper) discs separated by brine-soaked cloth into a "voltaic pile,"
he produced a steady electric current from inert materials alone, and
announced it to the Royal Society in a letter to Sir Joseph Banks on 20
March 1800. Both men, in the end, were partly right: Volta correctly
identified metal-metal contact (more precisely, as later electrochemistry
would clarify, the difference in the two metals' electrode potentials) as
the driving force, while Galvani's frog nerve-muscle preparation really
was, as he insisted, an exquisitely sensitive electrical detector --
just not the source.

*Connection:* Volta's pile -- a stack of identical two-metal cells wired
in series to boost the available voltage -- is the direct ancestor of
every battery in this package, and the potential difference it exploits
between two different metals is exactly what
``chemistrykit.electrochem.systems.standard_potentials.STANDARD_REDUCTION_POTENTIALS``
tabulates for each metal's reduction couple, and what
:func:`~chemistrykit.electrochem.cell_potential`
combines into a cell voltage (see 1836, below, for the two-metal cell
this package's own worked example uses).
:class:`chemistrykit.electrochem.systems.battery.ConstantCurrentBattery`
is this module's modern, quantitative descendant of the pile: a
load-driven source characterized by its terminal voltage and remaining
capacity rather than, as for Volta, an open scientific question of what
was flowing and why.

*References:* A. Volta, letter to Sir Joseph Banks, 20 March 1800,
published as "On the Electricity Excited by the Mere Contact of
Conducting Substances of Different Kinds," Phil. Trans. R. Soc. Lond. 90,
403-431 (1800); L. Galvani, *De Viribus Electricitatis in Motu Musculari
Commentarius* (Bologna, 1791).

.. minigallery:: ../../examples/electrochem/battery/plot_02_voltaic_pile.py

1800 -- Nicholson and Carlisle: the First Electrolysis
--------------------------------------------------------

Word of Volta's pile reached London within weeks. William Nicholson and
Anthony Carlisle built one of their own, and -- almost immediately,
apparently while just testing an electrical connection through a drop of
water meant to improve contact between wires -- noticed gas bubbling from
both submerged wire ends. Investigating deliberately, they found hydrogen
collecting at one wire and oxygen at the other, in a volume ratio of
roughly two to one: the pile had decomposed water into its elements, the
first chemical transformation ever driven by an electric current, and the
first genuinely new scientific result the pile made possible. They
published the observation within days, in Nicholson's own *Journal of
Natural Philosophy, Chemistry, and the Arts* in May 1800 -- electrolysis
as a phenomenon preceded, by more than three decades, Faraday's
quantitative laws of exactly how much substance a given charge transforms
(see 1834, below).

*Implementation:* :func:`~chemistrykit.electrochem.minimum_applied_voltage_electrolytic`,
applied to the H+/H2 and O2/H2O half-reactions already in
``STANDARD_REDUCTION_POTENTIALS``
via :func:`~chemistrykit.electrochem.cell_potential`,
reproduces the textbook 1.23 V theoretical decomposition voltage of
water, and :func:`~chemistrykit.electrochem.moles_from_charge`
recovers the 2:1 hydrogen-to-oxygen mole (and, by Avogadro's law, volume)
ratio Nicholson and Carlisle observed, directly from the two half-reactions'
electron counts (2 for H2, 4 for O2).

*References:* W. Nicholson, "Account of the New Electrical or Galvanic
Apparatus of Sig. Alex. Volta, and Experiments Performed with the Same,"
J. Nat. Philos. Chem. Arts 4, 179-187 (1800), reporting the electrolysis
observed with A. Carlisle.

.. minigallery:: ../../examples/electrochem/electrolysis/plot_02_water_electrolysis.py

1807 -- Davy and the Electrolytic Isolation of the Alkali Metals
--------------------------------------------------------------------

Humphry Davy, at the Royal Institution, reasoned that if electrolysis
could split water, a strong enough current applied to a molten (rather
than dissolved) alkali might split apart compounds no chemical reagent
had ever managed to decompose. Passing current through molten potash
(potassium hydroxide) in 1807, he watched small globules of a new,
silvery, intensely reactive metal appear at the cathode and immediately
catch fire in the moist air -- potassium, the first alkali metal ever
isolated, followed within days by sodium from molten soda. Davy announced
both in his Bakerian Lecture to the Royal Society that November. The
achievement was only possible *because* these metals are so
overwhelmingly difficult to reduce by ordinary chemical means: their
standard reduction potentials are among the most negative of any element,
meaning an enormous amount of energy must be supplied to force the
reverse (reducing) reaction -- energy no chemical reducing agent of the
day could provide, but a sufficiently large voltaic battery could.

*Connection:* the potassium and sodium half-reactions Davy's electrolysis
targeted are two of the most negative entries in
``chemistrykit.electrochem.systems.standard_potentials.STANDARD_REDUCTION_POTENTIALS``
(:math:`K^++e^-\to K`, :math:`-2.93` V; :math:`Na^++e^-\to Na`,
:math:`-2.71` V) -- exactly the quantitative statement of why these
metals resisted every purely chemical reduction attempt before Davy's,
and why :func:`~chemistrykit.electrochem.minimum_applied_voltage_electrolytic`
returns such a large forcing voltage when applied to either couple.

*References:* H. Davy, "The Bakerian Lecture, on Some New Phenomena of
Chemical Changes Produced by Electricity," Phil. Trans. R. Soc. Lond. 98,
1-44 (1808) (the lecture was read 19 November 1807; the page range is as
commonly cited in secondary literature and has not been independently
verified against the original volume).

.. minigallery:: ../../examples/electrochem/electrolysis/plot_03_davy_alkali_metals.py

1833 -- 1834 -- Faraday's Laws of Electrolysis and Electrochemical Nomenclature
------------------------------------------------------------------------------------

Michael Faraday spent much of the early 1830s putting electrolysis on a
quantitative footing for the first time, in a long series of papers to
the Royal Society titled, collectively, *Experimental Researches in
Electricity*. His two laws of electrolysis established, first, that the
mass of a substance liberated at an electrode is directly proportional to
the quantity of electric charge passed through the cell, and second, that
for a fixed charge, the mass liberated of different substances is
proportional to their equivalent weight -- together implying a single
universal constant of proportionality (later named the Faraday constant
in his honor) connecting charge to chemical change, independent of the
electrolyte, the electrodes, or the current used to pass it. Faraday was
also, by his own admission, unhappy with the vague and inconsistent
vocabulary electrochemistry had used up to that point, and turned for
help to William Whewell, the Cambridge polymath and philosopher of
science already known for coining "scientist" itself. Between them they
settled on the vocabulary still used unchanged nearly two centuries
later: "electrode" for the two conducting surfaces, "anode" and "cathode"
for the electrodes of entry and exit for the (then still hypothetical)
current, "electrolyte" for the decomposed substance, and "ion" (from the
Greek for "to go") for the electrolyte's migrating charged constituents.

.. math::

   m = \frac{QM}{nF}, \qquad Q = It

*Implementation:* :func:`~chemistrykit.electrochem.charge_from_current`
computes :math:`Q=It`; :func:`~chemistrykit.electrochem.moles_from_charge`
and :func:`~chemistrykit.electrochem.mass_from_charge`
implement Faraday's first law directly, and
:func:`~chemistrykit.electrochem.faradays_law_mass`
combines all three into the single-call form used throughout this
gallery -- its doctest verifies the defining linearity (doubling the
current exactly doubles the deposited mass) that was Faraday's original
experimental result.

*References:* M. Faraday, "Experimental Researches in Electricity,
Seventh Series," Phil. Trans. R. Soc. Lond. 124, 77-122 (1834) (the laws
of electrolysis); the electrochemical vocabulary is introduced in the
same series, with Faraday crediting William Whewell's assistance by name.

.. minigallery:: ../../examples/electrochem/electrolysis/plot_01_faradays_laws.py

1836 -- Daniell's Two-Fluid Cell
------------------------------------

Volta's simple pile had a serious practical flaw: as current was drawn
from it, hydrogen gas generated at the copper electrode formed an
insulating film that "polarized" the cell, causing its voltage to sag and
its current to fade within minutes. John Frederic Daniell, a professor of
chemistry at King's College London, solved the problem by physically
separating the two half-reactions into different compartments (or,
in his own construction, different fluids kept apart by a porous
earthenware barrier): a zinc electrode in dilute sulfuric acid on one
side, and a copper electrode in copper sulfate solution on the other. On
the copper side, rather than hydrogen gas, the reaction deposits solid
copper metal directly onto the electrode -- there is no gas to polarize
anything, and the cell delivers a genuinely steady voltage for as long as
its reactants last. Daniell announced the design in a letter to Faraday
read before the Royal Society in February 1836, and it became the first
electrochemical cell reliable and reproducible enough to serve as a
practical voltage reference and a working power source for early
telegraphy.

*Implementation:* :func:`~chemistrykit.electrochem.cell_potential`
and :func:`~chemistrykit.electrochem.standard_cell_potential`
combine exactly the Cu2+/Cu and Zn2+/Zn half-reactions Daniell paired,
both already in
``STANDARD_REDUCTION_POTENTIALS``
-- their doctests reproduce the textbook 1.10 V Daniell-cell potential
directly, and
:func:`~chemistrykit.electrochem.balance_redox_reaction`
gives the (here trivial, both sides already 2-electron) mass-balancing
multiples for the overall reaction.

*References:* J. F. Daniell, "On Voltaic Combinations," Phil. Trans. R.
Soc. Lond. 126, 107-124 (1836).

.. minigallery:: ../../examples/electrochem/standard_potentials/plot_01_daniell_cell.py

1839 -- 1842 -- Grove's Gas Battery: the First Fuel Cell
-------------------------------------------------------------

William Robert Grove, a Welsh lawyer turned experimental physicist,
reasoned that if a current could split water into hydrogen and oxygen,
the reverse reaction ought to be able to produce a current. In 1839 he
immersed two platinum strips in dilute sulfuric acid, their upper ends
sealed in tubes of hydrogen and of oxygen, and found a steady current
flowing between them; by 1842 he had built a "gaseous voltaic battery"
of such cells in series, strong enough to electrolyze water itself --
the gas battery driving the very decomposition it reversed. Grove also
recognized the essential role of the platinum surface, where gas,
electrolyte and electrode meet. His device is the direct ancestor of
the modern hydrogen fuel cell, whose electrical output is governed not
by the heat of combustion :math:`\Delta H` but by the Gibbs energy
:math:`\Delta G`, so that its efficiency is not bounded by the Carnot
factor of a heat engine:

.. math::

   E = -\frac{\Delta G}{nF}, \qquad \eta_{max} = \frac{\Delta G}{\Delta H}

*Implementation:* :func:`~chemistrykit.electrochem.reversible_cell_voltage`
converts :math:`\Delta G^\circ = -237.13` kJ/mol into the 1.23 V
reversible voltage of the hydrogen-oxygen cell (matching
:func:`~chemistrykit.electrochem.standard_cell_potential` for the O2/H2O
and H+/H2 couples),
:func:`~chemistrykit.electrochem.reversible_cell_voltage_at_temperature`
gives its temperature dependence, and
:func:`~chemistrykit.electrochem.fuel_cell_efficiency_limit` the 83%
thermodynamic efficiency limit.

*References:* W. R. Grove, "On Voltaic Series and the Combination of
Gases by Platinum," Phil. Mag. 14, 127-130 (1839); W. R. Grove, "On a
Gaseous Voltaic Battery," Phil. Mag. 21, 417-420 (1842).

.. minigallery:: ../../examples/electrochem/fuel_cell/plot_01_grove_gas_battery.py

1876 -- 1900 -- Kohlrausch's Laws of Electrolytic Conductivity
-------------------------------------------------------------------

Friedrich Kohlrausch turned the conductivity of salt solutions into a
precision measurement, using alternating current to avoid electrode
polarization. Comparing salts with a common ion, he found in 1876 that
the difference between, say, the potassium and sodium salts was the same
whatever the anion: at infinite dilution each ion migrates
independently and contributes its own fixed share,
:math:`\lambda_i^\circ`, to the molar conductivity -- the law of
independent migration of ions. His measurements on very dilute solutions
over the following decades, collected in his 1898 monograph with Ludwig
Holborn, established the second law: for strong electrolytes the molar
conductivity falls linearly with the square root of concentration,
allowing extrapolation to infinite dilution. Combining the two lets the
limiting conductivity of a *weak* electrolyte such as acetic acid, whose
curve cannot be extrapolated, be assembled from those of strong ones --
a value Arrhenius and Ostwald soon needed for their degree of
dissociation.

.. math::

   \Lambda_m^\circ = \sum_i \nu_i\lambda_i^\circ, \qquad
   \Lambda_m = \Lambda_m^\circ - K\sqrt{c}

*Implementation:* :func:`~chemistrykit.electrochem.limiting_molar_conductivity`
sums ionic contributions from
``chemistrykit.electrochem.systems.conductivity.LIMITING_IONIC_CONDUCTIVITIES``
(its doctest checks the common-ion difference that Kohlrausch observed),
:func:`~chemistrykit.electrochem.kohlrausch_molar_conductivity` implements
the square-root law, and :func:`~chemistrykit.electrochem.fit_kohlrausch_law`
performs the extrapolation to infinite dilution, returning a
:class:`~chemistrykit.electrochem.KohlrauschFit`.

*References:* F. Kohlrausch, Nachr. Königl. Ges. Wiss. Göttingen (1876)
(independent migration of ions); F. Kohlrausch and L. Holborn, *Das
Leitvermögen der Elektrolyte* (Teubner, Leipzig, 1898).

.. minigallery:: ../../examples/electrochem/conductivity/plot_01_kohlrausch_laws.py

1889 -- Nernst's Equation
------------------------------

Walther Nernst, applying the young science of chemical thermodynamics
(much of it his own) to the electrochemical cell, derived a single
equation relating a cell's actual potential under arbitrary conditions to
its standard potential and the concentrations (strictly, activities) of
the species involved -- explaining, for the first time on a quantitative
thermodynamic footing, why a cell's voltage sags as its reactants are
consumed and its products accumulate, and why two half-cells built from
identical chemistry but different concentrations (a "concentration cell")
can generate a voltage from concentration difference alone, with no
standard-potential difference at all. The equation earned Nernst the 1920
Nobel Prize in Chemistry, awarded specifically "in recognition of his
work in thermochemistry."

.. math::

   E = E^\circ - \frac{RT}{nF}\ln Q

*Implementation:* :func:`~chemistrykit.electrochem.nernst_potential`
implements exactly this equation; its doctest confirms it reduces to
:math:`E^\circ` at :math:`Q=1` and that increasing product-side activity
lowers the cell potential, per Le Chatelier's principle.
:func:`~chemistrykit.electrochem.concentration_cell_potential`
specializes it to :math:`E^\circ=0`, reproducing the textbook ~59 mV per
decade of concentration ratio for a one-electron couple at 25 degC.

*References:* W. Nernst, "Die elektromotorische Wirksamkeit der Ionen,"
Z. Phys. Chem. 4, 129-181 (1889).

.. minigallery:: ../../examples/electrochem/nernst/plot_01_nernst_equation.py

1897 -- Peukert's Law of Battery Capacity
------------------------------------------------

Studying the lead-acid cells used in the earliest electric vehicles and
telephone-exchange backup supplies, the German scientist Wilhelm Peukert
found an inconvenient empirical fact: a battery's delivered capacity is
not, in practice, independent of how fast it is discharged. Discharged
quickly, a real cell delivers noticeably *less* total charge before its
voltage collapses than the same cell discharged slowly -- an effect
absent from an idealized battery whose rated ampere-hour capacity is a
fixed number, but very much present in practice because a high discharge
rate lets the depleted-reactant layer at the electrode surface outrun the
rate at which fresh electrolyte can diffuse in to replace it. Peukert
captured the effect in a single empirical power law relating discharge
time to current, with an exponent (the "Peukert exponent," typically
1.1-1.3 for lead-acid chemistries) that quantifies exactly how much worse
than ideal a given battery's high-rate performance is -- a purely
empirical curve fit, with no claim to a first-principles diffusion
derivation, that nonetheless remains standard practice in battery
engineering more than a century later.

.. math::

   t = \frac{C_p}{I^k}, \qquad C_{eff}(I) = I\,t(I) = C_p I^{1-k}

*Implementation:* :func:`~chemistrykit.electrochem.peukert_discharge_time`
implements exactly this power law, and
:func:`~chemistrykit.electrochem.effective_capacity`
the delivered-capacity-vs-rate relationship it implies -- constant at
:math:`k=1` (the idealized, rate-independent battery), falling with
increasing current at the realistic :math:`k>1` lead-acid case.
:class:`~chemistrykit.electrochem.systems.battery.ConstantCurrentBattery`
combines this Peukert-law runtime with a constant ohmic voltage sag into
this package's simplified discharge-curve model (see its docstring, and
:class:`chemistrykit.electrochem.core.base_system.BatteryDischargeModel`,
for exactly what is and is not captured).

*References:* W. Peukert, "Über die Abhängigkeit der Kapazität von der
Entladestromstärke bei Bleiakkumulatoren," Elektrotechnische Zeitschrift
18 (1897); see also Linden & Reddy, *Handbook of Batteries*, 3rd ed.,
Ch. 3.3, for the modern textbook treatment this module follows.

.. minigallery:: ../../examples/electrochem/battery/plot_01_discharge_curves.py

1903 -- Cottrell and the Diffusion-Limited Current
------------------------------------------------------

Frederick Gardner Cottrell, working in Leipzig on the "residual current"
that flows after an electrode's potential is changed, treated it as a
problem in Fick's diffusion. If the potential is stepped so that the
electroactive species is consumed as soon as it reaches a planar
electrode, a depleted layer grows into the solution like
:math:`\sqrt{Dt}`, and the current -- the rate at which fresh material
diffuses in -- decays as :math:`t^{-1/2}`. The Cottrell equation was the
first quantitative link between an electrode current and mass transport,
and underlies chronoamperometry and the measurement of diffusion
coefficients to this day.

.. math::

   i(t) = nFAC\sqrt{\frac{D}{\pi t}}

*Implementation:* :func:`~chemistrykit.electrochem.cottrell_current`
implements the equation; its doctest confirms the constant
:math:`i\sqrt{t}` signature, and the unit tests check that integrating
it reproduces the accumulated charge :math:`2nFAC\sqrt{Dt/\pi}`.

*References:* F. G. Cottrell, "Der Reststrom bei galvanischer
Polarisation, betrachtet als ein Diffusionsproblem," Z. Phys. Chem. 42,
385-431 (1903).

.. minigallery:: ../../examples/electrochem/voltammetry/plot_01_cottrell_equation.py

1905 -- Tafel's Empirical Overpotential Law
--------------------------------------------------

Julius Tafel, studying the kinetics of hydrogen evolution at a mercury
cathode, found that once an electrode reaction is driven far enough from
equilibrium (a large enough "overpotential" applied beyond the reaction's
equilibrium potential), the logarithm of the resulting current density
grows *linearly* with the overpotential -- a simple empirical
straight-line relationship, found decades before anyone had a
microscopic kinetic theory of *why* it should hold. That explanation
would not arrive until the Butler-Volmer equation two decades later (see
1924/1930, below), whose high-overpotential limit reduces exactly to
Tafel's linear law -- but Tafel's own result, purely empirical and purely
kinetic, remains the standard practical tool for measuring an electrode
reaction's exchange current density from experimental current-voltage
data, independent of whichever microscopic model is eventually invoked
to explain it.

.. math::

   \eta = b\log_{10}\!\left(\frac{i}{i_0}\right), \qquad b_{anodic} = \frac{2.303RT}{\alpha nF}

*Implementation:* :func:`~chemistrykit.electrochem.tafel_slope`
computes exactly this slope from the charge-transfer coefficient, and
:func:`~chemistrykit.electrochem.tafel_overpotential`
the resulting linearized overpotential-current relationship;
:func:`~chemistrykit.electrochem.fit_tafel_plot`
performs the reverse operation Tafel's own experimentalists needed --
recovering the Tafel slope and exchange current density by linear
regression against measured (:math:`\eta`, :math:`\log_{10}i`) data.

*References:* J. Tafel, "Über die Polarisation bei kathodischer
Wasserstoffentwicklung," Z. Phys. Chem. 50, 641-712 (1905).

.. minigallery:: ../../examples/electrochem/butler_volmer/plot_01_tafel.py

1922 -- Heyrovský's Polarography
-------------------------------------

Jaroslav Heyrovský, in Prague, measured the current through a dropping
mercury electrode -- a capillary releasing a fresh, perfectly clean
mercury drop every few seconds -- as the applied voltage was slowly
increased. Each reducible species in the solution produced a sigmoidal
step, a "polarographic wave": the potential at half its height (the
half-wave potential) identifies the species, and the height of its
plateau is proportional to its concentration. With Masuzo Shikata he
automated the recording in the polarograph (1925), making polarography
the first widely used instrumental method of electroanalysis. Dionýz
Ilkovič derived the height of the plateau from diffusion to the growing
drop (1934), and Heyrovský and Ilkovič the shape of the reversible wave
(1935). Heyrovský received the 1959 Nobel Prize in Chemistry.

.. math::

   \bar{i}_d = 607\,n D^{1/2} m^{2/3} t^{1/6} C, \qquad
   E = E_{1/2} + \frac{RT}{nF}\ln\frac{i_d - i}{i}

*Implementation:* :func:`~chemistrykit.electrochem.ilkovic_diffusion_current`
implements the Ilkovič equation (maximum or drop-averaged, in its
traditional microampere units), and
:func:`~chemistrykit.electrochem.polarographic_wave_current` the
Heyrovský-Ilkovič wave, whose doctest confirms the current is exactly half
the plateau at the half-wave potential.

*References:* J. Heyrovský, "Elektrolysa se rtuťovou kapkovou kathodou,"
Chem. Listy 16, 256-264 (1922); J. Heyrovský and M. Shikata, Recl. Trav.
Chim. Pays-Bas 44, 496-498 (1925); D. Ilkovič, Collect. Czech. Chem.
Commun. 6, 498-513 (1934); J. Heyrovský and D. Ilkovič, Collect. Czech.
Chem. Commun. 7, 198-214 (1935).

.. minigallery:: ../../examples/electrochem/voltammetry/plot_02_heyrovsky_polarography.py

1923 -- Debye-Hückel Theory and Nonideal Electrolytes
------------------------------------------------------------

Peter Debye and Erich Hückel showed that an ion in solution is
surrounded, on average, by a diffuse cloud of counter-ions -- statistical
rather than a fixed structure, but strong enough to systematically screen
the ion's own electric field and lower its effective (thermodynamic)
activity below its raw concentration. Their theory gives a closed-form
prediction for exactly how much an ion's activity coefficient departs
from the ideal value of 1, as a function of the solution's total ionic
strength and the ion's own charge -- explaining, for the first time from
first principles, why real electrolyte solutions systematically deviate
from ideal (Nernstian, unit-activity-coefficient) behavior more strongly
as concentration and ionic charge increase.

*Implementation:* ``chemistrykit.electrochem.systems.nernst`` reuses
this package's existing Debye-Hückel machinery
(``chemistrykit.solutions.systems.activity``, developed for solution
equilibria) rather than reimplementing it: :func:`~chemistrykit.electrochem.activity_corrected_reaction_quotient`
builds a reaction quotient from Debye-Hückel-corrected activities instead
of raw concentrations, and
:func:`~chemistrykit.electrochem.nernst_potential_with_activity`
feeds that corrected quotient into the Nernst equation directly -- its
doctest confirms the correction vanishes in the dilute limit, recovering
the ideal Nernst potential, exactly as Debye-Hückel theory predicts.

*References:* P. Debye and E. Hückel, "Zur Theorie der Elektrolyte," Phys.
Z. 24, 185-206 (1923).

.. minigallery:: ../../examples/electrochem/nernst/plot_02_debye_huckel_activity.py

1924, 1930 -- Butler, Erdey-Grúz and Volmer: the Butler-Volmer Equation
------------------------------------------------------------------------

John Alfred Valentine Butler, in a pair of 1924 papers, worked out how
the Nernst equation's equilibrium electrode potential could arise
*kinetically*, from a genuinely dynamic balance of forward and reverse
electron-transfer rates at the electrode surface rather than from
equilibrium thermodynamics alone. Tibor Erdey-Grúz and Max Volmer
completed the picture in 1930, deriving the full current-overpotential
relationship in the exponential form still used today: two competing
exponential terms, one for the anodic (oxidation) and one for the
cathodic (reduction) partial current, that exactly cancel at zero
overpotential -- a genuinely dynamic equilibrium, with both partial
reactions still running at full (and equal) speed, rather than a
thermodynamic "nothing is happening." At high overpotential in either
direction, one exponential term dominates completely and the equation
reduces exactly to Tafel's earlier empirical straight line (see 1905,
above), finally giving Tafel's law a microscopic kinetic derivation.

.. math::

   i = i_0\left[\exp\!\left(\frac{\alpha nF\eta}{RT}\right)
                -\exp\!\left(-\frac{(1-\alpha)nF\eta}{RT}\right)\right]

*Implementation:* :func:`~chemistrykit.electrochem.butler_volmer_current_density`
implements exactly this equation; its doctest confirms the net current
vanishes exactly at zero overpotential, the equation's defining dynamic-
equilibrium property.
:func:`~chemistrykit.electrochem.exchange_current_density`
computes :math:`i_0` itself from a standard heterogeneous rate constant
and bulk concentrations, and the module's own numerical check (comparing
:func:`~chemistrykit.electrochem.tafel_overpotential`
against the full equation at 300 mV overpotential, agreeing to better
than 0.5%) is a direct demonstration of the high-overpotential
Tafel-law limit described above.

*References:* J. A. V. Butler, "Studies in Heterogeneous Equilibria,"
Trans. Faraday Soc. 19, 729-733 (1924); T. Erdey-Grúz and M. Volmer,
"Zur Theorie der Wasserstoff Überspannung," Z. Phys. Chem. A 150, 203-213
(1930).

.. minigallery:: ../../examples/electrochem/butler_volmer/plot_02_butler_volmer_equation.py

1948 -- Randles, Ševčík and the Theory of Linear-Sweep Voltammetry
------------------------------------------------------------------------

Instead of stepping the potential, John Randles in Birmingham and
Augustin Ševčík in Prague swept it linearly in time, recording the
current on a cathode-ray oscillograph. Both, independently and in the
same year, solved the resulting diffusion problem for a reversible
couple: the current rises to a peak and then decays as the depletion
layer grows, and the peak current is proportional to concentration and
to the square root of the scan rate. Their result is the foundation of
linear-sweep and cyclic voltammetry, today the most widely used
electroanalytical technique; a linear plot of :math:`i_p` against
:math:`\sqrt{v}` is the standard test that an electrode reaction is
diffusion-controlled.

.. math::

   i_p = 0.4463\,nFAC\sqrt{\frac{nFvD}{RT}}

*Implementation:* :func:`~chemistrykit.electrochem.randles_sevcik_peak_current`
implements the equation in SI units; its doctest checks the
:math:`\sqrt{v}` scaling, and the unit tests recover the familiar
:math:`2.69\times10^5` coefficient of the 25 degC textbook form.

*References:* J. E. B. Randles, "A cathode ray polarograph. Part II. The
current-voltage curves," Trans. Faraday Soc. 44, 327-338 (1948);
A. Ševčík, "Oscillographic polarography with periodical triangular
voltage," Collect. Czech. Chem. Commun. 13, 349-377 (1948).

.. minigallery:: ../../examples/electrochem/voltammetry/plot_03_randles_sevcik_cyclic_voltammetry.py

1953 -- The Stockholm Convention and the Modern Sign Convention
--------------------------------------------------------------------

For most of the nineteenth and early twentieth centuries, electrochemists
worked with two genuinely incompatible sign conventions for "the"
electrode potential of a half-reaction -- one (associated with Nernst and
later Latimer) reporting *oxidation* potentials, positive for a
spontaneous oxidation, and the other reporting *reduction* potentials,
positive for a spontaneous reduction -- with the same half-reaction
carrying opposite-signed tabulated values depending on which textbook
tradition a chemist had learned from. The International Union of Pure and
Applied Chemistry finally settled the ambiguity at its 1953 Stockholm
meeting, recommending unanimously that "the electrode potential" of a
half-reaction always mean its *reduction* potential -- the convention
this package, and every modern textbook table, now uses without comment.

*Implementation:* every value in
``chemistrykit.electrochem.systems.standard_potentials.STANDARD_REDUCTION_POTENTIALS``
is tabulated as a Stockholm-convention reduction potential, and
:func:`~chemistrykit.electrochem.cell_potential`'s
cathode-minus-anode combination rule, together with
:func:`~chemistrykit.electrochem.is_spontaneous`'s
sign convention (:math:`E_{cell}>0\iff` spontaneous), work correctly
precisely *because* both half-reactions are consistently expressed in
this single, internationally agreed sign convention rather than a
per-reaction mix of the two historical alternatives.

*References:* International Union of Pure and Applied Chemistry, "Report
of the Committee on Electrochemical Nomenclature," adopted at the 17th
IUPAC Conference, Stockholm, 1953; summarized in the electrochemistry
sections of the IUPAC *Green Book* (*Quantities, Units and Symbols in
Physical Chemistry*) in all subsequent editions.

.. minigallery:: ../../examples/electrochem/standard_potentials/plot_02_stockholm_sign_convention.py

See Also
--------

- :doc:`/api/electrochem`
- :doc:`/history/photochem_breakthroughs`
