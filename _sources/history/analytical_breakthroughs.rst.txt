Breakthroughs in Analytical Chemistry
========================================


.. include:: /_generated/nav/analytical.rst

.. epigraph::

   "It is a capital mistake to theorize before one has data. Insensibly
   one begins to twist facts to suit theories, instead of theories to
   suit facts." -- Sherlock Holmes, in Arthur Conan Doyle, "A Scandal in
   Bohemia," *The Strand Magazine*, 1891 -- quoted here for its literal
   truth in this domain: every technique below exists to turn a noisy
   measurement into a number with a defensible uncertainty attached.

Analytical chemistry is the discipline of answering "how much, and how
sure are we" -- separating a mixture into measurable components, titrating
a solution to a detectable endpoint, fitting a calibration line, and
deciding, by a defensible statistical rule rather than a hunch, whether a
stray data point should be thrown out. The systems in
:mod:`chemistrykit.analytical` retrace the handful of ideas that turned
each of those tasks from craft into quantitative method: chromatographic
separation and its plate theory, potentiometric and complexometric
titration, linear regression and detection limits, propagation of
uncertainty, and outlier rejection. This chronology traces the major
conceptual breakthroughs behind the package, with a pointer to the
corresponding implementation at each stop.

.. contents:: Timeline
   :local:
   :depth: 1

1805 -- Legendre, Gauss, and the Method of Least Squares
-----------------------------------------------------------

Adrien-Marie Legendre published the first description of fitting a
straight line (or any linear model) to noisy data by minimizing the sum of
squared residuals, as an appendix -- "Sur la Méthode des moindres
quarrés" -- to a memoir on determining the orbits of comets. Carl Friedrich
Gauss published his own, more fully probabilistic derivation of the same
method four years later, and claimed (very plausibly, on the strength of
his own earlier working notes and later astronomical use) to have been
using it in private practice since 1795, touching off a priority dispute
that has never been fully resolved either way in the historical record.
Whichever came first, the method of least squares gave quantitative
science its basic tool for extracting a best-fit relationship, and an
honest estimate of that relationship's uncertainty, from measurements that
scatter around a straight line rather than falling exactly on one -- the
single technique underneath essentially every calibration curve run in an
analytical laboratory since.

*Implementation:* ``chemistrykit.analytical.utils.regression.linear_fit()``
performs exactly this ordinary-least-squares fit
(``y = slope*x + intercept``), and additionally reports the residual
standard error :math:`s_{y/x}` that a calibration curve's detection limits
(below, 1983) are built from;
:func:`~chemistrykit.analytical.fit_calibration` wraps
it into a :class:`~chemistrykit.analytical.systems.calibration.LinearCalibration`
usable directly for predicting concentration from a measured signal.

*References:* A.-M. Legendre, *Nouvelles méthodes pour la détermination
des orbites des comètes* (Paris: Firmin Didot, 1805), appendix "Sur la
Méthode des moindres quarrés"; C. F. Gauss, *Theoria Motus Corporum
Coelestium in Sectionibus Conicis Solem Ambientium* (Hamburg: Perthes &
Besser, 1809), Book II, Sec. III.

.. minigallery:: ../../examples/analytical/calibration/plot_01_calibration_curve.py

1889 -- Nernst's Equation and Potentiometric Redox Chemistry
-----------------------------------------------------------------

Walther Nernst derived the general thermodynamic relationship between an
electrode's measured potential and the activities (in the dilute-solution
limit, concentrations) of the species in the redox couple that sets it:

.. math::

   E = E^\circ - \frac{RT}{nF}\ln Q

with `Q` the reaction quotient of the half-reaction and `n` the number of
electrons transferred. The equation gave electrochemistry its quantitative
bridge between a directly measurable quantity -- a voltage -- and the
otherwise inaccessible ratio of oxidized to reduced species in solution,
and it very quickly became the basis of a new titration methodology:
rather than watching a color-change indicator, an analyst could follow a
redox titration's progress by measuring the electrode potential directly
and locating the endpoint at the potential's point of steepest change,
sidestepping indicators unsuitable for some redox couples entirely.

*Implementation:* :class:`chemistrykit.analytical.systems.titration.RedoxTitration`
applies the Nernst equation to both half-reactions of a redox titration in
the standard large-equilibrium-constant approximation -- the analyte
couple's potential before the equivalence volume, the titrant couple's
potential after it, and the classical weighted-average result
:math:`E_{eq}=(n_1E^\circ_1+n_2E^\circ_2)/(n_1+n_2)` exactly at
equivalence -- with :meth:`~chemistrykit.analytical.TitrationCurve.find_equivalence_point`
(inherited from the shared
:class:`~chemistrykit.analytical.core.base_system.TitrationCurve` base)
locating the endpoint numerically, as the point of steepest potential
change, exactly the way a potentiometric titration is read in practice.

*References:* W. Nernst, "Die elektromotorische Wirksamkeit der Jonen,"
Z. Phys. Chem. 4 (1889), 129-181.

.. minigallery:: ../../examples/analytical/titration/plot_01_titration_curves.py

1901 -- 1906 -- Tsvet and the Invention of Chromatography
--------------------------------------------------------------

Mikhail Tsvet, studying plant pigments, packed a glass column with
powdered calcium carbonate, poured a petroleum-ether extract of leaf
pigments through it, and watched the pigments separate into a series of
distinctly colored bands as they traveled down the column at different
rates -- the different pigments adsorbing to, and desorbing from, the solid
packing with different affinities. Tsvet coined the name
"chromatography" (literally "color-writing") for the technique in his
1906 papers, and used it to demonstrate that what had been thought to be a
single pigment, chlorophyll, was in fact a mixture of several distinct
compounds. Tsvet's method was largely ignored for over three decades --
partly eclipsed by skepticism from established chemists of the era, and
by his own early death in 1919 -- before being independently rediscovered
and extended into the dominant separation technique of modern analytical
chemistry, beginning with Martin and Synge's work below.

*Implementation:* Tsvet's column is the physical apparatus that every
formula in ``chemistrykit.analytical.systems.chromatography`` is built
to quantify -- a mixture separating into discrete, differently-retained
bands as it migrates through a stationary phase --
:func:`~chemistrykit.analytical.simulate_chromatogram`
renders that outcome directly, as a sum of separately-retained Gaussian
elution peaks on a simulated detector trace, the modern instrumental
descendant of Tsvet's visually banded column.

*References:* M. Tswett, "Adsorptionsanalyse und chromatographische
Methode. Anwendung auf die Chemie des Chlorophylls," Ber. Dtsch. Bot.
Ges. 24 (1906), 384-393, and "Physikalisch-chemische Studien über das
Chlorophyll. Die Adsorptionen," same volume, 316-323.

.. minigallery:: ../../examples/analytical/chromatography/plot_01_van_deemter.py

1941 -- Martin, Synge, and the Theoretical-Plate Model
-------------------------------------------------------------

Archer Martin and Richard Synge introduced partition chromatography --
separating compounds by their relative solubility between two liquid
phases (one held stationary on an inert support) rather than by
adsorption onto a solid, as Tsvet's original method had -- and, in the same
paper, borrowed a piece of theory from an entirely different field to
describe how a chromatographic band broadens as it migrates: the
theoretical-plate model already used to describe the stepwise
equilibration of a fractional-distillation column. Treating a
chromatographic column as a stack of many small, discrete equilibration
stages ("theoretical plates") gave the field its first quantitative
measure of column efficiency, extractable directly from an eluted peak's
retention time and width:

.. math::

   N = 16\left(\frac{t_R}{w_{base}}\right)^2

Martin and Synge received the 1952 Nobel Prize in Chemistry for the
development of partition chromatography.

*Implementation:* :func:`~chemistrykit.analytical.theoretical_plates`
implements exactly this formula (and its equivalent full-width-at-half-
maximum form, :math:`N=5.545(t_R/w_{1/2})^2`, verified in its own
docstring to agree with the base-width form for a Gaussian peak of
consistent shape); :func:`~chemistrykit.analytical.plate_height`
converts a plate count into the column-length-per-plate `H` that the van
Deemter equation below predicts directly.

*References:* A. J. P. Martin and R. L. M. Synge, "A New Form of
Chromatogram Employing Two Liquid Phases," Biochem. J. 35 (1941),
1358-1368.

.. minigallery:: ../../examples/analytical/chromatography/plot_01_van_deemter.py

1945 -- Schwarzenbach and EDTA Complexometric Titration
--------------------------------------------------------------

Gerold Schwarzenbach, together with E. Kampitsch and R. Steiner, opened a
long series of papers (titled, collectively, "Komplexone") introducing
aminopolycarboxylic acid ligands -- ethylenediaminetetraacetic acid (EDTA)
foremost among them -- as general-purpose chelating titrants for metal
ions. EDTA's six potential donor atoms wrap around a single metal cation
to form an unusually stable 1:1 complex regardless of the metal's specific
charge or coordination preference, letting one titrant, buffered to an
appropriate pH, quantify essentially any of dozens of different metal ions
via a color-change or potentiometric endpoint -- founding complexometric
titration as a general analytical method rather than a family of
metal-specific reactions.

*Implementation:* :class:`chemistrykit.analytical.systems.titration.EDTATitration`
models exactly this 1:1 complexation :math:`M+Y\rightleftharpoons MY` via
its *conditional* (pH-corrected) formation constant :math:`K_f'`, solving
the resulting mass-action quadratic in closed form for the free-metal
concentration at each titrant volume, and verifying in its own docstring
that the exact solution approaches the textbook large-:math:`K_f'`
approximation :math:`pM\approx\frac12\log_{10}(K_f'/C_{M,eq})` at the
equivalence point as :math:`K_f'` grows.

*References:* G. Schwarzenbach, E. Kampitsch, and R. Steiner, "Komplexone
I. Über die Salzbildung der Nitrilotriessigsäure," Helv. Chim. Acta 28
(1945), 828-840, the first of the "Komplexone" series that introduced
EDTA and its relatives as general complexometric titrants over the
following years.

.. minigallery:: ../../examples/analytical/titration/plot_01_titration_curves.py

1950 -- 1991 -- Dixon's Q-test and Rorabacher's Revised Critical Values
------------------------------------------------------------------------------

W. J. Dixon developed a family of simple statistical tests for rejecting a
single suspect outlier from a small data set -- too small for the standard
large-sample outlier tests of the day to apply usefully -- based on the
ratio of the "gap" between the suspect value and its nearest neighbor to
the full range of the data:

.. math::

   Q = \frac{\text{gap}}{\text{range}}

A simplified version of the test, restricted to the single most common
case (rejecting the smallest or largest of `n` replicate measurements) and
aimed squarely at working analytical chemists rather than statisticians,
was popularized two years later by Robert Dean and Dixon himself, and
became -- and remains -- the standard quick test for a suspect replicate
measurement in analytical practice. The critical values in that original
1951 table, however, were computed under a less precise approximation than
later became available; D. B. Rorabacher recomputed and extended the
critical-value table in 1991 using more accurate Monte Carlo methods, and
it is Rorabacher's revised table, not the original Dean-Dixon values, that
this package (and most modern textbooks) actually use.

*Implementation:* :func:`~chemistrykit.analytical.dixon_q_test`
implements exactly this gap-over-range statistic and compares it against
``Q_CRITICAL_TABLE``, which
reproduces Rorabacher's revised (not the original Dean-Dixon) critical
values for sample sizes 3-10 at 90%, 95%, and 99% confidence.

*References:* W. J. Dixon, "Analysis of Extreme Values," Ann. Math.
Statist. 21 (1950), 488-506; R. B. Dean and W. J. Dixon, "Simplified
Statistics for Small Numbers of Observations," Anal. Chem. 23 (1951),
636-638; D. B. Rorabacher, "Statistical Treatment for Rejection of
Deviant Values: Critical Values of Dixon's 'Q' Parameter and Related
Subrange Ratios at the 95% Confidence Level," Anal. Chem. 63 (1991),
139-146.

.. minigallery:: ../../examples/analytical/qtest/plot_01_dixon_q_test.py

1956 -- van Deemter, Zuiderweg, and Klinkenberg: The van Deemter Equation
------------------------------------------------------------------------------

J. J. van Deemter, F. J. Zuiderweg, and A. Klinkenberg, working at Royal
Dutch Shell's laboratories, identified and separated the three physically
distinct mechanisms that broaden a chromatographic band as it migrates
through a packed column, and combined them into a single equation for
plate height `H` as a function of the mobile phase's linear velocity `u`:

.. math::

   H = A + \frac{B}{u} + Cu

`A` (eddy diffusion) reflects the multiple unequal flow paths a packed bed
offers; `B` (longitudinal molecular diffusion) dominates at low velocity,
where a band has time to spread by ordinary diffusion; and `C`
(resistance to mass transfer between the mobile and stationary phases)
dominates at high velocity, where equilibration between the two phases can
no longer keep pace with the flow. Because `B` and `C` push plate height in
opposite directions as `u` changes, the equation predicts -- and
experiments confirm -- a single optimum flow velocity,
:math:`u_{opt}=\sqrt{B/C}`, at which column efficiency is maximized; the
van Deemter equation remains the standard framework for understanding and
optimizing chromatographic column performance to this day.

*Implementation:* :func:`~chemistrykit.analytical.van_deemter_H`
implements exactly this equation;
:func:`~chemistrykit.analytical.optimum_flow_velocity`
and :func:`~chemistrykit.analytical.minimum_plate_height`
give the closed-form optimum velocity and minimum plate height from
:math:`dH/du=0`, each verified in its own docstring against a direct
numerical scan of :func:`~chemistrykit.analytical.van_deemter_H` itself.

*References:* J. J. van Deemter, F. J. Zuiderweg, and A. Klinkenberg,
"Longitudinal Diffusion and Resistance to Mass Transfer as Causes of
Nonideality in Chromatography," Chem. Eng. Sci. 5 (1956), 271-289.

.. minigallery:: ../../examples/analytical/chromatography/plot_01_van_deemter.py

1958 -- Golay and the Open-Tubular Column
-----------------------------------------------

Marcel Golay extended van Deemter's plate-height theory to a column
geometry packed columns cannot offer: a long, narrow, open capillary tube
with the stationary phase coated directly on its inner wall rather than
held on a packed solid support. Because an open tube has no packing
particles at all, it has no unequal packed-bed flow paths to speak of, and
Golay showed that the eddy-diffusion term simply vanishes -- the van
Deemter equation's `A` term drops to zero, leaving a column whose
efficiency, for a given length, can substantially exceed a packed
column's. The resulting "Golay equation" made possible the long, narrow,
extremely efficient capillary columns that displaced packed columns as
the standard for gas chromatography within a generation.

*Implementation:* an open-tubular column is not modeled as a separate
class in this package, but is exactly the :math:`A=0` special case of the
same :func:`~chemistrykit.analytical.van_deemter_H`
already used for the general (packed-column) case above -- the identical
equation, with one physically motivated term switched off, rather than a
separate model requiring its own implementation.

*References:* M. J. E. Golay, "Theory of Chromatography in Open and
Coated Tubular Columns with Round and Rectangular Cross-Sections," in
*Gas Chromatography 1958* (Amsterdam Symposium), ed. D. H. Desty
(London: Butterworths, 1958), 36-55.

.. minigallery:: ../../examples/analytical/chromatography/plot_01_van_deemter.py

1966 -- Ku's NBS Formalization of Uncertainty Propagation
-----------------------------------------------------------------

Harry Ku, at the National Bureau of Standards (now NIST), wrote what
became the standard expository reference for the propagation-of-error
formula already in wide informal use across the physical sciences by the
1960s -- consolidating scattered practice into a single systematic
statement of the first-order (linearized) propagation law,

.. math::

   \sigma_y=\sqrt{\sum_i\left(\frac{\partial f}{\partial x_i}\right)^2\sigma_{x_i}^2}

for a function :math:`y=f(x_1,\ldots,x_n)` of independent, uncorrelated
measured quantities, together with worked derivations of the specific
sum, product, and power shortcuts that follow from it and the conditions
under which the underlying linearization is (and is not) a good
approximation. Ku's paper, still cited in modern metrology guidance, is
the standard citation for the propagation-of-uncertainty formula as
practiced in analytical and physical chemistry laboratories today.

*Implementation:* :func:`~chemistrykit.analytical.propagate_uncertainty`
evaluates the general formula above directly via numerical partial
derivatives, for a function with no simple closed form;
:func:`~chemistrykit.analytical.propagate_sum`,
:func:`~chemistrykit.analytical.propagate_product`, and
:func:`~chemistrykit.analytical.propagate_power` give
the closed-form sum/product/power shortcuts of the same formula, each
verified in this package's tests to agree with the general numerical
form for the corresponding operation.

*References:* H. H. Ku, "Notes on the Use of Propagation of Error
Formulas," J. Res. Natl. Bur. Stand. Sect. C 70C (1966), 263-273.

.. minigallery:: ../../examples/analytical/uncertainty/plot_01_error_propagation.py

1983 -- Long, Winefordner, and the IUPAC LOD/LOQ Convention
--------------------------------------------------------------------

Gary Long and James Winefordner surveyed the many mutually inconsistent
definitions of an instrument's "limit of detection" then in circulation
across the analytical literature, and argued for standardizing on a
single statistically grounded convention, tied directly to a calibration
curve's own residual scatter about its fit rather than to an arbitrary,
method-specific rule of thumb: a signal is detectable once it exceeds the
blank by 3.3 standard deviations of the calibration residuals (limiting
the false-positive and false-negative rates to a commonly agreed
tolerance under a Gaussian noise model), and reliably quantifiable once it
exceeds the blank by 10 such standard deviations.

.. math::

   \text{LOD} = \frac{3.3\,s_{y/x}}{|m|}, \qquad \text{LOQ} = \frac{10\,s_{y/x}}{|m|}

with `m` the calibration curve's slope and :math:`s_{y/x}` its residual
standard error. Their recommendation was adopted by IUPAC and remains the
standard convention for reporting an analytical method's detection and
quantitation limits today.

*Implementation:* :meth:`~chemistrykit.analytical.LinearCalibration.lod`
and :meth:`~chemistrykit.analytical.LinearCalibration.loq`
implement exactly these two formulas, drawing the slope `m` and residual
standard error :math:`s_{y/x}` directly from the
:func:`~chemistrykit.analytical.fit_calibration`
least-squares fit (see 1805, above) -- so that LOQ is, by construction,
always exactly :math:`10/3.3` times LOD, a ratio fixed purely by Long and
Winefordner's convention and independent of any particular data set.

*References:* G. L. Long and J. D. Winefordner, "Limit of Detection: A
Closer Look at the IUPAC Definition," Anal. Chem. 55 (1983), 712A-724A.

.. minigallery:: ../../examples/analytical/calibration/plot_01_calibration_curve.py

See Also
--------

- :doc:`/api/analytical`
- :doc:`/history/crystal_breakthroughs`
