r"""
Volumetric titration: locating the equivalence point
======================================================

The volumetric analysis Gay-Lussac and Mohr systematized finds how much
acid a sample contains from the burette volume at which it has been
exactly neutralized -- the equivalence point, seen by eye as an
indicator's sudden colour change and here located numerically as the
steepest point of the pH curve
(:meth:`~chemistrykit.solutions.core.base_system.Titration.find_equivalence_point`).
Comparing a strong-acid/strong-base titration
(:class:`~chemistrykit.solutions.systems.titration.StrongAcidStrongBaseTitration`)
against a weak-acid/strong-base titration
(:class:`~chemistrykit.solutions.systems.titration.WeakAcidStrongBaseTitration`)
of the same concentrations: the weak-acid curve starts at a higher pH,
has a buffer region around its half-equivalence point (where pH = pKa),
and its equivalence point is basic rather than neutral -- all visible
directly from the curve shapes.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.titration import (
    StrongAcidStrongBaseTitration,
    WeakAcidStrongBaseTitration,
)

Ca, Va, Cb = 0.100, 0.050, 0.100
Ka = 1.8e-5  # acetic acid

strong = StrongAcidStrongBaseTitration(Ca=Ca, Va=Va, Cb=Cb)
weak = WeakAcidStrongBaseTitration(Ca=Ca, Va=Va, Ka=Ka, Cb=Cb)

Vb = np.linspace(1e-6, 0.09, 2000)
result_strong = strong.curve(Vb)
result_weak = weak.curve(Vb)

Vb_eq_strong = strong.find_equivalence_point(Vb)
Vb_eq_weak = weak.find_equivalence_point(Vb)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(result_strong.Vb * 1000.0, result_strong.pH, label="strong acid + strong base", color="steelblue")
ax.plot(result_weak.Vb * 1000.0, result_weak.pH, label="weak acid + strong base", color="darkorange")
ax.axvline(Vb_eq_strong * 1000.0, color="steelblue", linestyle=":", alpha=0.6)
ax.axvline(Vb_eq_weak * 1000.0, color="darkorange", linestyle=":", alpha=0.6)
ax.axhline(7.0, color="gray", linestyle="--", linewidth=0.7, label="pH 7")
ax.set_xlabel("volume of NaOH added (mL)")
ax.set_ylabel("pH")
ax.set_title("Titration of 50.0 mL of 0.100 M acid with 0.100 M NaOH")
ax.legend()
fig.tight_layout()

# %%
# The numerically detected equivalence points (steepest-ascent points)
# for each curve, compared against the exact stoichiometric volume
# (identical for both, since it only depends on moles of acid and base,
# not on Ka):

print(f"Stoichiometric equivalence volume: {strong.equivalence_volume() * 1000.0:.3f} mL")
print(f"Strong-acid titration: numeric equivalence at {Vb_eq_strong * 1000.0:.3f} mL, pH there = 7.00 (neutral)")
print(f"Weak-acid titration: numeric equivalence at {Vb_eq_weak * 1000.0:.3f} mL, pH there > 7 (basic)")

plt.show()
