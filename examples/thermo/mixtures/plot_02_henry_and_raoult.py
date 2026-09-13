r"""
Henry's law and Raoult's law as two limits of a real solution
==================================================================

A real (non-ideal) solute's partial vapor pressure is *not* a straight
line all the way from :math:`x_B=0` to :math:`x_B=1`. But it is tangent
to two different straight lines at its two ends: Henry's law
(:func:`~chemistrykit.thermo.systems.mixtures.henry_law_pressure`),
with slope :math:`K_H`, as the solute becomes infinitely dilute, and
Raoult's law (:func:`~chemistrykit.thermo.systems.mixtures.raoult_vapor_pressure`),
with slope :math:`P_B^*` (the pure-liquid vapor pressure), as the
"solute" becomes the only thing present. A simple quadratic
interpolation between the two limits, built only to match both slopes
exactly, illustrates the textbook picture without needing any real
non-ideal-mixture model.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.mixtures import henry_law_pressure, raoult_vapor_pressure

# A hypothetical solute B whose measured Henry's law constant (the
# initial slope at infinite dilution) differs substantially from its
# pure vapor pressure (the slope at x_B = 1) -- exactly the generic
# situation for a real, non-ideal solute/solvent pair.
P_B_star = 40.0  # kPa, pure-B vapor pressure -- the Raoult's-law slope
K_H = 120.0  # kPa, Henry's law constant -- the dilute-limit slope

x_B = np.linspace(0.0, 1.0, 200)
# The unique quadratic in x_B that equals 0 at x_B=0, has slope K_H at
# x_B=0, and equals P_B_star at x_B=1: an illustrative real-solution
# curve tangent to both limiting laws.
P_B_real = x_B * (K_H + (P_B_star - K_H) * x_B)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(x_B, P_B_real, color="black", linewidth=2, label="real (illustrative) partial pressure")
ax.plot(x_B, raoult_vapor_pressure(x_B, P_B_star), "--", color="steelblue", label=r"Raoult's law ($x_B P_B^*$)")
ax.plot(x_B, henry_law_pressure(x_B, K_H), "--", color="darkorange", label=r"Henry's law ($x_B K_H$)")
ax.set_xlabel("mole fraction of B")
ax.set_ylabel("partial pressure of B (kPa)")
ax.set_title("A real solute: Henry's law when dilute, Raoult's law when nearly pure")
ax.legend()
fig.tight_layout()

# %%
# Near x_B = 0 (B highly dilute in solvent A) the real curve hugs the
# Henry's-law line; near x_B = 1 (B nearly pure) it hugs the Raoult's-law
# line instead. Both are exact limiting laws of the same real solution,
# valid in opposite composition extremes -- neither is "wrong," they are
# just tangent lines to the same curve at its two ends.

for x_test, label in [(0.02, "dilute (x_B=0.02)"), (0.98, "nearly pure (x_B=0.98)")]:
    real = float(x_test * (K_H + (P_B_star - K_H) * x_test))
    henry = float(henry_law_pressure(x_test, K_H))
    raoult = float(raoult_vapor_pressure(x_test, P_B_star))
    print(f"{label}: real={real:.4f}, Henry's-law prediction={henry:.4f}, Raoult's-law prediction={raoult:.4f}")

plt.show()
