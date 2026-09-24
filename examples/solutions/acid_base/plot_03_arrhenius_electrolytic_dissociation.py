r"""
Arrhenius's electrolytic dissociation: how much of an electrolyte is ions?
============================================================================

Arrhenius proposed that an acid dissolved in water is split, at least in
part, into free ions even with no current flowing, and that the
*degree of dissociation* :math:`\alpha` distinguishes a strong
electrolyte (:math:`\alpha \approx 1`) from a weak one
(:math:`\alpha \ll 1`). The number of dissolved particles per formula
unit is then :math:`i = 1 + \alpha` for a 1:1 electrolyte -- the van 't
Hoff factor that explained anomalous colligative data.

Here :math:`\alpha` is computed with
:meth:`~chemistrykit.solutions.core.base_system.WeakElectrolyte.percent_dissociation`
for 0.10 M solutions of acids spanning twelve orders of magnitude in
:math:`K_a`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.acid_base import WeakAcid

Ca = 0.10
Ka_values = np.logspace(-10, 2, 200)
alpha = np.array([WeakAcid(Ca=Ca, Ka=Ka).percent_dissociation() / 100.0 for Ka in Ka_values])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
ax1.semilogx(Ka_values, alpha, color="steelblue")
for name, Ka in [("HCN", 6.2e-10), ("acetic", 1.8e-5), ("HF", 6.8e-4), ("HSO$_4^-$", 1.2e-2)]:
    a = WeakAcid(Ca=Ca, Ka=Ka).percent_dissociation() / 100.0
    ax1.plot(Ka, a, "o", color="darkorange")
    ax1.annotate(name, (Ka, a), textcoords="offset points", xytext=(5, -12), fontsize=8)
ax1.set_xlabel("acid dissociation constant $K_a$")
ax1.set_ylabel(r"degree of dissociation $\alpha$")
ax1.set_title("0.10 M acid: from weak to strong electrolyte")

ax2.semilogx(Ka_values, 1.0 + alpha, color="seagreen")
ax2.set_xlabel("acid dissociation constant $K_a$")
ax2.set_ylabel(r"van 't Hoff factor $i = 1 + \alpha$")
ax2.set_title("Dissolved particles per formula unit")
fig.tight_layout()

# %%
# Printed values for a few familiar acids at 0.10 M:

for name, Ka in [("hydrocyanic acid", 6.2e-10), ("acetic acid", 1.8e-5), ("hydrofluoric acid", 6.8e-4)]:
    acid = WeakAcid(Ca=Ca, Ka=Ka)
    print(f"{name:18s}: {acid.percent_dissociation():6.3f}% dissociated, i = {1 + acid.percent_dissociation() / 100:.4f}")

plt.show()
