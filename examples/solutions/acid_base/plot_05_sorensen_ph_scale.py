r"""
Sorensen's pH scale: a logarithm for hydrogen-ion concentration
=================================================================

Sorensen's 1909 definition :math:`\mathrm{pH} = -\log_{10}[H^+]`
(:func:`~chemistrykit.solutions.systems.acid_base.ph_from_h`) turns
hydrogen-ion concentrations spanning fourteen orders of magnitude into
numbers between 0 and 14. Below, the pH of several everyday solutions is
computed from their equilibrium :math:`[H^+]` and placed on the scale,
and the companion relation :math:`\mathrm{pH} + \mathrm{pOH} = pK_w = 14`
is checked with :func:`~chemistrykit.solutions.systems.acid_base.poh_from_oh`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.acid_base import WeakAcid, WeakBase, h_from_ph, ph_from_h, poh_from_oh

Kw = 1.0e-14
samples = {
    "0.1 M HCl": 0.1,
    "vinegar (0.8 M acetic)": WeakAcid(Ca=0.8, Ka=1.8e-5).h_concentration(),
    "0.01 M acetic acid": WeakAcid(Ca=0.01, Ka=1.8e-5).h_concentration(),
    "pure water": np.sqrt(Kw),
    "0.1 M ammonia": Kw / WeakBase(Cb=0.1, Kb=1.8e-5).oh_concentration(),
    "0.1 M NaOH": Kw / 0.1,
}

pH_axis = np.linspace(0, 14, 200)
fig, ax = plt.subplots(figsize=(8, 5))
ax.semilogy(pH_axis, h_from_ph(pH_axis), color="steelblue", label=r"$[H^+] = 10^{-\mathrm{pH}}$")
for name, h in samples.items():
    pH = float(ph_from_h(h))
    ax.plot(pH, h, "o", color="darkorange")
    ax.annotate(name, (pH, h), textcoords="offset points", xytext=(6, 4), fontsize=8)
ax.set_xlabel("pH")
ax.set_ylabel("[H+] (mol/L)")
ax.set_title("Sorensen's pH scale")
ax.legend()
fig.tight_layout()

# %%
# pH and pOH of each sample always sum to :math:`pK_w = 14` at 25 degC:

for name, h in samples.items():
    pH = float(ph_from_h(h))
    pOH = float(poh_from_oh(Kw / h))
    print(f"{name:24s}: [H+] = {h:.2e} M, pH = {pH:5.2f}, pOH = {pOH:5.2f}, sum = {pH + pOH:.2f}")

plt.show()
