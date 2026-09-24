r"""
Kasha's rule: emission comes from the lowest excited state
============================================================

Kasha (1950) stated that emission comes almost entirely from the lowest
excited state of a given multiplicity (:math:`S_1` for fluorescence),
because internal conversion from higher states such as :math:`S_2` takes
picoseconds or less, much faster than emission. With
:func:`~chemistrykit.photochem.kasha_emission_yields`, a molecule excited
to :math:`S_2` sends only a tiny fraction of its emission from
:math:`S_2`; the exception is a molecule like azulene, whose large
:math:`S_2`-:math:`S_1` gap slows internal conversion enough for
:math:`S_2` emission to compete.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem import kasha_emission_yields

kf2, kf1, knr1 = 1.0e8, 1.0e8, 1.0e8  # 1/s
k_ic21 = np.logspace(7, 14, 200)  # S2 -> S1 internal conversion rate, 1/s
phi2, phi1 = np.array([kasha_emission_yields(kf2, k, kf1, knr1) for k in k_ic21]).T

fig, ax = plt.subplots()
ax.semilogx(k_ic21, phi1 / (phi1 + phi2), label=r"share of emission from $S_1$")
ax.semilogx(k_ic21, phi2 / (phi1 + phi2), label=r"share of emission from $S_2$")
ax.axvspan(1e12, 1e14, color="0.9", label="typical $S_2\\to S_1$ internal conversion")
ax.set_xlabel(r"$k_{ic}(S_2\to S_1)$ (s$^{-1}$)")
ax.set_ylabel("Fraction of emitted photons")
ax.set_title("Kasha's rule: fast internal conversion funnels emission to $S_1$")
ax.legend(loc="center left")
fig.tight_layout()

# %%
for label, k in (("typical dye", 1e13), ("azulene-like (slow S2 -> S1)", 1e9)):
    p2, p1 = kasha_emission_yields(kf2, k, kf1, knr1)
    print(f"{label:30s}: Phi(S2 emission) = {p2:.2e}, Phi(S1 emission) = {p1:.3f}")

plt.show()
