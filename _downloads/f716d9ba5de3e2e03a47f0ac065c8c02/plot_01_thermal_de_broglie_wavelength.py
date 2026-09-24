r"""
De Broglie's matter waves: the thermal wavelength of a gas
=============================================================

De Broglie's relation :math:`\lambda=h/p`, averaged over a gas's thermal
spread of momenta, gives the thermal de Broglie wavelength

.. math::

    \Lambda=\frac{h}{\sqrt{2\pi mk_BT}},

computed by :func:`~chemistrykit.statmech.thermal_de_broglie_wavelength`.
Translational motion can be treated classically when :math:`\Lambda` is
much smaller than the mean spacing between molecules,
:math:`(V/N)^{1/3}=(k_BT/P)^{1/3}`. The plot compares the two for several
gases at 1 bar: only the lightest particles at the lowest temperatures
come close.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

from chemistrykit.constants import K_B
from chemistrykit.statmech import thermal_de_broglie_wavelength

species = {
    "electron": (sc.m_e, "black"),
    "He": (4.0026 * sc.atomic_mass, "crimson"),
    "H2": (2.016 * sc.atomic_mass, "darkorange"),
    "N2": (28.014 * sc.atomic_mass, "seagreen"),
    "Ar": (39.948 * sc.atomic_mass, "steelblue"),
}

T = np.logspace(0.0, 3.5, 300)
P = 1.0e5
spacing = (K_B * T / P) ** (1.0 / 3.0)

fig, ax = plt.subplots(figsize=(7, 5))
for name, (mass, color) in species.items():
    ax.loglog(T, thermal_de_broglie_wavelength(mass, T), color=color, label=rf"$\Lambda$ ({name})")
ax.loglog(T, spacing, color="gray", linestyle="--", linewidth=2, label=r"mean spacing $(k_BT/P)^{1/3}$ at 1 bar")
ax.set_xlabel("T (K)")
ax.set_ylabel("length (m)")
ax.set_title("Thermal de Broglie wavelength vs. intermolecular spacing")
ax.legend(fontsize=8)
fig.tight_layout()

# %%
# At room temperature Lambda is a small fraction of an angstrom for any
# molecule, some 100-1000 times smaller than the ~34 angstrom spacing of a
# gas at 1 bar -- the quantitative reason classical translational
# statistics work. An electron's wavelength, by contrast, is comparable to
# that spacing, and much larger than the spacing in a metal.

T_room = 298.15
d = (K_B * T_room / P) ** (1.0 / 3.0)
for name, (mass, _) in species.items():
    lam = thermal_de_broglie_wavelength(mass, T_room)
    print(f"{name:8s}: Lambda = {lam * 1e10:8.3f} angstrom, spacing / Lambda = {d / lam:8.1f}")

plt.show()
