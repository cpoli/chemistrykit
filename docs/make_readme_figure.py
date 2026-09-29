"""Regenerate the README hero figure: python docs/make_readme_figure.py

Writes docs/source/_static/images/readme_hero.png, which README.md embeds by
its raw.githubusercontent.com URL so it also renders on PyPI.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import chemistrykit as ck
from chemistrykit.spectro.visualizers.spectro_plots import plot_broadened_spectrum

OUT = Path(__file__).parent / "source" / "_static" / "images" / "readme_hero.png"

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)

# The Belousov-Zhabotinsky reaction's relaxation oscillations (Field-Noyes Oregonator).
result = ck.kinetics.Oregonator(f=1.0).integrate((0.0, 40.0), method="dopri5", rtol=1e-7, atol=1e-10, max_steps=2_000_000)
for series, label in zip(result.y.T, ("[HBrO$_2$]", "[Br$^-$]", "[Ce(IV)]"), strict=True):
    ax1.semilogy(result.t, series, label=label)
ax1.set_ylim(top=1e4)
ax1.set_xlabel("t (scaled)")
ax1.set_ylabel("concentration (scaled)")
ax1.set_title("An oscillating reaction: Belousov-Zhabotinsky")
ax1.legend(fontsize=8, ncol=3, loc="upper center")

# Ethanol's 1H NMR: the CH3 triplet and CH2 quartet of the n+1 rule.
ppm = np.linspace(0.8, 4.1, 6000)
for shift, neighbors, label in ((1.2, 2, "CH$_3$ triplet"), (3.7, 3, "CH$_2$ quartet")):
    multiplet = ck.spectro.first_order_multiplet(chemical_shift_ppm=shift, j_coupling_hz=7.0, n_neighbors=neighbors, spectrometer_frequency_mhz=60.0)
    plot_broadened_spectrum(multiplet, ppm, ax=ax2, shape="lorentzian", fwhm=0.005, label=label)
ax2.invert_xaxis()
ax2.set_xlabel("chemical shift (ppm)")
ax2.set_title("Ethanol's $^1$H NMR spectrum: spin-spin splitting")
ax2.legend(fontsize=8)

# van der Waals isotherms of CO2 around its critical point.
Tc, Pc = 304.13, 7.3773e6
vdw = ck.thermo.VanDerWaals.from_critical_constants(Tc, Pc)
Vc = 3.0 * vdw.b
Vm = np.linspace(1.4 * vdw.b, 8.0 * Vc, 600)
for Tr in (0.85, 0.9, 0.95, 1.0, 1.1, 1.2):
    ax3.plot(Vm / Vc, vdw.pressure(Vm, Tr * Tc) / Pc, "k-" if Tr == 1.0 else "-", label=rf"$T/T_c$ = {Tr}")
ax3.plot([1.0], [1.0], "ro")
ax3.set_xlim(0.4, 8.0)
ax3.set_ylim(0.0, 2.0)
ax3.set_xlabel(r"$V_m / V_c$")
ax3.set_ylabel(r"$P / P_c$")
ax3.set_title("van der Waals: the liquid-gas critical point")
ax3.legend(fontsize=7)

fig.savefig(OUT, dpi=110)
print(f"wrote {OUT}")
