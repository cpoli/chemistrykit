r"""
Stokes shift: fluorescence is emitted at longer wavelength than it is absorbed
================================================================================

Stokes (1852) found that quinine solution absorbs invisible ultraviolet
light and re-emits it as blue light -- always at longer wavelength than
the exciting light. The excited molecule relaxes vibrationally (and its
solvent reorganizes) before emitting, so emission starts from a lower
energy than absorption ended at.
:func:`~chemistrykit.photochem.stokes_shift` gives the gap in
wavenumbers. The bands below are model Gaussian spectra with maxima near
those of quinine sulfate in dilute acid (about 350 nm and 450 nm).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem import stokes_shift

lam_abs, lam_em = 350.0, 450.0  # nm
wavenumber = np.linspace(16000.0, 34000.0, 800)  # cm^-1
nu_abs, nu_em = 1e7 / lam_abs, 1e7 / lam_em
absorption = np.exp(-(((wavenumber - nu_abs) / 2000.0) ** 2))
emission = np.exp(-(((wavenumber - nu_em) / 2000.0) ** 2))

shift = stokes_shift(lam_abs, lam_em)
print(f"Stokes shift: {shift:.0f} cm^-1 ({lam_em - lam_abs:.0f} nm)")

fig, ax = plt.subplots()
ax.plot(1e7 / wavenumber, absorption, label="absorption")
ax.plot(1e7 / wavenumber, emission, label="fluorescence")
ax.annotate(
    "",
    xy=(lam_em, 1.03),
    xytext=(lam_abs, 1.03),
    arrowprops={"arrowstyle": "->"},
)
ax.text(0.5 * (lam_abs + lam_em), 1.06, f"Stokes shift = {shift:.0f} cm$^{{-1}}$", ha="center")
ax.set_ylim(0, 1.15)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Normalized intensity")
ax.set_title("Stokes shift of a quinine-like fluorophore")
ax.legend(loc="upper right")
fig.tight_layout()

plt.show()
