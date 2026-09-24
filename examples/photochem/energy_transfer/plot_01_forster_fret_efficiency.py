r"""
Förster resonance energy transfer: the inverse-sixth-power distance law
=========================================================================

Förster (1948) showed that an excited donor can hand its energy to a
nearby acceptor through dipole-dipole coupling at a rate
:math:`k_T=\tau_D^{-1}(R_0/r)^6`, so the transfer efficiency
:math:`E=1/(1+(r/R_0)^6)` falls from nearly 1 to nearly 0 over a narrow
range around the Förster distance :math:`R_0`
(:func:`~chemistrykit.photochem.forster_radius`,
:func:`~chemistrykit.photochem.forster_efficiency`). That steepness is
why FRET works as a "spectroscopic ruler" for distances of 1-10 nm.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem import forster_efficiency, forster_radius, forster_rate

R0 = forster_radius(kappa2=2 / 3, n=1.4, quantum_yield_donor=0.5, overlap_J=2.0e15) / 10.0  # nm
tau_D = 4.0e-9
print(f"Forster distance R0 = {R0:.2f} nm")

r = np.linspace(1.0, 12.0, 400)
E = forster_efficiency(r, R0)
fig, ax = plt.subplots()
ax.plot(r, E)
ax.axvline(R0, color="k", ls=":", label=rf"$R_0$ = {R0:.1f} nm, $E$ = 0.5")
ax.axvspan(0.5 * R0, 1.5 * R0, color="0.92", label=r"useful ruler range $0.5$-$1.5\,R_0$")
ax.set_xlabel("Donor-acceptor distance $r$ (nm)")
ax.set_ylabel("Transfer efficiency $E$")
ax.set_title(r"FRET efficiency $1/(1+(r/R_0)^6)$")
ax.legend()
fig.tight_layout()

# %%
# The donor lifetime shortens with transfer: tau_DA = 1/(1/tau_D + k_T),
# and 1 - tau_DA/tau_D reproduces E exactly.
for d in (0.5 * R0, R0, 1.5 * R0):
    tau_DA = 1.0 / (1.0 / tau_D + forster_rate(d, R0, tau_D))
    print(f"r = {d:4.2f} nm: E = {forster_efficiency(d, R0):.3f}, 1 - tau_DA/tau_D = {1 - tau_DA / tau_D:.3f}")

plt.show()
