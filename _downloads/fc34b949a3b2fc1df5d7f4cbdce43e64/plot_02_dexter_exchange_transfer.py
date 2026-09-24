r"""
Dexter exchange energy transfer: exponential fall-off with distance
=====================================================================

Dexter (1953) described a second energy-transfer mechanism in which the
donor and acceptor swap electrons, which needs their orbitals to overlap.
Its rate, :math:`k_{ET}=KJ\exp(-2r/L)`
(:func:`~chemistrykit.photochem.dexter_rate`), therefore dies off
exponentially within a few Å of contact, and because electron exchange
conserves total spin it can move triplet energy (triplet sensitization),
which dipole-dipole Förster coupling cannot. The log-scale plot shows the
straight line of an exponential law, with the slope set by :math:`L`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem import dexter_rate

r = np.linspace(0.0, 15.0, 300)  # edge-to-edge distance, Angstrom
K, J = 1.0e13, 1.0  # s^-1, normalized overlap
fig, ax = plt.subplots()
for L in (1.0, 2.0, 3.0):
    ax.semilogy(r, dexter_rate(r, K, J, L), label=f"$L$ = {L:.0f} Å")
ax.axhline(1.0e5, color="k", ls=":", label=r"triplet decay rate ~$10^5$ s$^{-1}$")
ax.set_ylim(1e-2, 2e13)
ax.set_xlabel("Donor-acceptor edge-to-edge distance $r$ (Å)")
ax.set_ylabel(r"$k_{ET}$ (s$^{-1}$)")
ax.set_title(r"Dexter exchange transfer: $k_{ET}=KJ\,e^{-2r/L}$")
ax.legend()
fig.tight_layout()

# %%
# Distance at which exchange transfer still outcompetes a 10 us triplet.
for L in (1.0, 2.0, 3.0):
    r_max = 0.5 * L * np.log(K * J / 1.0e5)
    print(f"L = {L:.0f} A: triplet-triplet transfer dominates out to r = {r_max:.1f} A")

plt.show()
