r"""
Lewis and Kasha: phosphorescence as slow emission from the triplet state
==========================================================================

Lewis and Kasha (1944) showed that phosphorescence comes from a triplet
state :math:`T_1`, whose spin-forbidden return to the singlet ground state
makes it far slower than fluorescence. Using realistic rate constants
(nanosecond fluorescence, millisecond phosphorescence),
:func:`~chemistrykit.photochem.jablonski_populations_analytic` shows the
two emissions separated by six orders of magnitude in time, and
:func:`~chemistrykit.photochem.phosphorescence_quantum_yield` gives the
fraction of absorbed photons re-emitted from the triplet.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem import (
    intersystem_crossing_yield,
    jablonski_populations_analytic,
    phosphorescence_quantum_yield,
)

kf, kic, kisc = 1.0e8, 2.0e7, 3.0e8  # S1 rates, 1/s (tau_S ~ 2.4 ns)
kp, kic_T = 1.0e2, 4.0e2  # T1 rates, 1/s (tau_T = 2 ms)

t = np.logspace(-11, -1, 600)
S1, T1, _ = jablonski_populations_analytic(kf, kic, kisc, kp, kic_T, 1.0, t)
fluorescence = kf * S1  # photons emitted per second (per excited molecule)
phosphorescence = kp * T1

fig, ax = plt.subplots()
ax.loglog(t, fluorescence, label=r"fluorescence $k_f[S_1]$ (singlet)")
ax.loglog(t, phosphorescence, label=r"phosphorescence $k_p[T_1]$ (triplet)")
ax.set_ylim(1e-6, 1e9)
ax.set_xlabel("Time after excitation (s)")
ax.set_ylabel("Emission rate (photons/s per molecule)")
ax.set_title("Triplet phosphorescence outlasts singlet fluorescence")
ax.legend()
fig.tight_layout()

# %%
# The phosphorescence yield is the product of two branchings: reaching
# the triplet, then emitting from it.
phi_isc = intersystem_crossing_yield(kisc, kf, kic)
phi_p = phosphorescence_quantum_yield(kisc, kf, kic, kp, kic_T)
print(f"Singlet lifetime: {1 / (kf + kic + kisc) * 1e9:.2f} ns, triplet lifetime: {1 / (kp + kic_T) * 1e3:.1f} ms")
print(f"Phi_isc = {phi_isc:.3f}, triplet emission probability = {kp / (kp + kic_T):.3f}, Phi_p = {phi_p:.3f}")

plt.show()
