r"""
The equipartition theorem, and why only one mode actually "freezes out"
===========================================================================

James Clerk Maxwell (1860) and, in fully general form, Ludwig Boltzmann
(1871) showed that every classical quadratic degree of freedom of a
system in thermal equilibrium contributes exactly :math:`\frac12k_B` to
its heat capacity -- the equipartition theorem. A rigid diatomic molecule
therefore "should" carry :math:`\frac32k_B` (translation) + :math:`k_B`
(2 rotational degrees of freedom) + :math:`k_B` (1 vibrational mode,
kinetic + potential) = :math:`\frac72k_B` at every temperature -- a
prediction 19th-century equipartition theory could not explain deviations
from (Maxwell himself flagged the puzzle sharply in 1879). Only
:class:`~chemistrykit.statmech.VibrationalPartitionFunctionHarmonic` is an
*exact* quantum treatment in this package, so it alone shows the real
resolution Einstein gave in 1907: a mode's contribution only "switches on"
once :math:`k_BT` exceeds its own quantum spacing.
:class:`~chemistrykit.statmech.RotationalPartitionFunctionLinear`, by
contrast, implements only the classical high-temperature limit (as its own
docstring flags), so it contributes its full :math:`k_B` at *every*
temperature here -- a deliberate approximation, harmless in practice
because real rotational spacings are tiny compared to vibrational ones (as
the numbers below confirm), but worth seeing explicitly rather than
assuming.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.statmech import IdealGasMolecule, VibrationalPartitionFunctionHarmonic

K_B = 1.380649e-23

# HCl-like molecule: moment of inertia and vibrational wavenumber matched
# to plot_01_vibrational_heat_capacity.py.
vib_mode = VibrationalPartitionFunctionHarmonic.from_wavenumber(2886.0)
hcl = IdealGasMolecule(
    mass=6.15e-26,
    volume=1.0e-3,
    moment_of_inertia=1.45e-46,
    vibrational_frequencies=[vib_mode.frequency],
)
theta_rot = hcl.rotational.rotational_temperature
theta_vib = vib_mode.vibrational_temperature
print(f"Rotational temperature Theta_rot = {theta_rot:.2f} K")
print(f"Vibrational temperature Theta_vib = {theta_vib:.1f} K")
print(f"Theta_vib / Theta_rot = {theta_vib / theta_rot:.0f}x higher -- why rotation looks 'always classical' here")

# %%
T = np.logspace(0.0, 4.3, 400)
Cv_total = np.array([hcl.heat_capacity_v(t, N=1.0) for t in T]) / K_B
Cv_rot = np.array([hcl.rotational.heat_capacity_v(t, N=1.0) for t in T]) / K_B
Cv_vib = np.array([hcl.vibrational_modes[0].heat_capacity_v(t, N=1.0) for t in T]) / K_B

fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogx(T, Cv_total, color="steelblue", label="total Cv/k_B")
ax.semilogx(T, 1.5 + Cv_rot, color="gray", linestyle="--", label="translation + rotation only")
for value in (1.5, 2.5, 3.5):
    ax.axhline(value, linestyle=":", linewidth=0.6, color="gray")
ax.axvline(theta_vib, color="darkorange", linestyle="--", linewidth=0.8, label="Theta_vib")
ax.set_xlabel("T (K)")
ax.set_ylabel("Cv / k_B (per molecule)")
ax.set_title("HCl-like diatomic: only the vibrational mode actually 'unfreezes'")
ax.legend()
fig.tight_layout()

# %%
# The classical rotational contribution never moves (it is pinned at
# exactly 1.0 k_B, its high-temperature limit, at every T shown), while
# the exact quantum vibrational contribution climbs from essentially zero
# to its own classical limit only once T approaches Theta_vib -- the
# single genuine step in the total curve, and the mode Einstein's 1907
# theory was built to explain:

for T_check, label in [(10.0, "T << Theta_vib"), (300.0, "T ~ room temperature"), (5000.0, "T ~ Theta_vib")]:
    Cv = hcl.heat_capacity_v(T_check, N=1.0) / K_B
    Cv_v = hcl.vibrational_modes[0].heat_capacity_v(T_check, N=1.0) / K_B
    print(f"{label:20s} (T={T_check:>7.1f} K): Cv_total/k_B = {Cv:.3f}  (Cv_vib/k_B = {Cv_v:.4f})")

plt.show()
