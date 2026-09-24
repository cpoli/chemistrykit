r"""
The Mark-Houwink equation: intrinsic viscosity vs. molar mass
================================================================

Mark (1938) and Houwink (1940) generalized Staudinger's linear rule to
:math:`[\eta]=KM^a`
(:func:`~chemistrykit.polymer.systems.solution_viscosity.mark_houwink_intrinsic_viscosity`),
the relation still used to turn a viscosity measurement into a molar
mass. Via Flory-Fox, :math:`[\eta]\propto R^3/M` with :math:`R\propto
M^\nu` gives :math:`a=3\nu-1`
(:func:`~chemistrykit.polymer.systems.solution_viscosity.mark_houwink_exponent_from_flory`):
:math:`a=0.5` in a theta solvent, :math:`a=0.8` in a good solvent, and
Staudinger's :math:`a=1` is recovered only for stiffer chains. Here
synthetic "measurements" with 3% seeded noise are fit back to recover
:math:`K` and :math:`a`
(:func:`~chemistrykit.polymer.systems.solution_viscosity.fit_mark_houwink`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.polymer.systems.solution_viscosity import (
    fit_mark_houwink,
    mark_houwink_exponent_from_flory,
    mark_houwink_intrinsic_viscosity,
)

rng = np.random.default_rng(0)
M = np.logspace(4, 6.5, 12)
cases = {
    "theta solvent (nu=1/2)": (0.08, mark_houwink_exponent_from_flory(0.5)),
    "good solvent (nu=3/5)": (0.01, mark_houwink_exponent_from_flory(0.6)),
}

fig, ax = plt.subplots(figsize=(6.5, 4.5))
for label, (K, a) in cases.items():
    eta = mark_houwink_intrinsic_viscosity(M, K, a) * (1 + 0.03 * rng.standard_normal(M.size))
    K_fit, a_fit = fit_mark_houwink(M, eta)
    print(f"{label}: true a = {a:.3f}, fitted a = {a_fit:.3f}; true K = {K:.3g}, fitted K = {K_fit:.3g}")
    ax.loglog(M, eta, "o", label=f"{label}: a = {a_fit:.2f}")
    ax.loglog(M, mark_houwink_intrinsic_viscosity(M, K_fit, a_fit), "k-", lw=0.8)

ax.set_xlabel("molar mass M (g/mol)")
ax.set_ylabel(r"$[\eta]$ (mL/g)")
ax.set_title(r"Mark-Houwink: $[\eta] = K M^a$")
ax.legend()
plt.tight_layout()
plt.show()
