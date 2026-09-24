r"""
The Lennard-Jones 12-6 potential: repulsion, dispersion, and the minimum
===========================================================================

Lennard-Jones's 1924 pair potential combines a steep :math:`r^{-12}`
repulsion with the :math:`r^{-6}` London-dispersion attraction:

.. math::

   U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right].

:class:`~chemistrykit.md.systems.lj_fluid.LennardJones` evaluates it and
its force. The left panel splits the potential into its two terms; the
right panel shows the potential and force together, crossing zero at
:math:`r=\sigma` and :math:`r=2^{1/6}\sigma`
(:attr:`~chemistrykit.md.systems.lj_fluid.LennardJones.r_min`)
respectively. The curves are drawn for argon, with
:math:`\epsilon/k_B = 120\ \text{K}` and :math:`\sigma = 3.4\ \text{\AA}`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.md.systems.lj_fluid import LennardJones

epsilon_K, sigma_A = 120.0, 3.4  # argon: epsilon/k_B in kelvin, sigma in angstrom
lj = LennardJones(epsilon=epsilon_K, sigma=sigma_A)
r = np.linspace(3.0, 9.0, 500)

repulsion = 4.0 * epsilon_K * (sigma_A / r) ** 12
attraction = -4.0 * epsilon_K * (sigma_A / r) ** 6

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
axes[0].plot(r, repulsion, color="firebrick", linestyle="--", label=r"repulsion $+4\epsilon(\sigma/r)^{12}$")
axes[0].plot(r, attraction, color="steelblue", linestyle="--", label=r"dispersion $-4\epsilon(\sigma/r)^{6}$")
axes[0].plot(r, lj.energy(r), color="black", label="sum U(r)")
axes[0].axhline(0.0, color="gray", linewidth=0.8)
axes[0].set_ylim(-250, 400)
axes[0].set_xlabel(r"r ($\mathrm{\AA}$)")
axes[0].set_ylabel(r"$U/k_B$ (K)")
axes[0].set_title("Two terms of the 12-6 potential")
axes[0].legend()

axes[1].plot(r, lj.energy(r), color="black", label="U(r)")
axes[1].plot(r, lj.force_scalar(r), color="darkorange", label="force -dU/dr")
axes[1].axhline(0.0, color="gray", linewidth=0.8)
axes[1].axvline(sigma_A, color="gray", linestyle=":", linewidth=0.8)
axes[1].axvline(lj.r_min, color="gray", linestyle="--", linewidth=0.8)
axes[1].annotate(r"$r=\sigma$", (sigma_A, 150), xytext=(3.45, 150))
axes[1].annotate(r"$r_{min}=2^{1/6}\sigma$", (lj.r_min, -200), xytext=(lj.r_min + 0.1, -200))
axes[1].set_ylim(-250, 250)
axes[1].set_xlabel(r"r ($\mathrm{\AA}$)")
axes[1].set_title("Potential and force")
axes[1].legend()
fig.tight_layout()

# %%
# The well depth equals epsilon and sits at :math:`2^{1/6}\sigma`, where
# the force vanishes:

print(f"r_min = {lj.r_min:.3f} angstrom, U(r_min)/k_B = {float(lj.energy(lj.r_min)):.1f} K")
print(f"U(sigma) = {float(lj.energy(sigma_A)):.1e} K, force at r_min = {float(lj.force_scalar(lj.r_min)):.1e}")

plt.show()
