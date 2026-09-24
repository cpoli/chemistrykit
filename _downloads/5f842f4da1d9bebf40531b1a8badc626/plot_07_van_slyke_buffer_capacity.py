r"""
Van Slyke's buffer capacity: how strongly a solution resists pH change
========================================================================

Donald Van Slyke defined the buffer value :math:`\beta = dC_b/d\mathrm{pH}`,
the strong base per litre needed to raise the pH by one unit. For a weak
acid/conjugate base pair of total concentration :math:`C`
(:func:`~chemistrykit.solutions.systems.acid_base.buffer_capacity`):

.. math::

   \beta = \ln 10\left([H^+] + \frac{K_w}{[H^+]} + \frac{C K_a [H^+]}{(K_a+[H^+])^2}\right).

The buffer term peaks at :math:`\mathrm{pH} = pK_a` with height
:math:`\ln 10\,C/4`, while water alone buffers only at the extremes of
the pH scale.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.acid_base import buffer_capacity

pKa = 4.76
pH = np.linspace(1, 13, 600)

fig, ax = plt.subplots(figsize=(7, 5))
for C, color in [(0.01, "seagreen"), (0.05, "darkorange"), (0.10, "steelblue")]:
    ax.plot(pH, buffer_capacity(pH, C=C, Ka=10**-pKa), color=color, label=f"acetate, C = {C} M")
ax.plot(pH, buffer_capacity(pH, C=0.0, Ka=10**-pKa), color="gray", linestyle="--", label="water only (C = 0)")
ax.axvline(pKa, color="gray", linestyle=":", linewidth=0.8)
ax.set_ylim(0, 0.08)
ax.set_xlabel("pH")
ax.set_ylabel(r"buffer capacity $\beta$ (mol/L per pH unit)")
ax.set_title("Van Slyke buffer capacity of acetate buffers")
ax.legend(fontsize=8)
fig.tight_layout()

# %%
# Checking the peak against :math:`\ln 10\,C/4`:

for C in (0.01, 0.05, 0.10):
    print(f"C = {C:.2f} M: beta(pKa) = {float(buffer_capacity(pKa, C, 10**-pKa)):.5f}, ln10*C/4 = {np.log(10) * C / 4:.5f}")

plt.show()
