r"""
The Daniell cell: a steady 1.10 V from zinc and copper
========================================================

Daniell's 1836 two-fluid cell separated a zinc electrode in zinc (or acid)
solution from a copper electrode in copper sulfate. Copper metal, not
hydrogen gas, is deposited at the cathode, so the cell does not polarize.
This example reproduces its standard potential with
:func:`~chemistrykit.electrochem.systems.standard_potentials.standard_cell_potential`,
balances the overall reaction with
:func:`~chemistrykit.electrochem.systems.standard_potentials.balance_redox_reaction`,
and uses the Nernst equation to show how gently the voltage falls as the
cell discharges -- the steadiness that made it the first practical voltage
standard.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.nernst import nernst_potential
from chemistrykit.electrochem.systems.standard_potentials import STANDARD_REDUCTION_POTENTIALS as T
from chemistrykit.electrochem.systems.standard_potentials import balance_redox_reaction, standard_cell_potential

# %%
# Zn(s) + Cu2+(aq) -> Zn2+(aq) + Cu(s)
E0 = standard_cell_potential("Cu2+/Cu", "Zn2+/Zn")
cathode_mult, anode_mult, n = balance_redox_reaction(T["Cu2+/Cu"], T["Zn2+/Zn"])
print(f"Daniell cell E0 = {E0:.2f} V; balancing multiples = ({cathode_mult}, {anode_mult}), n = {n}")

# %%
# As the cell discharges, Cu2+ is consumed and Zn2+ produced
# (starting from 1 M each). Q = [Zn2+]/[Cu2+].
extent = np.linspace(0.0, 0.99, 200)
Q = (1.0 + extent) / (1.0 - extent)
E = nernst_potential(E0, n, Q)

fig, ax = plt.subplots()
ax.plot(extent * 100, E)
ax.axhline(E0, color="gray", linestyle="--", label=r"$E^\circ$ = 1.10 V")
ax.set_xlabel(r"Cu$^{2+}$ consumed (%)")
ax.set_ylabel("Cell voltage (V)")
ax.set_ylim(0.9, 1.2)
ax.set_title("Daniell cell voltage during discharge")
ax.legend()
fig.tight_layout()
plt.show()
