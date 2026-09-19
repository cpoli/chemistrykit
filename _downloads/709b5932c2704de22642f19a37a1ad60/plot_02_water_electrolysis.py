r"""
Nicholson and Carlisle's electrolysis of water
==================================================

Days after learning of Volta's newly announced pile, William Nicholson
and Anthony Carlisle connected one to a bowl of water in May 1800 and
watched gas bubble from each wire -- the first electrolysis, and the
first new discovery the pile made possible. This example reconstructs
both halves of that observation from functions already in this package:
the theoretical minimum decomposition voltage of water
(:func:`~chemistrykit.electrochem.systems.electrolysis.minimum_applied_voltage_electrolytic`,
applied to the tabulated H+/H2 and O2/H2O standard potentials) comes out
at the textbook 1.23 V, and Faraday's law
(:func:`~chemistrykit.electrochem.systems.electrolysis.moles_from_charge`)
recovers the 2:1 hydrogen-to-oxygen ratio Nicholson and Carlisle actually
saw bubbling from the two wires, directly from electron-transfer
stoichiometry.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.electrolysis import (
    charge_from_current,
    minimum_applied_voltage_electrolytic,
    moles_from_charge,
)
from chemistrykit.electrochem.systems.standard_potentials import (
    STANDARD_REDUCTION_POTENTIALS,
    cell_potential,
    is_spontaneous,
)

# %%
# The overall water-splitting reaction, 2 H2O -> 2 H2 + O2, splits into a
# cathodic reduction (2 H+ + 2 e- -> H2) and an anodic oxidation
# (2 H2O -> O2 + 4 H+ + 4 e-, the reverse of the tabulated O2/H2O
# reduction) -- both half-reactions are already in this package's
# standard-potential table.
h2 = STANDARD_REDUCTION_POTENTIALS["H+/H2"]
o2 = STANDARD_REDUCTION_POTENTIALS["O2/H2O"]
E_water = cell_potential(cathode=h2, anode=o2)
V_min = minimum_applied_voltage_electrolytic(E_water)
print(f"Water-splitting cell potential (cathode H+/H2, anode O2/H2O): {E_water:.2f} V")
print(f"Spontaneous as written: {is_spontaneous(E_water)}")
print(f"Minimum applied voltage to force electrolysis: {V_min:.2f} V")
print("(the accepted standard decomposition voltage of water is 1.23 V)")

# %%
# Faraday's law applied to the two evolved gases: the same charge passed
# produces H2 (n=2 electrons per molecule) at twice the molar rate of O2
# (n=4 electrons per molecule) -- exactly the 2:1 volume ratio Nicholson
# and Carlisle observed bubbling from the two wires (equal moles of gas
# occupy equal volume, by Avogadro's law).
current, time_s = 0.5, 600.0  # A, s
Q = charge_from_current(current, time_s)
moles_H2 = moles_from_charge(Q, n=2)
moles_O2 = moles_from_charge(Q, n=4)
print(f"\nCharge passed in {time_s:.0f} s at {current} A: {Q:.1f} C")
print(f"Moles H2 evolved: {moles_H2:.3e} mol")
print(f"Moles O2 evolved: {moles_O2:.3e} mol")
print(f"H2:O2 mole (and volume) ratio: {moles_H2 / moles_O2:.2f}")

# %%
fig, ax = plt.subplots()
times = np.linspace(0.0, time_s, 100)
Q_t = charge_from_current(current, times)
ax.plot(times, moles_from_charge(Q_t, n=2) * 1e3, label="H2")
ax.plot(times, moles_from_charge(Q_t, n=4) * 1e3, label="O2")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Gas evolved (mmol)")
ax.set_title("Faraday's law: gas evolution during water electrolysis")
ax.legend()
fig.tight_layout()
plt.show()
