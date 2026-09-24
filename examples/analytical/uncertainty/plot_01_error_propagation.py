r"""
Propagation of uncertainty with the first-order (Ku) formula
===============================================================

The closed-form rules for sums (:func:`~chemistrykit.analytical.systems.uncertainty.propagate_sum`),
products/quotients (:func:`~chemistrykit.analytical.systems.uncertainty.propagate_product`),
and powers (:func:`~chemistrykit.analytical.systems.uncertainty.propagate_power`)
are each a special case of the general first-order propagation formula,
which :func:`~chemistrykit.analytical.systems.uncertainty.propagate_uncertainty`
evaluates directly via numerical partial derivatives -- useful when no
simple closed form is at hand.
"""

# %%
from chemistrykit.analytical.systems.uncertainty import (
    propagate_power,
    propagate_product,
    propagate_sum,
    propagate_uncertainty,
)

# %%
# A titration's net volume, Vb_final - Vb_initial, both read from a buret
# with the same reading uncertainty:
sigma_reading = 0.02  # mL
sigma_volume = propagate_sum([sigma_reading, sigma_reading])
print(f"Net volume uncertainty (two buret readings of {sigma_reading} mL each): {sigma_volume:.4f} mL")

# %%
# Molarity from mass, molar mass, and volume: M = m / (MW * V) -- a
# product/quotient of three measured quantities.
mass, sigma_mass = 0.2500, 0.0002  # g
molar_mass, sigma_molar_mass = 58.44, 0.01  # g/mol, NaCl
volume, sigma_volume_L = 0.1000, 0.0002  # L

moles = mass / molar_mass
molarity = moles / volume
sigma_molarity = propagate_product(
    [mass, 1.0 / molar_mass, 1.0 / volume],
    [sigma_mass, sigma_molar_mass / molar_mass**2, sigma_volume_L / volume**2],
)
print(f"\nMolarity = {molarity:.6f} +/- {sigma_molarity:.6f} mol/L")

# %%
# A cell's volume from a measured edge length, V = L^3:
L, sigma_L = 2.000, 0.005  # cm
V = L**3
sigma_V = propagate_power(L, sigma_L, 3.0)
print(f"\nVolume = {V:.4f} +/- {sigma_V:.4f} cm^3 (relative uncertainty tripled: {sigma_V / V:.4%} vs {sigma_L / L:.4%})")


# %%
# General numerical propagation reproduces the molarity closed-form
# result for an arbitrary function of the same three variables:
def molarity_func(m, mw, v):
    return m / (mw * v)


sigma_molarity_numeric = propagate_uncertainty(molarity_func, [mass, molar_mass, volume], [sigma_mass, sigma_molar_mass, sigma_volume_L])
print(f"\nMolarity uncertainty (closed form):  {sigma_molarity:.6f}")
print(f"Molarity uncertainty (general/numeric): {sigma_molarity_numeric:.6f}")
