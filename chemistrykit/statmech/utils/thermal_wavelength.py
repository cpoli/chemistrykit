r"""The thermal de Broglie wavelength.

.. math::

    \Lambda = \frac{h}{\sqrt{2\pi m k_BT}}

Sets the length scale below which a particle's translational motion must
be treated quantum-mechanically (its wave nature "smears out" a distance
comparable to Lambda); equivalently, twice the average thermal
wavelength of a free particle of mass `m` at temperature `T` (McQuarrie,
*Statistical Mechanics*, Ch. 4). Used by
:class:`chemistrykit.statmech.systems.partition_functions.TranslationalPartitionFunction`
and :class:`chemistrykit.statmech.systems.lattice_gas.LatticeGasAdsorption`.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import K_B, H

__all__ = ["thermal_de_broglie_wavelength"]


def thermal_de_broglie_wavelength(mass, temperature):
    r"""Return the thermal de Broglie wavelength :math:`\Lambda=h/\sqrt{2\pi mk_BT}`.

    Parameters
    ----------
    mass : float
        Particle mass, in kg.
    temperature : float
        Absolute temperature, in K.

    Returns
    -------
    float
        Wavelength, in m.

    Examples
    --------
    For argon at room temperature, Lambda is a small fraction of an
    angstrom -- far shorter than the mean interparticle spacing in a
    gas at atmospheric pressure, confirming that treating translational
    motion classically (as the ideal gas law does) is an excellent
    approximation there:

    >>> import scipy.constants as sc
    >>> wavelength = thermal_de_broglie_wavelength(mass=39.948 * sc.atomic_mass, temperature=298.15)
    >>> bool(wavelength < 1e-10)
    True
    """
    return H / np.sqrt(2.0 * np.pi * mass * K_B * temperature)
