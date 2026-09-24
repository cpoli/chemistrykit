"""chemistrykit.spectro: spectroscopy.

The Beer-Lambert absorbance law and one well-characterized instrumental
deviation from its linearity (stray light); rigid-rotor rotational
spectra with :math:`\\Delta J=\\pm1` selection rules and isotope shifts,
built on :mod:`chemistrykit.quantum.systems.rigid_rotor`; harmonic vs.
Morse (anharmonic) vibrational band positions, built on
:mod:`chemistrykit.quantum.systems.harmonic_oscillator`, plus a genuine
triatomic normal-mode calculation via the Wilson GF-matrix method;
Franck-Condon vibronic progressions for electronic (UV-Vis) spectra; a
first-order NMR multiplet simulator (chemical shifts and J-coupling
splitting patterns), Larmor frequencies, the Karplus relation, and
Fourier-transform NMR (free-induction decay to spectrum); the
Balmer-Rydberg formula for hydrogen-like atomic lines; and Lorentzian/Gaussian/Voigt lineshape utilities
shared across all of the above.
"""

__version__ = "0.1.0"

from chemistrykit.spectro.core.base_system import Spectrum
from chemistrykit.spectro.systems.atomic import rydberg_wavenumber
from chemistrykit.spectro.systems.beer_lambert import (
    absorbance,
    apparent_absorbance_with_stray_light,
    concentration_from_absorbance,
    transmittance,
)
from chemistrykit.spectro.systems.electronic import (
    franck_condon_factor,
    franck_condon_progression,
    franck_condon_spectrum,
    huang_rhys_factor,
)
from chemistrykit.spectro.systems.nmr import (
    chemical_shift_ppm,
    fid_to_spectrum,
    first_order_multiplet,
    free_induction_decay,
    karplus_coupling,
    larmor_frequency,
    multi_coupling_multiplet,
    multiplicity,
    pascals_triangle_intensities,
)
from chemistrykit.spectro.systems.rotational import energy_to_wavenumber, isotope_shift_ratio, rotational_line_wavenumbers, rotational_spectrum
from chemistrykit.spectro.systems.vibrational import (
    NormalModeResult,
    TriatomicNormalModes,
    anharmonicity_from_overtones,
    harmonic_fundamental_wavenumber,
    morse_transition_wavenumbers,
)
from chemistrykit.spectro.utils.lineshapes import broaden_stick_spectrum, gaussian, lorentzian, voigt

__all__ = [
    "__version__",
    "Spectrum",
    "absorbance",
    "transmittance",
    "concentration_from_absorbance",
    "apparent_absorbance_with_stray_light",
    "energy_to_wavenumber",
    "rotational_line_wavenumbers",
    "rotational_spectrum",
    "isotope_shift_ratio",
    "harmonic_fundamental_wavenumber",
    "morse_transition_wavenumbers",
    "anharmonicity_from_overtones",
    "TriatomicNormalModes",
    "NormalModeResult",
    "huang_rhys_factor",
    "franck_condon_factor",
    "franck_condon_progression",
    "franck_condon_spectrum",
    "multiplicity",
    "pascals_triangle_intensities",
    "first_order_multiplet",
    "multi_coupling_multiplet",
    "larmor_frequency",
    "chemical_shift_ppm",
    "karplus_coupling",
    "free_induction_decay",
    "fid_to_spectrum",
    "rydberg_wavenumber",
    "gaussian",
    "lorentzian",
    "voigt",
    "broaden_stick_spectrum",
]
