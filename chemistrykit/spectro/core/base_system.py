r"""The shared :class:`Spectrum` result container.

:mod:`chemistrykit.spectro` mixes several genuinely different physical
models -- rigid-rotor rotational transitions, harmonic/Morse vibrational
bands, Franck-Condon vibronic progressions, and first-order NMR
multiplets -- that don't share a common polymorphic solve interface the
way :mod:`chemistrykit.quantum`'s :class:`~chemistrykit.quantum.core.base_system.QuantumSystem`/
:class:`~chemistrykit.quantum.core.base_system.VariationalSolver` do
(following :mod:`chemistrykit.solutions` and
:mod:`chemistrykit.structure`'s precedent for non-uniform model
families, each is implemented directly as classes/functions in its own
``systems/`` module). What *is* shared, and genuinely common
machinery -- the same role
:class:`chemistrykit.quantum.core.base_system.EigenstateResult` plays for
every :class:`~chemistrykit.quantum.core.base_system.VariationalSolver` --
is the output shape every one of those models eventually produces: a
discrete "stick spectrum" of peak positions and relative intensities,
optionally broadened into a continuous lineshape via
:mod:`chemistrykit.spectro.utils.lineshapes`. :class:`Spectrum` is that
stable container.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from chemistrykit.spectro.utils.lineshapes import broaden_stick_spectrum

__all__ = ["Spectrum"]


@dataclass
class Spectrum:
    """A discrete "stick" spectrum: peak positions and relative intensities.

    Parameters
    ----------
    positions : ndarray, shape (n_peaks,)
        Peak positions (e.g. wavenumber in cm^-1, or chemical shift in
        ppm -- the unit is model-specific and documented by whichever
        function constructs the `Spectrum`).
    intensities : ndarray, shape (n_peaks,)
        Relative intensity (peak area) of each position, same order as
        `positions`.
    labels : sequence of str, optional
        Human-readable label for each peak (e.g. ``"J=0->1"``).
    extra : dict, optional
        Free-form slot for additional diagnostics a model chooses to
        attach.

    Examples
    --------
    >>> import numpy as np
    >>> spectrum = Spectrum(positions=np.array([100.0, 200.0]), intensities=np.array([1.0, 0.5]))
    >>> x = np.linspace(0, 300, 3001)
    >>> broadened = spectrum.broaden(x, shape="gaussian", fwhm=5.0)
    >>> round(float(np.trapezoid(broadened, x)), 2)
    1.5
    """

    positions: np.ndarray
    intensities: np.ndarray
    labels: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    def __post_init__(self):
        self.positions = np.asarray(self.positions, dtype=np.float64)
        self.intensities = np.asarray(self.intensities, dtype=np.float64)
        if self.positions.shape != self.intensities.shape:
            raise ValueError("positions and intensities must have the same shape")

    def broaden(self, x, shape: str = "gaussian", fwhm: float = 1.0, fwhm_lorentzian=None) -> np.ndarray:
        """Broaden every stick into a continuous spectrum via :func:`chemistrykit.spectro.utils.lineshapes.broaden_stick_spectrum`.

        Parameters
        ----------
        x : array-like of float
            Grid to evaluate the broadened spectrum on.
        shape : {"gaussian", "lorentzian", "voigt"}, default "gaussian"
        fwhm : float, default 1.0
        fwhm_lorentzian : float, optional
            Required for `shape="voigt"`.

        Returns
        -------
        ndarray, shape matching `x`
        """
        return broaden_stick_spectrum(self.positions, self.intensities, x, shape=shape, fwhm=fwhm, fwhm_lorentzian=fwhm_lorentzian)

    def normalized(self) -> Spectrum:
        """Return a copy with intensities rescaled so the maximum is 1.

        Returns
        -------
        Spectrum

        Examples
        --------
        >>> spectrum = Spectrum(positions=np.array([1.0, 2.0]), intensities=np.array([4.0, 2.0]))
        >>> spectrum.normalized().intensities.tolist()
        [1.0, 0.5]
        """
        peak = np.max(np.abs(self.intensities))
        return Spectrum(positions=self.positions.copy(), intensities=self.intensities / peak, labels=list(self.labels), extra=dict(self.extra))
