# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A unified numerical toolkit for computational chemistry: 14 domain subpackages (`analytical`, `crystal`, `electrochem`, `kinetics`, `md`, `photochem`, `polymer`, `quantum`, `solutions`, `spectro`, `statmech`, `structure`, `surface`, `thermo`) sharing common ODE integrators, chemical constants, periodic-table data, and a consistent NumPy-based API. Conventionally imported as `ck`.

chemistrykit is a sibling package to `../physicskit`, cloning its architecture and engineering conventions exactly and replacing physics content with chemistry content — every model is implemented from first-principles formulas and small, self-contained numerics (no RDKit, ASE, PySCF, OpenMM, or SMILES/PDB parsing dependency), exactly how physicskit implements e.g. tight-binding models or N-body dynamics by hand rather than importing a materials-science package. When a convention here seems underspecified, check how physicskit (the literal template) handles the equivalent case, rather than inventing a new one.

## Commands

```bash
# setup
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install

# lint / format
ruff check .
ruff format .

# unit tests -- everything, or a single subpackage / file / test
MPLBACKEND=Agg pytest -q
pytest chemistrykit/kinetics/tests
pytest chemistrykit/kinetics/tests/test_rate_laws.py::test_first_order_matches_exponential_decay -q

# docstring examples (every subpackage's doctests)
MPLBACKEND=Agg pytest --doctest-modules chemistrykit --ignore-glob="*/tests/*"

# type check (advisory only in CI, not blocking -- see below)
mypy chemistrykit

# docs (re-executes examples/*/plot_*.py via sphinx-gallery -- the only way
# to catch a broken example; do this if you touch docs/ or examples/)
pip install -e ".[docs]"
cd docs && make html
```

All of lint, unit tests, doctests, and (advisory) mypy run in CI (`.github/workflows/ci.yml`) across Python 3.10-3.12 on Linux and macOS. Unlike physicskit, `pytest-xdist` is not a test dependency here, so there's no `-n auto` parallel flag -- just `pytest -q`.

## Architecture

**One subpackage per chemistry domain**, each living at `chemistrykit/<name>/` with its own `tests/` directory (`chemistrykit/<name>/tests/`, not a top-level `tests/`). New chemistry belongs in the subpackage it fits best; a genuinely new domain gets its own subpackage. Every subpackage follows the same internal layout: `core/` (`base_system.py` — ABCs, result dataclasses, and any shared `@njit` kernels), `systems/` (concrete models, one module per family), `utils/` (domain-specific numerics supporting `systems/` but not a model itself), `visualizers/` (matplotlib/plotly plotting), `tests/`. Read the target subpackage's own `__init__.py` and `core/base_system.py` before assuming further detail.

**The numba first-class-function pattern**, documented in physicskit's `classical/core/base_system.py` and referenced directly from this codebase's own docstrings (e.g. `chemistrykit/kinetics/core/base_system.py`), applies where used — but only in the two subpackages with genuinely performance-critical inner loops: `chemistrykit.kinetics` (integrating reaction networks via `chemistrykit.integrators`) and `chemistrykit.md` (pairwise force evaluation in molecular dynamics). A module-level factory function closes over a system's numeric parameters and returns a standalone `@njit` dispatcher (conventionally `self._rhs_njit`); base classes always integrate using that attribute, never a bound Python method, since Numba's nopython mode can only call a genuine njit dispatcher from inside another njit function. Most other subpackages have no `@njit` code at all — don't add it speculatively.

`chemistrykit.integrators` (top-level, shared) mirrors physicskit's — `rk4_integrate`, `leapfrog_integrate`/`velocity_verlet_integrate`, `yoshida4_integrate`, adaptive `dopri5_integrate` — with the same `f(state_or_pos, t, params) -> ndarray` calling convention, letting one compiled callback be reused across systems with different parameter values. `chemistrykit.kinetics.core.integrators` and `chemistrykit.md.core` wrap this shared module rather than reimplementing per-domain integrators, exactly as physicskit's `classical/core/integrators.py` does for its own njit-callback calling convention.

`chemistrykit.constants` is the single source of truth for *actual SI values* (from `scipy.constants`, CODATA/2019 SI): `R`, `NA`, `K_B`, `H`, `HBAR`, `C`, `ELEMENTARY_CHARGE`, `FARADAY`, `ATOMIC_MASS_UNIT`, standard pressure/temperature. Unlike physicskit, **no chemistrykit subpackage adopts a rescaled unit system** — there is no chemistry analogue of setting ħ=1 — so every subpackage computes with the constant's actual numeric value rather than a domain-specific natural-unit convention. `chemistrykit.periodic_table` (re-exported from `constants`) provides a small built-in periodic table (atomic number, symbol, atomic mass) as plain data, not a cheminformatics dependency. Only reach for a `constants` conversion helper (energy↔temperature via `K_B`, atm↔Pa) where one is explicitly provided; a conversion that secretly requires picking a solvent/activity-coefficient model (e.g. concentration→mole fraction) has no generic helper by design — see the relevant subpackage (`chemistrykit.solutions`, `chemistrykit.thermo`) instead.

Where subpackages share base classes, the common flow is: `core/base_system.py` defines the domain's ABC(s) — e.g. `chemistrykit.kinetics` splits into `RateLaw` (a single elementary reaction with a closed-form integrated profile, no numerical integration needed) and `ReactionNetwork` (coupled species integrated numerically), mirroring physicskit's `chaos.core.base_system.DynamicalSystem` (a bare vector field) rather than `classical`'s Hamiltonian/Lagrangian hierarchy, which has no chemistry analogue — plus a small `@dataclass` result container per model family (mirroring physicskit's `SimulationResult` pattern — never bare tuples, so visualizers have a stable interface); concrete models in `systems/` implement the ABC; `visualizers/` consumes the result dataclass for plotting.

## Conventions

- Chemistry single-letter variable names (`t`, `A`, `B`, `k`, etc.) are intentional and preserved — see the deliberate `ruff` ignores in `pyproject.toml` (`E741` and others) rather than "fixing" them.
- Every public function/class needs a NumPy-style docstring (`Parameters`, `Returns`, and an `Examples` section with a runnable doctest where it adds real value). Doctests are checked in CI — an example that doesn't actually execute correctly is worse than no example.
- New models should cite their source formula (docstring or comment) so it can be independently verified.
- Tests prefer closed-form/analytically-verifiable assertions (`pytest.approx` against a known formula, e.g. first-order kinetics matching the exponential-decay solution) over snapshot-testing plot output; a visualizer needs only a smoke test (right return type/shape; for animations, that `anim.save()` to a temp file succeeds).
- `mypy` is configured but not yet fully clean and runs advisory/non-blocking in CI. New code should type-check where practical; fixing unrelated pre-existing errors is not required.
- `docs/source/history/` documents each subpackage's foundational chemistry breakthroughs linked to the corresponding implementation — worth checking when adding a major new model to understand the expected historical framing.
- Every breakthrough in `docs/source/history/` links to one or more gallery examples via `.. minigallery::`, and no example is linked from two breakthroughs. Each example's title and content must show that breakthrough's subject, so a reader can tell at a glance why it illustrates that entry. Split an example that covers several breakthroughs into focused ones, one per breakthrough, and delete the combined file. When the match is unclear, change the example, not the breakthrough's title or text.

# Output Guidelines
- **Be Concise:** Provide direct code and answers first. Omit setup text, fluff, and conversational responses.
- **Code Generation:** Output only updated code blocks, functions, or unified diffs. Never output full files unless instructed.
- **Explanations:** Limit inline code comments. Explain high-level logic in 1–2 bullet points after the code block.
- **Format:** Use bulleted lists and tables for comparisons instead of paragraph blocks.
