# Contributing to chemistrykit

Thanks for considering a contribution. chemistrykit is organized as one
subpackage per chemistry domain (`chemistrykit/<name>/`), each with its own
`systems/` (model classes), `core/` (ABCs and, where needed, numba-
accelerated kernels), `utils/` (supporting numerics), `visualizers/`
(matplotlib/plotly plotting), and `tests/` directory. New chemistry belongs
in the subpackage it fits best; a genuinely new domain gets its own
subpackage following the same layout.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install
```

## Before opening a PR

```bash
ruff check .                 # lint
ruff format .                # format
MPLBACKEND=Agg pytest -q     # unit tests
MPLBACKEND=Agg pytest --doctest-modules chemistrykit \
    --ignore-glob="*/tests/*"   # docstring examples
```

All three run in CI (`.github/workflows/ci.yml`) on every PR, across
Python 3.10-3.12 on Linux and macOS. `mypy` also runs in CI but is
currently advisory (non-blocking) — see "Type checking" below.

If you touch anything under `docs/` or add/modify an example in
`examples/`, also build the docs locally before opening a PR (this
re-executes every changed `examples/*/plot_*.py` script via
sphinx-gallery, which is the only way to catch a broken example):

```bash
pip install -e ".[docs]"
cd docs && make html
```

## Code style

- Follow the existing style in the subpackage you're editing — chemistry
  notation (`k`, `A`, `B`, `T`, single-letter variables) is intentional and
  allowed (see the `ruff` ignores in `pyproject.toml` for the specific,
  deliberate exceptions).
- Every public function/class gets a NumPy-style docstring with a
  `Parameters`, `Returns`, and (where it adds real value beyond what the
  signature already says) an `Examples` section with a runnable doctest.
  Doctests are checked in CI — a docstring example that doesn't actually
  run is worse than no example.
- Prefer closed-form / analytically-verifiable results in tests
  (`pytest.approx` against a known formula) over snapshot-testing plot
  output.
- New chemistry should cite where the formula comes from (a textbook
  name/edition or a paper), either in the docstring or as a comment, so a
  reader can verify it independently.

## Type checking

`mypy` is configured in `pyproject.toml` but not yet fully clean across
the codebase (mostly matplotlib/numpy stub gaps around animation objects
and array-typed arguments) — it currently runs in CI as an advisory,
non-blocking job. New code should type-check cleanly where practical;
fixing pre-existing errors in code you're not otherwise touching is
welcome but not required.

## Tests

Tests live alongside each subpackage (`chemistrykit/<name>/tests/`), not in
a top-level `tests/` directory. A new model needs:

- At least one test that checks a closed-form / analytically-known result,
  not just "it runs without crashing."
- A visualizer needs only a smoke test (it returns the right type/shape,
  and — for animations — that `anim.save()` to a temp file succeeds).

## History entries

`docs/source/history/*_breakthroughs.rst` is a curated chronology per
domain, not a general-purpose list of "interesting chemistry history."
Each entry exists because it connects to something this package actually
implements. Keep it that way:

- A new entry must cite the concrete class/method it connects to
  (`:meth:`.../:class:`...` cross-reference) and belongs in the same PR
  as the implementation it describes, not added on its own.
- One entry per PR. If you're adding several related pieces of
  chemistry, open separate PRs (or ask first) rather than batching a set
  of history entries together.
- Changes to `docs/source/history/**` require a maintainer review
  (enforced via `.github/CODEOWNERS`) even if the rest of the PR is
  otherwise approved.

## Reporting bugs / requesting features

Open a GitHub issue. For a chemistry bug specifically, include the formula
or reference you expected the code to match, and (if possible) the
numeric discrepancy — "this doesn't look right" is much slower to act on
than "the half-life should be within a few % of the textbook worked
example and I'm getting X."

## Code of Conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md).
