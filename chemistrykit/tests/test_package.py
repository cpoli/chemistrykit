"""Package-level smoke tests: import surface and version consistency."""

import chemistrykit as ck


def test_top_level_imports_expose_declared_subpackages():
    for name in ck.__all__:
        assert hasattr(ck, name), f"chemistrykit.{name} is in __all__ but not importable"


def test_version_is_a_string_and_matches_kinetics_version():
    assert isinstance(ck.__version__, str)
    assert ck.__version__ == ck.kinetics.__version__


def test_constants_and_integrators_are_modules():
    import types

    assert isinstance(ck.constants, types.ModuleType)
    assert isinstance(ck.integrators, types.ModuleType)


def test_kinetics_subpackage_is_declared_in_all():
    assert "kinetics" in ck.__all__
