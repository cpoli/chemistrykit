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


def test_thermo_subpackage_is_declared_in_all():
    assert "thermo" in ck.__all__


def test_solutions_subpackage_is_declared_in_all():
    assert "solutions" in ck.__all__


def test_md_subpackage_is_declared_in_all():
    assert "md" in ck.__all__


def test_statmech_subpackage_is_declared_in_all():
    assert "statmech" in ck.__all__


def test_quantum_subpackage_is_declared_in_all():
    assert "quantum" in ck.__all__


def test_spectro_subpackage_is_declared_in_all():
    assert "spectro" in ck.__all__


def test_structure_subpackage_is_declared_in_all():
    assert "structure" in ck.__all__


def test_domain_subpackage_versions_match_top_level_version():
    assert ck.thermo.__version__ == ck.__version__
    assert ck.solutions.__version__ == ck.__version__
    assert ck.md.__version__ == ck.__version__
    assert ck.statmech.__version__ == ck.__version__
    assert ck.quantum.__version__ == ck.__version__
    assert ck.spectro.__version__ == ck.__version__
    assert ck.structure.__version__ == ck.__version__
