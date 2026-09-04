"""The package shim, the version, and the typing stubs."""

import pathlib
from importlib import metadata

import sylvester
import sylvester._sylvester as extension


def test_the_version_matches_the_distribution():
    assert sylvester.__version__ == metadata.version("sylvester")


def test_the_shim_and_the_extension_export_the_same_names():
    exported = {name for name in dir(extension) if not name.startswith("_")}
    assert set(sylvester.__all__) == exported | {"__version__"}


def test_every_exported_name_is_in_the_stub():
    stub = pathlib.Path(sylvester.__file__).with_name("__init__.pyi").read_text()
    for name in sylvester.__all__:
        assert name in stub, name


def test_the_package_ships_the_typing_marker():
    assert pathlib.Path(sylvester.__file__).with_name("py.typed").is_file()


def test_every_class_names_the_package_it_lives_in():
    for name in sylvester.__all__:
        value = getattr(sylvester, name)
        if isinstance(value, type):
            assert value.__module__ == "sylvester", name
