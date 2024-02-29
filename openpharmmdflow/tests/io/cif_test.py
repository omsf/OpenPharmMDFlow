from importlib import resources

import pytest

from openpharmmdflow import from_cif


@pytest.fixture()
def AWUDEB_cif():
    with resources.files("openpharmmdflow.tests.io.data") as d:
        yield str(d / "AWUDEB.cif")


@pytest.fixture()
def HXACAN_cif():
    with resources.files("openpharmmdflow.tests.io.data") as d:
        yield str(d / "HXACAN.cif")


@pytest.mark.parametrize(
    "cif_file, should_load", [(HXACAN_cif, True), (AWUDEB_cif, False)]
)
def test_cif_loading_defaults(cif_file, should_load):
    if should_load:
        from_cif(cif_file)
    else:
        with pytest.raises(RuntimeError):
            from_cif(cif_file)
