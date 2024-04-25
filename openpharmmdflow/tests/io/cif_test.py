from importlib import resources

import pytest
from openff.toolkit.utils.exceptions import MoleculeParseError

from openpharmmdflow import from_cif

IO_DATA_PATH = resources.files("openpharmmdflow.tests.io.data")


def AWUDEB_cif():
    return IO_DATA_PATH / "AWUDEB.cif"


def HXACAN_cif():
    return IO_DATA_PATH / "HXACAN.cif"


def TRICKY_cif():
    # From https://www.crystallography.net/cod/2102215.html
    return IO_DATA_PATH / "2102215.cif"


@pytest.mark.parametrize(
    "cif_file, should_load",
    [(HXACAN_cif(), True), (AWUDEB_cif(), True), (TRICKY_cif(), False)],
)
def test_cif_loading_defaults(cif_file, should_load):
    print(cif_file)
    if should_load:
        from_cif(cif_file)
    else:
        with pytest.raises(MoleculeParseError):
            from_cif(cif_file)
