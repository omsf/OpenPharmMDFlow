import sys
from importlib import resources

import pytest
from openff.toolkit.utils.exceptions import MoleculeParseError

from openpharmmdflow import from_cif
from openpharmmdflow.io.load import load_file

IO_DATA_PATH = resources.files("openpharmmdflow.tests.io.data")


def AWUDEB_sdf():
    return IO_DATA_PATH / "AWUDEB.sdf"


def HXACAN_sdf():
    return IO_DATA_PATH / "HXACAN.sdf"


def AWUDEB_cif():
    return IO_DATA_PATH / "AWUDEB.cif"


def HXACAN_cif():
    return IO_DATA_PATH / "HXACAN.cif"


def TRICKY_cif():
    # From https://www.crystallography.net/cod/2102215.html
    return IO_DATA_PATH / "2102215.cif"


@pytest.mark.skipif(sys.platform.startswith("darwin"), reason="no cif reader for osx")
@pytest.mark.parametrize(
    "file, should_load",
    [
        (HXACAN_cif(), True),
        (AWUDEB_cif(), True),
        (TRICKY_cif(), False),
        (HXACAN_sdf(), True),
        (AWUDEB_sdf(), True),
    ],
)
def test_cif_loading_defaults(file, should_load):
    if should_load:
        load_file(file)
    else:
        with pytest.raises(MoleculeParseError):
            load_file(file)


def test_pdb_loading():
    load_file("examples/large_molecules/CGMD/humera/20240627_humira_noglyco_ph6.pdb")
