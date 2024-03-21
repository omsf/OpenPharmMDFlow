"""
Module for cif file operations
"""
import os
import shutil
import subprocess
import tempfile

import openff.toolkit


def from_cif(path: str | os.PathLike) -> openff.toolkit.Molecule:
    """
    Convert a Crystallographic Information File (CIF) to an OpenFF Toolkit Molecule.

    Parameters
    ----------
    path : str or os.PathLike
        The path to the CIF file.

    Returns
    -------
    openff.toolkit.Molecule
        An OpenFF Toolkit Molecule object representing the molecular structure.

    Raises
    ------
    RuntimeError
        If the required external tool 'codcif2sdf' is not found.

    Notes
    -----
    This function utilizes the 'codcif2sdf' external tool to convert a CIF file
    to a Structure Data File (SDF), and then reads the SDF file into an OpenFF
    Toolkit Molecule object. Bond order is guessed.

    Examples
    --------
    >>> from openpharmmdflow import from_cif
    >>> molecule = from_cif("/path/to/example.cif")
    >>> molecule.visualize()
    """

    if shutil.which("codcif2sdf") is None:
        raise RuntimeError("codcif2sdf is required to read cif")
    with tempfile.TemporaryDirectory() as td:
        out_path = f"{td}/out.sdf"
        with open(f"{td}/out.sdf", "w") as f:
            cmd_output = subprocess.run(
                ["codcif2sdf", f"{path}"],
                check=True,
                env={"PATH": "/usr/bin/"},
                stdout=f,
            )
            molecule = openff.toolkit.Molecule.from_file(out_path)

    return molecule
