"""
Module for cif file operations
"""

import os
import shutil
import subprocess
import tempfile

import openff.toolkit

from openpharmmdflow.apptainer.apptainer import GOODBOY


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
        If the required external tool 'apptainer' is not found.

    Notes
    -----
    This function utilizes the 'codcif2sdf' external tool inside an apptainer image
    to convert a CIF file to a Structure Data File (SDF), and then reads the SDF file
    into an OpenFF Toolkit Molecule object. Bond order is guessed.

    Examples
    --------
    >>> from openpharmmdflow import from_cif
    >>> molecule = from_cif("/path/to/example.cif")
    >>> molecule.visualize()
    """

    if shutil.which("apptainer"):
        CONTAINER_EXE_NAME = "apptainer"

    elif shutil.which("singularity"):
        CONTAINER_EXE_NAME = "singularity"

    else:
        raise RuntimeError("apptainer or singularity is required to read cif")

    # Download the apptainer image if it isn't in the cache
    container_path = GOODBOY.fetch("cod-tools.sif")

    with tempfile.TemporaryDirectory() as td:
        out_path = f"{td}/out.sdf"
        with open(f"{td}/out.sdf", "w") as f:
            cmd_output = subprocess.run(
                [CONTAINER_EXE_NAME, "exec", container_path, "codcif2sdf", f"{path}"],
                check=True,
                stdout=f,
            )
            molecule = openff.toolkit.Molecule.from_file(out_path)

    return molecule
