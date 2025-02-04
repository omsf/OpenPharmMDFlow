"""
Wrapper for loading componets
"""

from pathlib import Path
from typing import Union

from openff.toolkit import Molecule
from openff.toolkit import Topology

from openpharmmdflow import from_cif


def _to_path(path: str | Path) -> Path:
    if isinstance(path, str):
        return Path(path)
    elif isinstance(path, Path):
        return path
    else:
        raise TypeError("Input must be a string or pathlib.Path object.")


def _path_exists(path: Path) -> bool:
    return path.exists()


def load_file(path: str | Path, **openff_kwargs) -> Molecule:
    path = _to_path(path)
    if not _path_exists(path):
        raise OSError(f"File at {path} does not exist.")

    if path.suffix == ".cif":
        return from_cif(path)

    if path.suffix == ".pdb":
        return Topology.from_pdb(path, **openff_kwargs)

    return Molecule.from_file(path, **openff_kwargs)
