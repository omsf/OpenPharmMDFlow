"""
Module for cif file operations
"""
import os

import openff.toolkit


def from_cif(path: str | os.PathLike) -> openff.toolkit.Topology:
    topology = openff.toolkit.Topology()
    return topology
