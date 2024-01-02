#!/usr/bin/env python
import os

from pymatgen.io.cif import CifParser
from pymatgen.io.cif import CifWriter

cifdir = os.path.join(os.getcwd(), "cifs")
symdir = os.path.join(os.getcwd(), "symmed")
os.makedirs(symdir, exist_ok=True)

# Tunable for symm tolerance -- High tol. for MM calculations
symprec = 1e-1

for file in os.listdir(cifdir):
    if ".cif" in file:
        ciffile = os.path.join(cifdir, file)
        structure = CifParser(ciffile).get_structures()[0]
        writer = CifWriter(structure, symprec=1e-1)
        writer.write_file(os.path.join(symdir, file))
