"""
Utility functions for Small Molecule Pipeline
"""

def write_residue_names(pipeline, molecule_to_resname=None):
    """
    Assign residue names to atoms in molecules for proper PDB output.
    
    Parameters
    ----------
    pipeline : SmallMoleculePipeline
        The pipeline object with loaded topology
    molecule_to_resname : dict, optional
        Mapping of molecule names to residue names. If None, uses first 3 letters uppercase.
        
    Example
    -------
    molecule_to_resname = {
        "ibuprofen": "IBU",
        "CAP": "CAP",
    }
    write_residue_names(smp, molecule_to_resname)
    """
    if molecule_to_resname is None:
        molecule_to_resname = {}
    
    # Print loaded molecule names to confirm
    for mol in pipeline.topology.molecules:
        print(f"Loaded molecule: {mol.name}")
    
    # Assign residue names to atoms
    for molecule in pipeline.topology.molecules:
        resname = molecule_to_resname.get(molecule.name, molecule.name[:3].upper())
        for atom in molecule.atoms:
            if not hasattr(atom, "metadata") or atom.metadata is None:
                atom.metadata = {}
            atom.metadata["residue_name"] = resname
    
    # Verify assignment worked
    molecule_list = list(pipeline.topology.molecules)
    for i, molecule in enumerate(molecule_list[:3]):
        print(f"Molecule {i}: {molecule.name} -> First atom residue_name = {molecule.atoms[0].metadata['residue_name']}")
