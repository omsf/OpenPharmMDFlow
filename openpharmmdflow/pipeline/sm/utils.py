"""
Utility functions for Small Molecule Pipeline
"""

def write_residue_names(pipeline, molecule_to_resname=None):
    """
    Assign residue names to atoms in molecules for proper PDB output.
    
    Parameters
    ----------
    pipeline : SmallMoleculePipeline
        The pipeline object with loaded topology (works with both topology and solvated_topology)
    molecule_to_resname : dict, optional
        Mapping of molecule names to residue names. If None, uses first 3 letters uppercase.
        Standard names: water -> "WAT", sodium -> "SOD", chloride -> "CLA"
        
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
    
    # Use solvated_topology if available, otherwise topology
    topology = getattr(pipeline, 'solvated_topology', pipeline.topology)
    
    # Print loaded molecule names to confirm
    print("Assigning residue names to molecules:")
    
    # Assign residue names to atoms
    for molecule in topology.molecules:
        # Get SMILES to identify molecule type
        smiles = molecule.to_smiles()
        
        # Default residue names for common molecules
        default_names = {
            "[H][O][H]": "WAT",  # Water
            "[Na+]": "SOD",      # Sodium ion
            "[Cl-]": "CLA",      # Chloride ion
        }
        
        # Determine residue name
        if molecule.name in molecule_to_resname:
            resname = molecule_to_resname[molecule.name]
        elif smiles in default_names:
            resname = default_names[smiles]
        else:
            resname = molecule.name[:3].upper() if molecule.name else "UNK"
        
        # Assign to all atoms in the molecule
        for atom in molecule.atoms:
            if not hasattr(atom, "metadata") or atom.metadata is None:
                atom.metadata = {}
            atom.metadata["residue_name"] = resname
        
        print(f"  {molecule.name} (SMILES: {smiles}) -> {resname}")
    
    # Verify assignment worked
    molecule_list = list(topology.molecules)
    print(f"\nVerification (showing first 3 molecules):")
    for i, molecule in enumerate(molecule_list[:3]):
        first_atom_resname = molecule.atoms[0].metadata.get('residue_name', 'NONE')
        print(f"  Molecule {i}: {molecule.name} -> residue_name = {first_atom_resname}")
