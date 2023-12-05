#!/usr/bin/env python
# coding: utf-8
# atoms for large Z!!


import spglib as spg
import ase
import os
from ase.io import read,write
import numpy as np

from ase import neighborlist
from scipy.optimize import minimize

from openmm.app import *
from openmm import *
from openmm.unit import *
from pyxtal import pyxtal

from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit.topology import Molecule,Topology
from openmmforcefields.generators import SMIRNOFFTemplateGenerator,GAFFTemplateGenerator, EspalomaTemplateGenerator
from pymatgen.core.structure import IMolecule, IStructure, Structure#, Molecule
from pymatgen.io import cif

import glob
import pandas as pd

from pyxtal.io import read_cif,search_molecules_in_crystal
import subprocess 

atom_dict = {
    'P' : [15],
    'C' : [6],
    'H' : [1],
    'O' : [8],
    'N' : [7],
    'S' : [16],
    'Cl' : [17],
    'F' : [9]
}


def prep_molecule(workingdir,molecule,count):
    xyzfile = os.path.join(workingdir,'molec-{}.xyz'.format(count))
    molfile = os.path.join(workingdir,'molec-{}.mol'.format(count))
    with open(xyzfile,'w') as f:
        f.write('{}\n\n'.format(len(molecule)))
        for a in molecule:
            ID = a.as_dict()['name']
            xyz = a.as_dict()['xyz']
            f.write('{:s}\t{:f}\t{:f}\t{:f}\n'.format(ID,xyz[0],xyz[1],xyz[2]))
    Command = 'obabel ' + xyzfile + ' -O ' + molfile
    Command = Command.split()
    subprocess.call(Command, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    molecule = IMolecule.from_file(molfile)
    return molecule

def reorder(molec,template):
    #new = os.path.join('critic2',str(molec) +'-reordered.xyz')
    new = str(molec) +'-reordered.xyz'
    critic_inp = os.path.join('critic2','critic.cri')
    with open(critic_inp,'w') as f:
        f.write('MOLREORDER {:s} {:d}.xyz WRITE {:s}'.format(template,molec,new))
    Command = 'critic2 {}'.format(critic_inp)
    Command = Command.split()
    subprocess.call(Command, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    molfile = new.replace('.xyz','mol')
    Command = 'obabel ' + new + ' -O ' + molfile
    Command = Command.split()
    subprocess.call(Command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return read(os.path.join('critic2',new))

def clean_molecule(m):
    m = Chem.AddHs(m)
    for atom in m.GetAtoms():
        if atom.GetNumRadicalElectrons() > 0:
            print('FOUND FREE RADICAL, ATTEMPTING FIX...')
            atom.SetNumRadicalElectrons(0)
            if atom.GetSymbol() == 'O':
                atom.SetFormalCharge(-1)
            for bond in atom.GetBonds():
                if bond.GetBeginAtom().GetSymbol() == 'O' or bond.GetBeginAtom().GetSymbol() == 'N':
                    bond.GetBeginAtom().SetFormalCharge(1)
                    bond.SetBondType(Chem.BondType.DOUBLE)
    AllChem.EmbedMolecule(m)
    return(m)

def residue_setup(tagged_residues,images,ss_dim,file=None,forcefield='GAFF'):
    #print(len(tagged_residues) * images)
    parm_residue = [molfile for molfile in os.listdir(parmtarget) if molfile.endswith('.mol')]
    rd_mol = []
    for mol in parm_residue:
        m = Chem.MolFromMolFile(os.path.join(parmtarget,mol),removeHs=False)
        m = clean_molecule(m)
        molecule = Molecule.from_rdkit(m) 
        molecule.assign_partial_charges("gasteiger")
        #molecule.assign_partial_charges(partial_charge_method='am1bcc')
        rd_mol.append(molecule)
            
    if forcefield == 'GAFF':
        generator = GAFFTemplateGenerator()
    elif forcefield == 'smirnoff':
        generator = SMIRNOFFTemplateGenerator()
    elif forcefield == 'custom':
        # print('CUSTOM GEN', file)
        generator = SMIRNOFFTemplateGenerator(forcefield=file)
        generator.generate_residue_template(molecule)
    
    topology = Topology()
    count = 0
    for s in range(ss_dim):
        for i in range(images):
            for a in tagged_residues:
                molecule = rd_mol[a[1]]
                topology.add_molecule(molecule)
                generator.add_molecules(molecule)
                count += 1
    # print("N MOLECS: ", count)
    return topology, generator

def structure_setup(topology,generator,structure,**kwargs):
    topology.box_vectors = structure.lattice.matrix * angstrom
    topology.is_periodic = True


    omm_top = topology.to_openmm()
    # ff = ForceField()

    #Check if there are custom FF parameters
    for arg,val in kwargs.items():
        # print(arg,val)
        if arg == 'custom':
            iscustom = val
        if arg == 'ff_list':
            ff_list = val

    if not iscustom:
        ff = ForceField()
        ff.registerTemplateGenerator(generator.generator)
    else:
        ff = ForceField((',').join(ff_list))
        ff.registerTemplateGenerator(generator.generator)


    
    system = ff.createSystem(omm_top)
    context_coords = [structure.sites[i].coords/10 for i in range(0, len(structure.sites))]
    context_lattice = structure.lattice.matrix / 10
    system.setDefaultPeriodicBoxVectors(context_lattice[0],context_lattice[1],context_lattice[2])
    
    return system,context_coords,context_lattice,omm_top

def center_molecules(molecs,lattice):
    # lattice = structure.lattice.matrix
    for m in molecs:
        abc_com = np.matmul(np.linalg.inv(lattice).T, m.center_of_mass)
        center = np.zeros(3)
        for a in range(len(abc_com)):
            while abc_com[a] < 0:
                abc_com[a] += 1
                center[a] += 1
            while abc_com[a] > 1:
                abc_com[a] -= 1
                center[a] -= 1
        for xyz in m:
            abc = np.matmul(np.linalg.inv(lattice).T, xyz.coords)
            abc += center
            xyz.coords = np.matmul(lattice.T, abc)

    return molecs

class ProcessCell:
    def __init__(self,s):
        self.structure = s

    def standard_cell(self):
        numbers = [self.structure.sites[x].specie.number for x in range(len(self.structure.sites))]
        lattice = self.structure.lattice.matrix
        positions = self.structure.frac_coords
        cell = (lattice, positions, numbers)
        self.dataset = spg.get_symmetry_dataset(cell, symprec=1E-3)
        self.standard_structure = Structure(self.dataset['std_lattice'], self.dataset['std_types'], self.dataset['std_positions'])

    def make_asym(self):
        asym_ids = set(self.dataset['equivalent_atoms'])
        asym_coords = [self.dataset['std_positions'][i] for i in asym_ids]
        asym_atoms = [self.dataset['std_types'][i] for i in asym_ids]
        asym_lattice = self.dataset['std_lattice']
        self.asym = Structure(asym_lattice, asym_atoms, asym_coords)

        #Maybe remove supercell stuff...
        self.supercell_dim(lmin=20)

        #return self.standard_structure.make_supercell(self.ss)

    def ordered_crystal(self):
        test = {'count': [], 'parmtag': []}
        self.operations = [(r, t) for r, t in zip(self.dataset['rotations'], self.dataset['translations'])]

        #Search all chemically unique molecules in standardized cell for parameter assignment and ordering
        #Asym cell here???
        molecs = search_molecules_in_crystal(self.standard_structure)
        molecs = center_molecules(molecs, self.standard_structure.lattice.matrix)
        unique_num_atoms = set([len(m) for m in molecs]) #Chemically unique
        unique_molecules = []
        count = 0
        for a in unique_num_atoms:
            parmmolec = molecs[[len(m) for m in molecs].index(a)]
            molecule = prep_molecule(parmtarget, parmmolec, count)
            unique_molecules.append(molecule)
            test['parmtag'].append(count)
            test['count'].append(a)
            count += 1

        #Reorder all molecules according to their assigned parameters
        asym_molecules = search_molecules_in_crystal(self.asym)
        asym_molecules = center_molecules(asym_molecules,self.standard_structure.lattice.matrix)
        xtal_abc = []
        xtal_species = []
        ID = 0
        reordered_residue = []
        for m in asym_molecules:
            parmtag = test['count'].index(len(m))
            m.to(os.path.join('critic2', '{}.xyz'.format(ID)))
            reordered_molecule = reorder(ID, os.path.join(parmtarget, 'molec-{}.xyz'.format(parmtag)))
            tagged = [reordered_molecule, parmtag]
            reordered_residue.append(tagged)
            ID += 1

        e_a = 1 / (self.ss[0])
        e_b = 1 / (self.ss[1])
        e_c = 1 / (self.ss[2])

        e_T = np.zeros((3,3))
        e_T[0,0] = e_a
        e_T[1,1] = e_b
        e_T[2,2] = e_c
        for molec in reordered_residue:
            coords = molec[0].get_positions()
            species = molec[0].get_chemical_symbols()
            for t in self.translators:
                for o in self.operations:
                    for c in range(len(species)):
                        xyz = coords[c]
                        ID = species[c]
                        xtal_species.append(ID)
                        abc = (np.matmul(np.linalg.inv(self.standard_structure.lattice.matrix).T, xyz))
                        abc = np.matmul(o[0], abc) + o[1]  + t
                        abc = np.matmul(e_T.T,abc)
                        xtal_abc.append(abc)
        self.standard_structure.make_supercell(self.ss)
        ordered_cell = Structure(self.standard_structure.lattice.matrix, xtal_species, xtal_abc)
        return ordered_cell,reordered_residue

    def supercell_dim(self,lmin):
        abc = self.standard_structure.lattice.abc
        self.ss = [1,1,1]
        for v in abc:
            index = abc.index(v)
            i = 1
            N = v
            while N < lmin:
                i += 1
                N = v * i
            self.ss[index] = i

        self.translators = []
        for i in range(self.ss[0]):
            for j in range(self.ss[1]):
                for k in range(self.ss[2]):
                    self.translators.append([i,j,k])
                    k+=1
                j+=1
            i+=1




#parm_target
parmtarget = os.path.join(os.getcwd(),'parmmols')
mintarget = os.path.join(os.getcwd(),'minimized')
reordertarget = os.path.join(os.getcwd(),'critic2')
os.makedirs(parmtarget,exist_ok=True)
os.makedirs(mintarget,exist_ok=True)
os.makedirs(reordertarget, exist_ok=True)

#Troubleshooting
df = pd.DataFrame(columns=['structure','E(kJ/mol)','Z'])

if __name__ == '__main__':

    #Determine the ewald cutoffs by looking at all crystals:
    for f in os.listdir('cif'):
        print(f)
        test = {'count': [], 'parmtag': []}
        file = os.path.join('cif',f)
        xtalID = file.replace('cif/', '')

        structure = Structure.from_file(file)
        processed = ProcessCell(structure)
        processed.standard_cell()
        processed.make_asym()
        ordered,reordered = processed.ordered_crystal()
        ordered.to(
            'ordered.cif'
        )
        ss_dim = processed.ss[0]*processed.ss[1]*processed.ss[2]
        images = len(processed.operations)
        # Using custom parameters
        # topology, generator = residue_setup(reordered, images, ss_dim, file='FFs/ABX.offxml', forcefield='custom')
        # [system, context_coords, context_lattice, omm_top] = structure_setup(topology, generator, ordered,
        #                                                                     custom=True,ff_list=['FFs/ABX.offxml'])

        #Default SMIRNOFF
        topology, generator = residue_setup(reordered, images, ss_dim, forcefield='smirnoff')
        [system, context_coords, context_lattice, omm_top] = structure_setup(topology, generator, ordered,
                                                                            custom=False)

        Z = len([i for i in omm_top.chains()])
        #minimize E @ constant cell volume. Re-introduce NPT simulation for cell relaxations?
        nonb = system.getForces()[0]
        nonb.setNonbondedMethod(method=4) #This needs attention
        nonb.setCutoffDistance(0.9*nanometer)

        system.getForces()[0] = nonb

        integrator = VerletIntegrator(0.)
        integrator.setConstraintTolerance(1e-8)

        simulation = Simulation(omm_top,system,integrator)
        simulation.context.setPeriodicBoxVectors(context_lattice[0], context_lattice[1], context_lattice[2])
        simulation.context.setPositions(context_coords)

        s1 = simulation.context.getState(getPositions=True,getEnergy=True)
        startcoords = s1.getPositions().value_in_unit(angstrom)
        simulation.minimizeEnergy()
        s2 = simulation.context.getState(getPositions=True,getEnergy=True)
        mincoords = s2.getPositions().value_in_unit(angstrom)
        E_kJ = s2.getPotentialEnergy().value_in_unit(kilojoule_per_mole)

        df.loc[len(df.index)] = [f,E_kJ,Z]

        lattice = np.zeros((3, 3))

        lattice[0, 0] = system.getDefaultPeriodicBoxVectors()[0].x * 10
        lattice[0, 1] = system.getDefaultPeriodicBoxVectors()[0].y * 10
        lattice[0, 2] = system.getDefaultPeriodicBoxVectors()[0].z * 10

        lattice[1, 0] = system.getDefaultPeriodicBoxVectors()[1].x * 10
        lattice[1, 1] = system.getDefaultPeriodicBoxVectors()[1].y * 10
        lattice[1, 2] = system.getDefaultPeriodicBoxVectors()[1].z * 10

        lattice[2, 0] = system.getDefaultPeriodicBoxVectors()[2].x * 10
        lattice[2, 1] = system.getDefaultPeriodicBoxVectors()[2].y * 10
        lattice[2, 2] = system.getDefaultPeriodicBoxVectors()[2].z * 10

        abccoords = []
        for p in s2.getPositions():
            xyz = np.zeros(3)
            xyz[0] = p.x * 10
            xyz[1] = p.y * 10
            xyz[2] = p.z * 10
            abccoords.append(np.matmul(np.linalg.inv(lattice).T, xyz))

        atoms = []
        for a in omm_top.atoms():
            atoms.append(a.name[0])

        Structure(lattice, atoms, abccoords).to(os.path.join(mintarget,xtalID))

        # #remove all parameter files
        filelist = glob.glob(os.path.join(reordertarget,"*"))
        for f in filelist:
            os.remove(f)
        filelist = glob.glob(os.path.join(parmtarget, "*"))
        for f in filelist:
            os.remove(f)

    df.to_csv('MM_data.csv')


