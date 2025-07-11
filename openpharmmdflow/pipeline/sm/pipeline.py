"""
Small Molecule Pipeline
"""

import os

from openff.interchange import Interchange
from openff.interchange.components._packmol import pack_box
from openff.interchange.components._packmol import solvate_topology
from openff.toolkit import ForceField
from openff.toolkit import Molecule

from openpharmmdflow.io.load import load_file
from openpharmmdflow.pipeline.sm.pipeline_settings import SmallMoleculePipelineConfig
from openpharmmdflow.pipeline.sm.simulation import create_simulation
from openpharmmdflow.pipeline.sm.simulation import run_simulation

# TODO: Use snakemake to manage pipeline?
# TODO: Use decorators for DAG/deps?


class SmallMoleculePipeline:
    # TODO: Track the stage
    # TODO: serialize
    # TODO: factory model?
    # TODO: if prep config, store nested configs more flat?
    def __init__(self, config: SmallMoleculePipelineConfig):
        self.config = config
        self.inputs = (
            config.inputs if isinstance(config.inputs, list) else [config.inputs]
        )
        self.prep_config = config.prep_config if config.prep_config else None
        self.pack_config = config.pack_config
        self.solvate_config = config.solvate_config if config.solvate_config else None
        self.parameterize_config = config.parameterize_config
        self.bespoke_ff = None
        self.simulate_config = config.simulate_config

    def load(self):
        # TODO raise error if duplcate name
        self.loaded_mols = {}
        for input in self.inputs:
            mol = load_file(input.path)
            mol.name = input.name
            self.loaded_mols[mol.name] = mol

    def prep(self):
        # run bespokefit here
        from openpharmmdflow.bespokefit import build_bespoke_workflow_factory
        from openpharmmdflow.bespokefit import run_bespokefit

        if self.prep_config:
            self.factory = build_bespoke_workflow_factory(
                self.prep_config.bespokefit_config.bespoke_workflow_factory_config
            )
            self.bespoke_ff = run_bespokefit(
                self.prep_config.bespokefit_config,
                self.loaded_mols[self.prep_config.bespokefit_config.mol_to_bespoke],
                self.factory,
            )
        else:
            print("Nothing to prep")

    def pack(self):
        # build the box here
        # TODO: use mBuild for lattice tooling
        # Right now random packing is supported
        try:
            self.topology = pack_box(
                molecules=[
                    self.loaded_mols[mol_name]
                    for mol_name in self.pack_config.molecule_names
                ],
                number_of_copies=self.pack_config.number_of_copies,
                target_density=self.pack_config.target_density,
                box_shape=self.pack_config.box_shape,
            )
        # older versions of interchange use "mass_density" instead of target_density
        except TypeError:
            self.topology = pack_box(
                molecules=[
                    self.loaded_mols[mol_name]
                    for mol_name in self.pack_config.molecule_names
                ],
                number_of_copies=self.pack_config.number_of_copies,
                mass_density=self.pack_config.target_density,
                box_shape=self.pack_config.box_shape,
            )

    def solvate(self):
        # Solvate the box
        self.solvated_topology = solvate_topology(
            self.topology,
            nacl_conc=self.solvate_config.nacl_conc,
            padding=self.solvate_config.padding,
            box_shape=self.solvate_config.box_shape,
            target_density=self.solvate_config.target_density,
            tolerance=self.solvate_config.tolerance,
        )
        # create a water molecule
        self.water = Molecule.from_smiles("O")
        self.water.generate_conformers(n_conformers=1)
        # TODO See if we can get the number of waters, Na, Cl out of solvate_topology
        # TODO Maybe import water + ions instead of making them here?
        # find out the number of waters in the system
        self.n_water = len(
            [
                m
                for m in self.solvated_topology.molecules
                if m.to_smiles() == "[H][O][H]"
            ]
        )
        # create a Sodium ion
        self.sodium_ion = Molecule.from_smiles("[Na+]")
        # find out the number of sodium ions in the system
        self.n_sodium_ion = len(
            [m for m in self.solvated_topology.molecules if m.to_smiles() == "[Na+]"]
        )
        # create a Chlorine ion
        self.chlorine_ion = Molecule.from_smiles("[Cl-]")
        # find out the number of sodium ions in the system
        self.n_chlorine_ion = len(
            [m for m in self.solvated_topology.molecules if m.to_smiles() == "[Cl-]"]
        )

    def parameterize(self):
        # TODO test to make sure we use the FF we expect to use
        # Use bespoke ff if we made one
        self.force_field = (
            self.bespoke_ff if self.bespoke_ff else self.parameterize_config.force_field
        )
        # Now if force_field is a path or string, we need to turn it into a ForceField object
        if not isinstance(self.force_field, ForceField):
            self.force_field = ForceField(self.force_field)
        self.components_intrcg = Interchange.from_smirnoff(
            force_field=self.force_field, topology=self.topology
        )
        # if there are waters built during the solvate step combine the components topology
        # with the water topology
        # TODO better way to check if we ran solvate?
        # TODO let user specify ff for solute
        if hasattr(self, "water"):
            self.water_intrcg = Interchange.from_smirnoff(
                force_field=ForceField("openff_unconstrained-2.0.0.offxml"),
                topology=[self.water] * self.n_water
                + [self.sodium_ion] * self.n_sodium_ion
                + [self.chlorine_ion] * self.n_chlorine_ion,
            )
            # combine is still experimental
            os.environ["INTERCHANGE_EXPERIMENTAL"] = "1"
            self.interchange = self.components_intrcg.combine(self.water_intrcg)
            self.interchange.positions = self.solvated_topology.get_positions()
            self.interchange.box = self.solvated_topology.box_vectors

        else:
            self.interchange = self.components_intrcg

    def simulate(self):
        self.simulation = create_simulation(self.simulate_config, self.interchange)
        run_simulation(self.simulate_config, self.simulation)

    def analyize(self):
        # run analysis here
        pass
