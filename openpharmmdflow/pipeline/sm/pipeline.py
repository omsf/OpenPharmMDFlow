"""
Small Molecule Pipeline
"""

from openff.interchange import Interchange
from openff.interchange.components._packmol import pack_box
from openff.toolkit import ForceField

from openpharmmdflow.bespokefit import build_bespoke_workflow_factory
from openpharmmdflow.bespokefit import run_bespokefit
from openpharmmdflow.io.load import load_file
from openpharmmdflow.pipeline.sm import SmallMoleculePipelineInputConfig
from openpharmmdflow.pipeline.sm import create_simulation
from openpharmmdflow.pipeline.sm import run_simulation

# TODO: Use snakemake to manage pipeline?
# TODO: Use decorators for DAG/deps?


class SmallMoleculePipeline:
    # TODO: Track the stage
    # TODO: serialize
    # TODO: factory model?
    # TODO: if prep config, store nested configs more flat?
    def __init__(self, config: SmallMoleculePipelineInputConfig):
        self.config = config
        self.inputs = (
            config.inputs if isinstance(config.inputs, list) else [config.inputs]
        )
        self.prep_config = config.prep_config if config.prep_config else None
        self.pack_config = config.pack_config
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
        self.topology = pack_box(
            molecules=[
                self.loaded_mols[mol_name]
                for mol_name in self.pack_config.molecule_names
            ],
            number_of_copies=self.pack_config.number_of_copies,
            mass_density=self.pack_config.mass_density,
            box_shape=self.pack_config.box_shape,
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
        self.interchange = Interchange.from_smirnoff(
            force_field=self.force_field, topology=self.topology
        )

    def simulate(self):
        self.simulation = create_simulation(self.simulate_config, self.interchange)
        run_simulation(self.simulate_config, self.simulation)

    def analyize(self):
        # run analysis here
        pass
