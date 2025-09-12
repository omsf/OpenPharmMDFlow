"""
Settings for running the Small Molecule Pipeline
"""

from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from openff.interchange.components._packmol import RHOMBIC_DODECAHEDRON
from openff.interchange.components._packmol import UNIT_CUBE
from openff.models.types import FloatQuantity
from openff.models.types import Quantity
from openff.toolkit import ForceField
from pydantic.v1 import BaseModel


class SmallMoleculePipelineInputConfig(BaseModel):
    # TODO: support SMILES
    name: str
    path: Path


class BespokeWorkflowFactoryConfig(BaseModel):
    initial_force_field: str = "openff-2.2.0.offxml"
    qc_method: str = "gfn2xtb"
    qc_basis: str | None = None
    qc_program: str = "xtb"
    qc_spec_name: str = "xtb"
    qc_spec_description: str = "gfn2xtb"


class BespokeExecutorConfig(BaseModel):
    n_fragmenter_workers: int = 1
    n_optimizer_workers: int = 1
    n_qc_compute_workers: int = 2
    n_bespoke_workers: int = 1
    n_fragment_workers: int = 1


class BespokefitConfig(BaseModel):
    # Right now this works to optimize one small molecule
    # we could hack together something that works but it would
    # be best to run this out of band and just load in the
    # bespoke forcefeild
    bespoke_workflow_factory_config: BespokeWorkflowFactoryConfig
    bespoke_executor_config: BespokeExecutorConfig
    save_bespoke_ff: bool = True
    mol_to_bespoke: str


class SmallMoleculePipelinePrepConfig(BaseModel):
    bespokefit_config: BespokefitConfig


class SmallMoleculePipelinePackConfig(BaseModel):
    # TODO: add validator for len(molecule_names) == len(number_of_copies)
    # TODO: add validator for box_shape

    # for packing
    molecule_names: list[str]
    number_of_copies: list[int]
    target_density: FloatQuantity["g/cm**3"]
    box_shape: np.ndarray = UNIT_CUBE

    class Config:
        arbitrary_types_allowed = True


class SmallMoleculePipelineSolvateConfig(BaseModel):
    # for solvating
    # TODO fix this so the defaults work
    nacl_conc: Quantity = (Quantity(0.1, "mole / liter"),)
    padding: Quantity = (Quantity(1.2, "nanometer"),)
    box_shape: np.ndarray = (RHOMBIC_DODECAHEDRON,)
    target_density: Quantity = (Quantity(0.9, "gram / milliliter"),)
    tolerance: Quantity = (Quantity(2.0, "angstrom"),)

    class Config:
        arbitrary_types_allowed = True


class SmallMoleculePipelineParameterizeConfig(BaseModel):
    force_field: ForceField | str | Path | None = "openff-2.2.1.offxml"

    class Config:
        arbitrary_types_allowed = True


class SmallMoleculePipelineSimulateConfig(BaseModel):
    save_frequency_steps: int = 500
    save_data_frequency_steps: int = 10
    trajectory_name: str = "trajectory.pdb"
    temp_k: float = 300
    time_step_fs: int = 1
    pressure_bar: float = 1
    n_steps: int = 5000


class SmallMoleculePipelineAnalyizeConfig(BaseModel):
    pass


class SmallMoleculePipelineConfig(BaseModel):
    work_dir: Path
    inputs: list[SmallMoleculePipelineInputConfig] | SmallMoleculePipelineInputConfig
    prep_config: SmallMoleculePipelinePrepConfig | None
    pack_config: SmallMoleculePipelinePackConfig
    solvate_config: SmallMoleculePipelineSolvateConfig | None
    parameterize_config: SmallMoleculePipelineParameterizeConfig
    simulate_config: SmallMoleculePipelineSimulateConfig
    analyize_config: SmallMoleculePipelineAnalyizeConfig
