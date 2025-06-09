"""
Settings for running the Small Molecule Pipeline
"""

from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from openff.interchange.components._packmol import RHOMBIC_DODECAHEDRON
from openff.interchange.components._packmol import UNIT_CUBE
from openff.units import Quantity
from openff.toolkit import ForceField
from pydantic import BaseModel, model_validator

from typing import TYPE_CHECKING, Any
from openff.units import Unit, unit
from openff.utilities import has_package, requires_package
import openmm


class MissingUnitError(ValueError):
    """
    Exception for data missing a unit tag.
    """


class UnitValidationError(ValueError):
    """
    Exception for bad behavior when validating unit-tagged data.
    """


def _is_openmm_quantity(obj: object) -> bool:
    if has_package("openmm"):
        import openmm.unit

        return isinstance(obj, openmm.unit.Quantity)

    else:
        return "openmm.unit.quantity.Quantity" in str(type(object))


@requires_package("openmm.unit")
def _from_omm_quantity(val: "openmm.unit.Quantity") -> Quantity:
    """
    Convert float or array quantities tagged with SimTK/OpenMM units to a Pint-compatible quantity.
    """
    unit_: openmm.unit.Unit = val.unit
    val_ = val.value_in_unit(unit_)
    if type(val_) in {float, int}:
        unit_ = val.unit
        return float(val_) * Unit(str(unit_))
    # Here is where the toolkit's ValidatedList could go, if present in the environment
    elif (type(val_) in {tuple, list, np.ndarray}) or (type(val_).__module__ == "openmm.vec3"):
        array = np.asarray(val_)
        return array * Unit(str(unit_))
    elif isinstance(val_, (float, int)) and type(val_).__module__ == "numpy":
        return val_ * Unit(str(unit_))
    else:
        raise UnitValidationError(
            "Found a openmm.unit.Unit wrapped around something other than a float-like "
            f"or numpy.ndarray-like. Found a unit wrapped around type {type(val_)}."
        )


class _FloatQuantityMeta(type):
    def __getitem__(self, t):
        return type("FloatQuantity", (FloatQuantity,), {"__unit__": t})



class FloatQuantity(float, metaclass=_FloatQuantityMeta):
    """A model for unit-bearing floats."""

    @classmethod
    def __get_validators__(cls):
        yield cls.validate_type

    @model_validator(mode='before')
    @classmethod
    def validate_type(cls, val):
        """Process a value tagged with units into one tagged with "OpenFF" style units."""
        unit_ = getattr(cls, "__unit__", Any)
        if unit_ is Any:
            if isinstance(val, (float, int)):
                # TODO: Can this exception be raised with knowledge of the field it's in?
                raise MissingUnitError(f"Value {val} needs to be tagged with a unit")
            elif isinstance(val, Quantity):
                return Quantity(val)
            elif _is_openmm_quantity(val):
                return _from_omm_quantity(val)
            else:
                raise UnitValidationError(f"Could not validate data of type {type(val)}")
        else:
            unit_ = Unit(unit_)
            if isinstance(val, Quantity):
                # some custom behavior could go here
                assert unit_.dimensionality == val.dimensionality
                # return through converting to some intended default units (taken from the class)
                val._magnitude = float(val.m)
                return val.to(unit_)

            if _is_openmm_quantity(val):
                return _from_omm_quantity(val).to(unit_)
            if isinstance(val, int) and not isinstance(val, bool):
                # coerce ints into floats for a FloatQuantity
                return float(val) * unit_
            if isinstance(val, float):
                return val * unit_
            if isinstance(val, str):
                # could do custom deserialization here?
                val = Quantity(val).to(unit_)
                val._magnitude = float(val._magnitude)
                return val
            if "unyt" in str(val.__class__):
                if val.value.shape == ():
                    # this is a scalar force into an array by unyt's design
                    if "float" in str(val.value.dtype):
                        return float(val.value) * unit_
                    elif "int" in str(val.value.dtype):
                        return int(val.value) * unit_

            raise UnitValidationError(f"Could not validate data of type {type(val)}")

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
    mass_density: FloatQuantity["g/cm**3"]
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
    pdb_stride: int = 500
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
