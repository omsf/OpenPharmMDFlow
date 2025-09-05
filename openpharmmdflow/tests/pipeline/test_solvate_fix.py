"""
Test for the solvate box_vectors/padding fix
"""

import pytest
import tempfile
from pathlib import Path

from openpharmmdflow.pipeline.sm.pipeline_settings import *
from openpharmmdflow.pipeline.sm.pipeline import SmallMoleculePipeline


def test_solvate_with_predefined_box_vectors():
    """Test that solvate method works when topology has box_vectors from pack step"""
    
    # Create a minimal SDF file for testing
    test_sdf_content = """
  Mrv2411 01010000002D          

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test molecule file
        mol_path = Path(tmpdir) / "test_mol.sdf"
        mol_path.write_text(test_sdf_content)
        
        # Set up pipeline configuration 
        settings = {
            "work_dir": tmpdir,
            "inputs": [SmallMoleculePipelineInputConfig(name="test_mol", path=mol_path)],
            "pack_config": SmallMoleculePipelinePackConfig(
                molecule_names=["test_mol"],
                number_of_copies=[1],
                target_density=0.1,
            ),
            "parameterize_config": SmallMoleculePipelineParameterizeConfig(),
            "simulate_config": SmallMoleculePipelineSimulateConfig(),
            "analyize_config": SmallMoleculePipelineAnalyizeConfig(),
            "solvate_config": SmallMoleculePipelineSolvateConfig(
                nacl_conc=Quantity(0.00, "mole / liter"),
                padding=Quantity(1.2, "nanometer"),  # This should get set to None
                box_shape=UNIT_CUBE,
                target_density=Quantity(0.9, "gram / milliliter"),
                tolerance=Quantity(2.0, "angstrom"),
            ),
        }
        
        sm_config = SmallMoleculePipelineConfig(**settings)
        
        # Initialize and run pipeline
        smp = SmallMoleculePipeline(sm_config)
        smp.load()
        smp.pack()
        
        # Verify that topology has box_vectors after packing
        assert smp.topology.box_vectors is not None, "Topology should have box_vectors after packing"
        
        # This should not raise PACKMOLValueError anymore
        smp.solvate()
        
        # Verify solvated topology was created
        assert hasattr(smp, 'solvated_topology'), "Should have solvated_topology after solvate()"
        assert smp.solvated_topology is not None, "solvated_topology should not be None"


def test_solvate_without_config_raises_error():
    """Test that solvate method raises error when solvate_config is None"""
    
    test_sdf_content = """
  Mrv2411 01010000002D          

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test molecule file
        mol_path = Path(tmpdir) / "test_mol.sdf"
        mol_path.write_text(test_sdf_content)
        
        # Set up pipeline configuration WITHOUT solvate_config
        settings = {
            "work_dir": tmpdir,
            "inputs": [SmallMoleculePipelineInputConfig(name="test_mol", path=mol_path)],
            "pack_config": SmallMoleculePipelinePackConfig(
                molecule_names=["test_mol"],
                number_of_copies=[1],
                target_density=0.1,
            ),
            "parameterize_config": SmallMoleculePipelineParameterizeConfig(),
            "simulate_config": SmallMoleculePipelineSimulateConfig(),
            "analyize_config": SmallMoleculePipelineAnalyizeConfig(),
            "solvate_config": None,  # Explicitly set to None
        }
        
        sm_config = SmallMoleculePipelineConfig(**settings)
        
        # Initialize and run pipeline
        smp = SmallMoleculePipeline(sm_config)
        smp.load()
        smp.pack()
        
        # This should raise ValueError
        with pytest.raises(ValueError, match="solvate_config must be provided"):
            smp.solvate()
