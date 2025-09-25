"""
OpenMM simulation wrapper following OpenMM best practices with optimized trajectory formats
"""

import os
import time
from typing import Union

import numpy as np
import openmm
import openmm.app
import openmm.unit
from openff.interchange import Interchange


def create_simulation(
    simulate_config,
    interchange: Interchange,
) -> openmm.app.Simulation:
    """Create OpenMM simulation following standard patterns with optimized reporters"""

    # Extract configuration parameters
    save_frequency_steps = simulate_config.save_frequency_steps
    save_data_frequency_steps = simulate_config.save_data_frequency_steps
    simulation_chunk_step_size = getattr(
        simulate_config, "simulation_chunk_step_size", save_data_frequency_steps
    )
    temp = simulate_config.temp_k
    time_step_fs = simulate_config.time_step_fs
    pressure_bar = simulate_config.pressure_bar
    output_directory = simulate_config.output_directory

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Use LangevinMiddleIntegrator (more stable than LangevinIntegrator)
    integrator = openmm.LangevinMiddleIntegrator(
        temp * openmm.unit.kelvin,
        time_step_fs / openmm.unit.picosecond,  # Friction coefficient
        time_step_fs * openmm.unit.femtoseconds,
    )

    # Create simulation from interchange
    simulation = interchange.to_openmm_simulation(
        combine_nonbonded_forces=True,
        integrator=integrator,
    )

    # Add barostat for NPT ensemble
    ensemble = getattr(simulate_config, "ensemble", "npt")
    if ensemble.lower() == "npt":
        barostat = openmm.MonteCarloBarostat(
            pressure_bar * openmm.unit.bar,
            temp * openmm.unit.kelvin,
            25,  # Frequency of volume moves
        )
        simulation.system.addForce(barostat)

        # Required after adding forces
        simulation.context.reinitialize(preserveState=True)

    # Set up optimized trajectory reporters
    setup_trajectory_reporters(simulation, simulate_config, output_directory)

    return simulation


def setup_trajectory_reporters(simulation, simulate_config, output_directory):
    """Set up optimized trajectory reporters for large simulations"""
    save_frequency_steps = simulate_config.save_frequency_steps
    save_data_frequency_steps = simulate_config.save_data_frequency_steps
    simulation_chunk_step_size = getattr(
        simulate_config, "simulation_chunk_step_size", save_data_frequency_steps
    )

    print(f"📁 Setting up trajectory reporters in: {output_directory}")

    # 1. DCD Reporter (Primary trajectory - binary, compact)
    dcd_file = os.path.join(output_directory, "trajectory.dcd")
    dcd_reporter = openmm.app.DCDReporter(
        dcd_file, save_frequency_steps, enforcePeriodicBox=True
    )
    simulation.reporters.append(dcd_reporter)
    print(f"   🔸 DCD trajectory: {dcd_file} (every {save_frequency_steps} steps)")

    # 2. State Data Reporter (Thermodynamic data - CSV format)
    data_file = os.path.join(output_directory, "simulation_data.csv")
    state_reporter = openmm.app.StateDataReporter(
        data_file,
        save_data_frequency_steps,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,  # Performance monitoring
        separator=",",
    )
    simulation.reporters.append(state_reporter)
    print(
        f"   🔸 Thermodynamic data: {data_file} (every {save_data_frequency_steps} steps)"
    )

    # 3. Checkpoint Reporter (For resuming simulations)
    checkpoint_file = os.path.join(output_directory, "checkpoint.chk")
    checkpoint_multiplier = getattr(
        simulate_config, "checkpoint_frequency_multiplier", 10
    )
    checkpoint_frequency = save_frequency_steps * checkpoint_multiplier
    checkpoint_reporter = openmm.app.CheckpointReporter(
        checkpoint_file, checkpoint_frequency
    )
    simulation.reporters.append(checkpoint_reporter)
    print(f"   🔸 Checkpoint: {checkpoint_file} (every {checkpoint_frequency} steps)")

    # 4. HDF5 Reporter (Optional - if MDTraj available, best for analysis)
    try:
        from mdtraj.reporters import HDF5Reporter

        h5_file = os.path.join(output_directory, "trajectory.h5")

        # Check if HDF5 file exists (indicates resumption)
        if os.path.exists(h5_file):
            print(
                f"   ⚠️ HDF5 trajectory exists, skipping to avoid file lock: {h5_file}"
            )
        else:
            h5_reporter = HDF5Reporter(
                h5_file,
                save_frequency_steps,
                coordinates=True,
                time=True,
                cell=True,
                potentialEnergy=True,
                temperature=True,
            )
            simulation.reporters.append(h5_reporter)
            print(f"   ✅ HDF5 trajectory: {h5_file} (analysis-optimized)")
    except ImportError:
        print("   ⚠️ MDTraj not available - HDF5 trajectory disabled")
        print("   💡 Install with: micromamba install -c conda-forge mdtraj")

    # 5. PDB Reporter (Sparse frames for visualization only)
    pdb_file = os.path.join(output_directory, "trajectory_vis.pdb")

    # Use pdb_frequency_multiplier if available, otherwise default to 100
    pdb_multiplier = getattr(simulate_config, "pdb_frequency_multiplier", 100)
    pdb_frequency = save_frequency_steps * pdb_multiplier

    pdb_reporter = openmm.app.PDBReporter(
        pdb_file, pdb_frequency, enforcePeriodicBox=True
    )
    simulation.reporters.append(pdb_reporter)
    print(f"   🔸 PDB visualization: {pdb_file} (every {pdb_frequency} steps)")


def run_simulation(
    simulate_config,
    simulation: openmm.app.Simulation,
):
    """Run simulation following OpenMM best practices with improved monitoring"""
    n_steps = simulate_config.n_steps
    save_data_frequency_steps = simulate_config.save_data_frequency_steps
    simulation_chunk_step_size = getattr(
        simulate_config, "simulation_chunk_step_size", save_data_frequency_steps
    )
    temp_k = simulate_config.temp_k
    time_step_fs = simulate_config.time_step_fs

    print("🚀 Starting simulation setup...")
    start_time = time.time()

    # Energy minimization
    print("⚡ Minimizing energy...")
    try:
        simulation.minimizeEnergy(maxIterations=1000)
        minimized_state = simulation.context.getState(getEnergy=True)
        pe_initial = minimized_state.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoules_per_mole
        )
        print(f"✅ Energy minimization completed: PE = {pe_initial:.1f} kJ/mol")
    except Exception as e:
        print(f"⚠️ Energy minimization failed: {e}")
        print("   Continuing with original structure...")

    # Set initial velocities
    simulation.context.setVelocitiesToTemperature(temp_k * openmm.unit.kelvin)
    print(f"🌡️ Initial velocities set to {temp_k} K")

    # Compute virtual sites if present
    simulation.context.computeVirtualSites()

    print(f"🏃 Starting MD simulation: {n_steps:,} steps")
    print(f"   Timestep: {time_step_fs} fs")
    print(f"   Total time: {n_steps * time_step_fs / 1e6:.2f} ns")
    print()
    print("Step,      Time (ps), Volume (nm³), Temperature (K), PE (kJ/mol)")
    print("-" * 65)

    # Standard OpenMM simulation loop - run in chunks for monitoring
    chunk_size = min(
        simulation_chunk_step_size, 50000
    )  # Use new parameter for chunk size

    for start_step in range(0, n_steps, chunk_size):
        end_step = min(start_step + chunk_size, n_steps)
        steps_to_run = end_step - start_step

        try:
            # Run chunk of simulation
            simulation.step(steps_to_run)

            # Progress monitoring at configured intervals using new parameter
            if start_step % simulation_chunk_step_size == 0 or end_step == n_steps:
                state = simulation.context.getState(getEnergy=True)

                # Extract state information
                step_num = end_step
                time_ps = step_num * time_step_fs / 1000.0

                # Get box volume
                box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
                volume_nm3 = np.linalg.det(
                    box_vectors.value_in_unit(openmm.unit.nanometer)
                )

                # Get energies
                pe = state.getPotentialEnergy().value_in_unit(
                    openmm.unit.kilojoules_per_mole
                )

                # Calculate temperature from integrator (most reliable method)
                # This avoids unit conversion issues with kinetic energy
                temp_k_current = simulation.integrator.getTemperature().value_in_unit(
                    openmm.unit.kelvin
                )

                print(
                    f"{step_num:8d}, {time_ps:10.2f}, {volume_nm3:10.3f}, {temp_k_current:12.1f}, {pe:12.1f}"
                )

        except Exception as e:
            print(f"❌ Simulation error at step {start_step}: {e}")

            # Save emergency checkpoint
            output_dir = getattr(simulate_config, "output_directory", ".")
            emergency_checkpoint = os.path.join(
                output_dir, f"emergency_checkpoint_{start_step}.chk"
            )
            simulation.saveCheckpoint(emergency_checkpoint)
            print(f"💾 Emergency checkpoint saved: {emergency_checkpoint}")

            raise

    end_time = time.time()
    elapsed_time = end_time - start_time

    print("\n" + "=" * 65)
    print("🎉 Simulation completed successfully!")
    print(f"⏱️ Total time: {elapsed_time:.2f} seconds ({elapsed_time/3600:.2f} hours)")
    print(f"⚡ Performance: {n_steps/elapsed_time:.0f} steps/second")

    # Calculate simulation rate metrics
    simulated_time_ns = n_steps * time_step_fs / 1e6
    ns_per_day = simulated_time_ns * 86400 / elapsed_time
    print(f"📊 Sampling rate: {ns_per_day:.2f} ns/day")
    print(f"📈 Simulated time: {simulated_time_ns:.2f} ns")

    return elapsed_time, ns_per_day


def resume_simulation(
    simulate_config,
    simulation: openmm.app.Simulation,
    checkpoint_file: str | None = None,
):
    """Resume simulation from checkpoint"""
    output_dir = getattr(simulate_config, "output_directory", ".")

    if checkpoint_file is None:
        checkpoint_file = os.path.join(output_dir, "checkpoint.chk")

    if not os.path.exists(checkpoint_file):
        raise FileNotFoundError(f"Checkpoint file not found: {checkpoint_file}")

    print(f"🔄 Resuming simulation from: {checkpoint_file}")

    # Load checkpoint
    simulation.loadCheckpoint(checkpoint_file)

    # Get current state info
    state = simulation.context.getState(getEnergy=True)
    pe = state.getPotentialEnergy()
    print(f"   Resumed state: PE = {pe}")

    # Continue simulation
    run_simulation(simulate_config, simulation)
