print("importing packages")
from sys import stdout
from pathlib import Path
import os
import json
from openmm import *
from openmm.app import *
from openmm.unit import *
print("imports done")

# Simulation parameters (can be modified between runs)
dcd_write_period = 1000
log_write_period = 1000
nvt_steps = 10000
npt_steps = 10000

# Use these for faster testing if needed
dcd_write_period = 100
log_write_period = 100
nvt_steps = 1000
npt_steps = 1000

# Directories and file paths
out_dir = Path("A40")
dcd_path = out_dir / Path("output.dcd")
md_log_path = out_dir / Path("md_log.txt")
checkpoint_path = out_dir / Path("checkpoint.chk")
nvt_checkpoint_path = out_dir / Path("nvt_checkpoint.chk")
npt_checkpoint_path = out_dir / Path("npt_checkpoint.chk")
nvt_final_pdb_path = out_dir / Path("nvt_final.pdb")
npt_final_pdb_path = out_dir / Path("npt_final.pdb")
progress_file = out_dir / Path("progress.json")

# Create output directory
print("Setting up dirs")
out_dir.mkdir(parents=True, exist_ok=True)

# Check if we're resuming from a previous run
resuming = False
completed_nvt_steps = 0
completed_npt_steps = 0


class CheckpointReporter(object):
    def __init__(self, file, reportInterval, progress_file, phase):
        self._file = file
        self._reportInterval = reportInterval
        self._progress_file = progress_file
        self._phase = phase

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, False)

    def report(self, simulation, state):
        simulation.saveCheckpoint(self._file)
        completed_steps = completed_nvt_steps + simulation.currentStep if self._phase == 'nvt' else nvt_steps
        completed_npt_steps_now = completed_npt_steps + simulation.currentStep if self._phase == 'npt' else completed_npt_steps

        with open(self._progress_file, 'w') as f:
            json.dump({
                'completed_nvt_steps': completed_steps,
                'completed_npt_steps': completed_npt_steps_now,
                'minimization_completed': True
            }, f)




if progress_file.exists():
    print("Found progress file, checking simulation status")
    with open(progress_file, 'r') as f:
        progress_data = json.load(f)
        completed_nvt_steps = progress_data.get('completed_nvt_steps', 0)
        completed_npt_steps = progress_data.get('completed_npt_steps', 0)

    if completed_nvt_steps > 0 or completed_npt_steps > 0:
        resuming = True
        print(f"Resuming from previous run: NVT steps = {completed_nvt_steps}, NPT steps = {completed_npt_steps}")

# Calculate remaining steps
remaining_nvt_steps = max(0, nvt_steps - completed_nvt_steps)
remaining_npt_steps = max(0, npt_steps - completed_npt_steps)
total_steps = remaining_nvt_steps + remaining_npt_steps

print(f"Planning to run {remaining_nvt_steps} additional NVT steps and {remaining_npt_steps} additional NPT steps")

if not resuming:
    print("Loading pdb")
    pdb = PDBFile("../../inputs/mAb/lai_2022_mab3_nogly_ph6.pdb")
    print("Loading forcefield")
    forcefield = ForceField("amber14/protein.ff14SB.xml", "amber14/tip3pfb.xml")
    print("Setting up modeller")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, padding=1.0 * nanometer, boxShape="dodecahedron")
    print("Creating system")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
    )
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)
    platform = Platform.getPlatformByName('CUDA')
    simulation = Simulation(modeller.topology, system, integrator, platform=platform)
    simulation.context.setPositions(modeller.positions)

    print("Minimizing energy")
    simulation.minimizeEnergy()

    # Save checkpoint after minimization
    print("Saving checkpoint after minimization")
    simulation.saveCheckpoint(str(checkpoint_path))

    # Save minimized structure as PDB
    print("Saving minimized structure as PDB")
    minimized_pdb_path = out_dir / Path("minimized.pdb")
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(str(minimized_pdb_path), 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)

    # Update progress file
    with open(progress_file, 'w') as f:
        json.dump({
            'completed_nvt_steps': 0,
            'completed_npt_steps': 0,
            'minimization_completed': True
        }, f)

else:
    # Resume from checkpoint
    print("Loading forcefield to recreate system")
    forcefield = ForceField("amber14/protein.ff14SB.xml", "amber14/tip3pfb.xml")

    # We need to recreate the system, but we'll load coordinates from checkpoint
    print("Loading pdb to get topology")
    pdb = PDBFile("../../inputs/mAb/lai_2022_mab3_nogly_ph6.pdb")
    print("Setting up modeller")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, padding=1.0 * nanometer, boxShape="dodecahedron")

    # For NVT phase
    if completed_nvt_steps < nvt_steps:
        print("Resuming NVT phase")
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0 * nanometer,
            constraints=HBonds,
        )
        integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)
        platform = Platform.getPlatformByName('CUDA')
        simulation = Simulation(modeller.topology, system, integrator, platform=platform)

        # Load from the appropriate checkpoint
        if completed_nvt_steps == 0:
            print("Loading from post-minimization checkpoint")
            simulation.loadCheckpoint(str(checkpoint_path))
        else:
            print(f"Loading from NVT checkpoint with {completed_nvt_steps} steps completed")
            simulation.loadCheckpoint(str(nvt_checkpoint_path))

    # For NPT phase
    elif completed_npt_steps < npt_steps:
        print("Resuming NPT phase")
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0 * nanometer,
            constraints=HBonds,
        )
        system.addForce(MonteCarloBarostat(1 * bar, 300 * kelvin))
        integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)
        platform = Platform.getPlatformByName('CUDA')
        simulation = Simulation(modeller.topology, system, integrator, platform=platform)

        print(f"Loading from NPT checkpoint with {completed_npt_steps} steps completed")
        simulation.loadCheckpoint(str(npt_checkpoint_path))

# Set up reporters
simulation.reporters.append(DCDReporter(dcd_path, dcd_write_period, append=resuming))
simulation.reporters.append(
    StateDataReporter(
        md_log_path.as_posix(),
        log_write_period,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        totalSteps=total_steps,
        speed=True,
        remainingTime=True,
        elapsedTime=True,
        progress=True,
        append=resuming
    )
)

# Run NVT if needed
if remaining_nvt_steps > 0:
    print(f"Running NVT for {remaining_nvt_steps} steps")

    # Add a checkpoint reporter for regular checkpoints during long NVT runs
    if remaining_nvt_steps > 10000:
        checkpoint_frequency = 10000
    else:
        checkpoint_frequency = min(1000, remaining_nvt_steps)


    simulation.reporters.append(CheckpointReporter(
        str(nvt_checkpoint_path),
        checkpoint_frequency,
        str(progress_file),
        'nvt'
    ))

    # Run the simulation
    simulation.step(remaining_nvt_steps)

    # Save the final NVT structure as PDB
    print("Saving final NVT structure as PDB")
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(str(nvt_final_pdb_path), 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)

    # Update progress file after NVT completion
    with open(progress_file, 'w') as f:
        json.dump({
            'completed_nvt_steps': nvt_steps,
            'completed_npt_steps': completed_npt_steps,
            'minimization_completed': True
        }, f)

    # Save final NVT state as checkpoint
    simulation.saveCheckpoint(str(nvt_checkpoint_path))

    # Now we need to add barostat for NPT if NPT steps remain
    if remaining_npt_steps > 0:
        system.addForce(MonteCarloBarostat(1 * bar, 300 * kelvin))
        simulation.context.reinitialize(preserveState=True)

# Run NPT if needed
if remaining_npt_steps > 0:
    print(f"Running NPT for {remaining_npt_steps} steps")

    # Add checkpoint reporter for NPT phase
    if remaining_npt_steps > 10000:
        checkpoint_frequency = 10000
    else:
        checkpoint_frequency = min(1000, remaining_npt_steps)

    simulation.reporters.append(CheckpointReporter(
        str(npt_checkpoint_path),
        checkpoint_frequency,
        str(progress_file),
        'npt'
    ))

    # Run the simulation
    simulation.step(remaining_npt_steps)

    # Save the final NPT structure as PDB
    print("Saving final NPT structure as PDB")
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(str(npt_final_pdb_path), 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)

    # Update progress file after NPT completion
    with open(progress_file, 'w') as f:
        json.dump({
            'completed_nvt_steps': nvt_steps,
            'completed_npt_steps': npt_steps,
            'minimization_completed': True
        }, f)

    # Save final state
    simulation.saveCheckpoint(str(npt_checkpoint_path))

print("Simulation completed successfully!")
