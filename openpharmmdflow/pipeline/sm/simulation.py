"""
Wrapper for running openmm simulation
"""

import time

import numpy as np
import openmm
import openmm.app
import openmm.unit
from openff.interchange import Interchange


# TODO: Resume simulation
def create_simulation(
    interchange: Interchange,
    pdb_stride: int = 500,
    trajectory_name: str = "trajectory.pdb",
) -> openmm.app.Simulation:
    integrator = openmm.LangevinIntegrator(
        300 * openmm.unit.kelvin,
        1 / openmm.unit.picosecond,
        1 * openmm.unit.femtoseconds,
    )

    barostat = openmm.MonteCarloBarostat(
        1.0 * openmm.unit.bar, 293.15 * openmm.unit.kelvin, 25
    )

    simulation = interchange.to_openmm_simulation(
        combine_nonbonded_forces=True,
        integrator=integrator,
    )

    simulation.system.addForce(barostat)

    # https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#why-does-it-ignore-changes-i-make-to-a-system-or-force
    simulation.context.reinitialize(preserveState=True)

    # https://github.com/openmm/openmm/issues/3736#issuecomment-1217250635
    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(300 * openmm.unit.kelvin)
    simulation.context.computeVirtualSites()

    pdb_reporter = openmm.app.PDBReporter(trajectory_name, pdb_stride)
    state_data_reporter = openmm.app.StateDataReporter(
        "data.csv",
        10,
        step=True,
        potentialEnergy=True,
        temperature=True,
        density=True,
    )
    simulation.reporters.append(pdb_reporter)
    simulation.reporters.append(state_data_reporter)

    return simulation


def run_simulation(simulation: openmm.app.Simulation, n_steps: int = 5000):
    print("Starting simulation")
    start_time = time.process_time()

    print("Step, volume (nm^3)")

    for step in range(n_steps):
        simulation.step(1)
        if step % 500 == 0:
            box_vectors = simulation.context.getState().getPeriodicBoxVectors()
            print(step, np.linalg.det(box_vectors._value).round(3))

    end_time = time.process_time()
    print(f"Elapsed time: {(end_time - start_time):.2f} seconds")
