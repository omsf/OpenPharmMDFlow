import os
import subprocess
import sys
from sys import stdout

import martini_openmm as martini
from mdtraj.reporters import XTCReporter
from openmm import *
from openmm.app import *
from openmm.unit import *


# Parameters
mab = "minimization-vac.gro"
water = "../water.gro"
default_nmab = 1
nmab = int(sys.argv[1]) if len(sys.argv) > 1 else default_nmab
default_box_size = 20
box_size = int(sys.argv[2]) if len(sys.argv) > 2 else default_box_size


# run_martinize.sh
command = [
    "martinize2",
    "-f",
    "antibody.pdb",
    "-o",
    "system.top",
    "-x",
    "antibody-CG.pdb",
    "-merge",
    "A,B,C,D",
    "-cys",
    "auto",
    "-dssp",
    "mkdssp",
]

subprocess.run(command, check=True)

# Write system.top file
system_top_content = f"""\
#include "../martini_v3.0.0.itp"

#include "molecule_0.itp"

[ system ]
; name
Martini system from antibody.pdb

[ molecules ]
; name        number
molecule_0    1
"""
with open("system.top", "w") as f:
    f.write(system_top_content)

# Run gmx commands
subprocess.run(
    [
        "gmx",
        "editconf",
        "-f",
        "antibody-CG.pdb",
        "-d",
        "2.0",
        "-bt",
        "cubic",
        "-o",
        "antibody-CG.gro",
    ]
)
subprocess.run(
    [
        "gmx",
        "grompp",
        "-p",
        "system.top",
        "-f",
        "../minimization.mdp",
        "-c",
        "antibody-CG.gro",
        "-o",
        "minimization-vac.tpr",
    ]
)
subprocess.run(["gmx", "mdrun", "-deffnm", "minimization-vac", "-v"], check=True)

# Insert mab molecules
subprocess.run(
    [
        "gmx",
        "insert-molecules",
        "-ci",
        mab,
        "-nmol",
        str(nmab),
        "-o",
        "mabs.gro",
        "-box",
        str(box_size),
    ],
    check=True,
)

# Solvate
subprocess.run(
    [
        "gmx",
        "solvate",
        "-cp",
        "mabs.gro",
        "-cs",
        water,
        "-radius",
        "0.21",
        "-o",
        "solvated.gro",
    ],
    check=True,
)

# Count water molecules
with open("solvated.gro") as f:
    nwater = sum(1 for line in f if line.startswith("W"))

# Update system.top
system_top_content = f"""\
#include "../martini_v3.0.0.itp"
#include "../martini_v3.0.0_solvents_v1.itp"

#include "molecule_0.itp"

[ system ]
; name
Martini system from antibody.pdb

[ molecules ]
; name        number
molecule_0    {nmab}
W    {nwater}
"""
with open("system.top", "w") as f:
    f.write(system_top_content)

# Copy solvated.gro to system.gro
os.rename("solvated.gro", "system.gro")


def run(epsilon_r):
    platform = None  # Platform.getPlatformByName("OpenCL")
    properties = {"Precision": "double"}

    conf = GromacsGroFile("system.gro")
    box_vectors = conf.getPeriodicBoxVectors()

    # get any defines
    defines = {}
    try:
        with open("defines.txt") as def_file:
            for line in def_file:
                line = line.strip()
                defines[line] = True
    except FileNotFoundError:
        pass

    top = martini.MartiniTopFile(
        "system.top",
        periodicBoxVectors=box_vectors,
        defines=defines,
        epsilon_r=epsilon_r,
    )
    system = top.create_system(nonbonded_cutoff=1.1 * nanometer)

    integrator = LangevinIntegrator(310 * kelvin, 10.0 / picosecond, 20 * femtosecond)
    integrator.setRandomNumberSeed(0)

    simulation = Simulation(top.topology, system, integrator, platform, properties)

    simulation.context.setPositions(conf.getPositions())

    # Print the platform being used
    platform_name = simulation.context.getPlatform().getName()
    print(f"Using platform: {platform_name}")
    ################################################################################
    ### Minimization ###

    simulation.reporters.append(PDBReporter("mini.pdb", 1000))
    simulation.reporters.append(
        StateDataReporter(
            "mini.log",
            5000,
            step=True,
            speed=True,
            potentialEnergy=True,
            temperature=True,
            volume=True,
        )
    )
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=5000, tolerance=1.0)

    energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print("System minimized at", energies, "\n")

    ################################################################################
    ### NVT equilibration ###

    simulation.context.setVelocitiesToTemperature(310 * kelvin)
    print("Running NVT equilibration...")
    simulation.step(50000)  # 1ns

    ################################################################################
    ### NPT equilibration ###

    system.addForce(MonteCarloBarostat(1 * bar, 310 * kelvin))
    # to update the simulation object to take in account the new system
    simulation.context.reinitialize(True)
    print("Running NPT equilibration...")
    simulation.step(50000)  # 1ns

    # save the equilibration results to file
    simulation.saveState("equi.state")
    simulation.saveCheckpoint("equi.chk")

    ################################################################################
    ### Production run ###

    # Set up the reporters to report energies every 1000 steps.
    simulation.reporters.append(
        StateDataReporter(
            "prod.log",
            1000,
            step=True,
            speed=True,
            potentialEnergy=True,
            totalEnergy=True,
            density=True,
            temperature=True,
            volume=True,
        )
    )
    # save the trajectory in XTC format
    xtc_reporter = XTCReporter("prod.xtc", 1000)
    simulation.reporters.append(xtc_reporter)

    # run simulation
    print("Running simulation...")
    simulation.step(250000)  # 5ns


run(15)
