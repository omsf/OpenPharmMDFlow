[![CI](https://github.com/omsf/OpenPharmMDFlow/actions/workflows/ci.yaml/badge.svg)](https://github.com/omsf/OpenPharmMDFlow/actions/workflows/ci.yaml)
[![Documentation Status](https://readthedocs.org/projects/openpharmmdflow/badge/?version=latest)](https://openpharmmdflow.readthedocs.io/en/latest/?badge=latest)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/omsf/OpenPharmMDFlow/main.svg)](https://results.pre-commit.ci/latest/github/omsf/OpenPharmMDFlow/main)
[![codecov](https://codecov.io/gh/omsf/OpenPharmMDFlow/graph/badge.svg?token=SM16TD8IZ1)](https://codecov.io/gh/omsf/OpenPharmMDFlow)

# OpenPharmMDFlow

Right now the project is at [phase 0](https://github.com/orgs/omsf/projects/2), where most of the work will be done by creating issues which we will organize on our phase 0 project board.

## Install

See our docs [here](https://openpharmmdflow.readthedocs.io/en/latest/install.html).

### Small Molecule Pipeline Example

See our small molecule pipeline example notebook [here](https://github.com/omsf/OpenPharmMDFlow/blob/main/examples/pipeline/Small_Molecule_Pipeline.ipynb)
or view the notebook rendered as html in our docs page [here](https://openpharmmdflow.readthedocs.io/en/latest/small_molecule_pipeline_example.html).

## Manifest

```
├── inputs
│   ├── mAb  # Inputs for mAb workflows (excipients go here too)
│   └── small_molecule  # Inputs for small molecule workflows
└── README.md
```



## Small Molecule Molecular Mechanics Optimizer
`MM_opt_cif.py` uses openFF, openMM, and other cheminformatics tools to optimize a crystal in `.cif` format with molecular mechanics. Care has NOT been taken to make the script robust ; toggles for custom force fields (if available) are hard coded in lines 313-320 in `__main__`.

Representative `.cifs` are in the `cif/` directory. These outline a set of easy (ACETAC01, HXACAN), medium (RESORA A/B, SEDNOG), and difficult (UNISUG) cases. Difficult cases are determined by higher Z' and inhomogeneity. This is a charge transfer cocrystal that is NOT correctly represented by the forcefield (see opted structure).

The script does several steps including the identification of molecular units, assignment of force field parameters, standardization of the unit cell, and expansion of the unit cell (supercell) for nonbonded potential. Final sstructures are in `minimized/` then subsequently symmetrized with the `symmetrize.py` script.

`symmetrize.py` uses `spglib` to symmetrize and standardize the final crystal.
