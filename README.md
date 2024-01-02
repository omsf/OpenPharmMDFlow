# OpenPharmMDFlow

Right now the project is at [phase 0](https://github.com/orgs/omsf/projects/2), where most of the work will be done by creating issues which we will organize on our phase 0 project board.

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

