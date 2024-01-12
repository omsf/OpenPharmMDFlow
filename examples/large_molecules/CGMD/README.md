# CGMD
CG models for Large Molecule Assets

![NivoCG Image](anitbody-CG.png)

## Install
First install martini_openmm by following instructions in `https://github.web.bms.com/PDModeling/martini_openmm/blob/master/tutorial/README.md`
Then,

```
conda install -c salilab dssp
conda install -c conda-forge gromacs
pip install vermouth==0.9.3
```

# Run
```
mkdir {new molecule name}
cd {new molecule name}
```
copy the new molecule pdb file into the folder. ** CALL THE PDB FILE antibody.pdb **

### For single antibody in vacuum
```
bash ../run_martinize.sh
python ../generate_gro_file.py
python ../run_CG_simulation.py
```

### For single antibody in a 20 nm water box
```
bash ../run_martinize.sh
python ../generate_gro_file.py
bash ../create_simulation.sh
python ../run_CG_simulation.py
```
You should see an output like below (rendered in VMD)
![PAD4 gif](antibody_in_water.gif)


### For 2 antibodies in a 20 nm water box
```
bash ../run_martinize.sh
python ../generate_gro_file.py
bash ../create_simulation.sh 2 20
python ../run_CG_simulation.py
```

## Run on cluster
```
git clone https://github.web.bms.com/PDModeling/CGMD.git
cd CGMD
docker build -t gromacs_openmm:latest .
docker tag gromacs_openmm bmset/gromacs_openmm:latest
docker push bmset/gromacs_openmm:latest
mkdir {new molecule name}
cd {new molecule name}
cp {path to antibody.pdb} ./
sbatch ../slurm_submit.sh
```
