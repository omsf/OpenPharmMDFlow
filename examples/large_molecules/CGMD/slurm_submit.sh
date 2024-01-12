#!/bin/bash
#SBATCH --job-name=openmm_cg_cuda
#SBATCH --partition gpu
#SBATCH --gres=gpu:tesla:1
#SBATCH --exclude=gpu-dy-p38xlarge-1,gpu-dy-p38xlarge-2
#SBATCH --output=/scratch/job.%j.out

whoami
id -a
docker run --rm --runtime=nvidia -v $(dirname `pwd`):$(dirname `pwd`) -w `pwd` bmset/gromacs_openmm:latest bash ../run_martinize.sh
docker run --rm --runtime=nvidia -v $(dirname `pwd`):$(dirname `pwd`) -w `pwd` bmset/gromacs_openmm:latest python ../generate_gro_file.py
docker run --rm --runtime=nvidia -v $(dirname `pwd`):$(dirname `pwd`) -w `pwd` bmset/gromacs_openmm:latest bash ../create_simulation.sh 2 20
docker run --rm --runtime=nvidia -v $(dirname `pwd`):$(dirname `pwd`) -w `pwd` bmset/gromacs_openmm:latest python ../run_CG_simulation.py

