#!/bin/bash -x
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --time=2:00:00
#SBATCH --partition=mem512
##SBATCH --partition=devel
#SBATCH --account=jibg31


OMP_NUM_THREADS=2 srun python Driver_pf3km.py 

