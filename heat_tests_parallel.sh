#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:10:00
#SBATCH --partition=batch
#SBATCH --account=slts
#SBATCH --exclusive
#SBATCH --output=parallel_080322.out
#SBATCH --error=parallel_080322.err
export OMP_NUM_THREADS=24

echo test1
srun -n 2 python test1.py

#srun python drainageslab.py

#srun python inflowslab1
echo
echo parkinglot
srun -n 2 python parkinglot.py

echo
echo profileclm
srun -n 2 python profileclm.py

echo
echo slab
srun -n 2 python slab.py

echo
echo terrainfollwing
srun -n 2 python terrainfollowing.py
