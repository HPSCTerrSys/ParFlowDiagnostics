#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:15:00
#SBATCH --partition=batch
#SBATCH --account=slts
#SBATCH --exclusive
#SBATCH --output=serial_compare.out
#SBATCH --error=serial.err
export OMP_NUM_THREADS=24

echo test1
srun python test1.py

#srun python drainageslab.py

#srun python inflowslab1
echo
echo parkinglot
srun python parkinglot.py

echo
echo profileclm
srun python profileclm.py

echo
echo slab
srun python slab.py

echo
echo terrainfollwing
srun python terrainfollowing.py
