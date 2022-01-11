#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --account=slts
#SBATCH --exclusive

export PARFLOW_DIR="/p/project/cslts/slts31/parflow_ben_original/run"

tclsh test1.tcl
srun $PARFLOW_DIR/bin/parflow testOne

tclsh drainageslab.tcl
srun $PARFLOW_DIR/bin/parflow drainageslab1

tclsh inflowslab.tcl
srun $PARFLOW_DIR/bin/parflow inflowslab1

tclsh parkinglot.tcl
srun $PARFLOW_DIR/bin/parflow parkinglot1

tclsh profileclm.tcl
srun $PARFLOW_DIR/bin/parflow profileclm1

tclsh slab.tcl
srun $PARFLOW_DIR/bin/parflow slab1

tclsh terrainfollowing.tcl
srun $PARFLOW_DIR/bin/parflow terrainfollowing1
