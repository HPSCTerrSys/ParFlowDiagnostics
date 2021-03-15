import os

os.environ['PARFLOW_DIR'] = '/p/project/cslts/slts31/parflow_ben_original/run'

cmd='tclsh test1.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow testOne'
os.system(cmd)

cmd='tclsh drainageslab.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow drainageslab1'
os.system(cmd)

cmd='tclsh inflowslab.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow inflowslab1'
os.system(cmd)

cmd='tclsh parkinglot.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow parkinglot1'
os.system(cmd)

cmd='tclsh profileclm.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow profileclm1'
os.system(cmd)

cmd='tclsh slab.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow slab1'
os.system(cmd)

cmd='tclsh terrainfollowing.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow terrainfollowing1'
os.system(cmd)