#  This runs the tilted-v catchment problem
#  similar to that in Kollet and Maxwell (2006) AWR

#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin 
package require parflow
namespace import Parflow::*

pfset FileVersion 4

pfset Process.Topology.P 1
pfset Process.Topology.Q 1
pfset Process.Topology.R 1

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

pfset ComputationalGrid.NX                100
pfset ComputationalGrid.NY                1
pfset ComputationalGrid.NZ                300

pfset ComputationalGrid.DX	         1.0
pfset ComputationalGrid.DY               1.0
pfset ComputationalGrid.DZ	            .05

#---------------------------------------------------------
# Domain Geometry 
#---------------------------------------------------------
pfset GeomInput.Names                 "solidinput1 slabin clayin"

pfset GeomInput.solidinput1.InputType  SolidFile
pfset GeomInput.solidinput1.GeomNames  domain
pfset GeomInput.solidinput1.FileName   domain.pfsol

pfset Geom.domain.Patches             "z-upper x-lower y-lower \
                                      x-upper y-upper z-lower"

				  
pfset GeomInput.slabin.InputType  SolidFile
pfset GeomInput.slabin.GeomNames  slab
pfset GeomInput.slabin.FileName   slab.pfsol

pfset Geom.slab.Patches             "z-upper x-lower y-lower \
                                      x-upper y-upper z-lower"

pfset GeomInput.clayin.InputType       Box
pfset GeomInput.clayin.GeomName        clay

#-----------------------------------------------------------------------------
# Concen_Region Geometry
#-----------------------------------------------------------------------------
pfset Geom.clay.Lower.X   8.0
pfset Geom.clay.Lower.Y   0.0
pfset Geom.clay.Lower.Z   5.8

pfset Geom.clay.Upper.X   50.0
pfset Geom.clay.Upper.Y   1.0
pfset Geom.clay.Upper.Z    6.2

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------

pfset Geom.Perm.Names                 "domain slab clay"

pfset Geom.domain.Perm.Type            Constant
pfset Geom.domain.Perm.Value           10.

pfset Geom.clay.Perm.Type            Constant
pfset Geom.clay.Perm.Value           0.025

pfset Geom.slab.Perm.Type            Constant
pfset Geom.slab.Perm.Value           0.001

pfset Perm.TensorType               TensorByGeom

pfset Geom.Perm.TensorByGeom.Names  "domain"

pfset Geom.domain.Perm.TensorValX  1.0d0
pfset Geom.domain.Perm.TensorValY  1.0d0
pfset Geom.domain.Perm.TensorValZ  1.0d0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------

pfset SpecificStorage.Type            Constant
pfset SpecificStorage.GeomNames       "domain"
pfset Geom.domain.SpecificStorage.Value 0.0e-5

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------

pfset Phase.Names "water"

pfset Phase.water.Density.Type	        Constant
pfset Phase.water.Density.Value	        1.0

pfset Phase.water.Viscosity.Type	Constant
pfset Phase.water.Viscosity.Value	1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------

pfset Contaminants.Names			""

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------

pfset Geom.Retardation.GeomNames           ""

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------

pfset Gravity				1.0

#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------

# run for 2 hours @ 6min timesteps
# 
pfset TimingInfo.BaseUnit        0.1
pfset TimingInfo.StartCount      0
pfset TimingInfo.StartTime       0.0
pfset TimingInfo.StopTime        2.0
pfset TimingInfo.DumpInterval    0.1
pfset TimeStep.Type              Constant
pfset TimeStep.Value             0.1
 
#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames          "domain"


pfset Geom.domain.Porosity.Type          Constant
pfset Geom.domain.Porosity.Value         0.1

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------

pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "domain"

pfset Geom.domain.RelPerm.Alpha         1.0
pfset Geom.domain.RelPerm.N             2. 

#pfset Geom.slab.RelPerm.Alpha         1.0
#pfset Geom.slab.RelPerm.N             3. 

#pfset Geom.clay.RelPerm.Alpha         1.0
#pfset Geom.clay.RelPerm.N             3. 

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "domain"

pfset Geom.domain.Saturation.Alpha        1.0
pfset Geom.domain.Saturation.N            2.
pfset Geom.domain.Saturation.SRes         0.2
pfset Geom.domain.Saturation.SSat         1.0

#pfset Geom.slab.Saturation.Alpha        1.0
#pfset Geom.slab.Saturation.N            3.
#pfset Geom.slab.Saturation.SRes         0.3
#pfset Geom.slab.Saturation.SSat         1.0

#pfset Geom.clay.Saturation.Alpha        1.0
#pfset Geom.clay.Saturation.N            3.
#pfset Geom.clay.Saturation.SRes         0.3
#pfset Geom.clay.Saturation.SSat         1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant rainrec"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

# rainfall and recession time periods are defined here
# rain for 1 hour, recession for 2 hours

pfset Cycle.rainrec.Names                 "rain rec"
pfset Cycle.rainrec.rain.Length           4
pfset Cycle.rainrec.rec.Length            16
pfset Cycle.rainrec.Repeat                -1
 
#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames                   [pfget Geom.domain.Patches]

pfset Patch.x-lower.BCPressure.Type		      FluxConst
pfset Patch.x-lower.BCPressure.Cycle		      "constant"
pfset Patch.x-lower.BCPressure.alltime.Value	      0.0

pfset Patch.y-lower.BCPressure.Type		      FluxConst
pfset Patch.y-lower.BCPressure.Cycle		      "constant"
pfset Patch.y-lower.BCPressure.alltime.Value	      0.0

pfset Patch.z-lower.BCPressure.Type		      FluxConst
pfset Patch.z-lower.BCPressure.Cycle		      "constant"
pfset Patch.z-lower.BCPressure.alltime.Value	      0.0

pfset Patch.x-upper.BCPressure.Type		      FluxConst
pfset Patch.x-upper.BCPressure.Cycle		      "constant"
pfset Patch.x-upper.BCPressure.alltime.Value	      0.0

pfset Patch.y-upper.BCPressure.Type		      FluxConst
pfset Patch.y-upper.BCPressure.Cycle		      "constant"
pfset Patch.y-upper.BCPressure.alltime.Value	      0.0

## overland flow boundary condition with very heavy rainfall then slight ET
pfset Patch.z-upper.BCPressure.Type		      OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle		      "rainrec"
pfset Patch.z-upper.BCPressure.rain.Value	      -0.01
pfset Patch.z-upper.BCPressure.rec.Value	      0.00

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type "Constant"
pfset TopoSlopesX.GeomNames "domain"
pfset TopoSlopesX.Geom.domain.Value 0.1

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------


pfset TopoSlopesY.Type "Constant"
pfset TopoSlopesY.GeomNames "domain"
pfset TopoSlopesY.Geom.domain.Value 0.00

#---------------------------------------------------------
# Mannings coefficient 
#---------------------------------------------------------

pfset Mannings.Type "Constant"
pfset Mannings.GeomNames "domain"
pfset Mannings.Geom.domain.Value 1.e-6

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.water.Type                         Constant
pfset PhaseSources.water.GeomNames                    domain
pfset PhaseSources.water.Geom.domain.Value        0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------

pfset KnownSolution                                    NoKnownSolution


#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------

pfset Solver                                             Richards
pfset Solver.MaxIter                                     2500
pfset Solver.Drop                                      1E-20
pfset Solver.AbsTol                                     1E-9

pfset Solver.Nonlinear.MaxIter                           300
pfset Solver.Nonlinear.ResidualTol                       1e-7
pfset Solver.Nonlinear.EtaChoice                         Walker1 
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          0.001
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-16
pfset Solver.Nonlinear.StepTol				 1e-20
pfset Solver.Nonlinear.Globalization                     LineSearch
pfset Solver.Linear.KrylovDimension                      20
pfset Solver.Linear.MaxRestart                           2

#pfset Solver.Linear.Preconditioner                            PFMGOctree
#pfset Solver.Linear.Preconditioner.PFMGOctree.MaxIter          1
#pfset Solver.Linear.Preconditioner.PFMGOctree.MaxLevels       10
#pfset Solver.Linear.Preconditioner.PFMGOctree.MaxPreRelax      1
#pfset Solver.Linear.Preconditioner.PFMGOctree.MaxPostRelax     1
#pfset Solver.Linear.Preconditioner.PFMGOctree.BoxSizePowerOf2  1

pfset Solver.Linear.Preconditioner                            PFMG
#pfset Solver.Linear.Preconditioner.PCMatrixType               FullJacobian

pfset Solver.PrintSubsurf				False
pfset Solver.PrintMask                                  True
pfset Solver.PrintSlopes                                True
pfset Solver.PrintMannings                              True
 
#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------

# set water table to be at the bottom of the domain, the top layer is initially dry
pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.Value                      0.0

pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   z-lower

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------
pfwritedb slab3.1m
