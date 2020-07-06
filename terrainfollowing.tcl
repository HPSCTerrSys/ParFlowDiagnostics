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

pfset ComputationalGrid.NX                10
pfset ComputationalGrid.NY                1
pfset ComputationalGrid.NZ                10

pfset ComputationalGrid.DX	         1.0
pfset ComputationalGrid.DY               1.0
pfset ComputationalGrid.DZ	            .1

#---------------------------------------------------------
# Domain Geometry 
#---------------------------------------------------------
pfset GeomInput.Names                 "domain_input"


pfset GeomInput.domain_input.InputType       Box
pfset GeomInput.domain_input.GeomName        domain

pfset Geom.domain.Patches             "x-lower x-upper y-lower y-upper z-lower z-upper"

#-----------------------------------------------------------------------------
# slab 1
#-----------------------------------------------------------------------------

pfset Geom.domain.Lower.X   0.0
pfset Geom.domain.Lower.Y   0.0
pfset Geom.domain.Lower.Z   0.0

pfset Geom.domain.Upper.X   10.0
pfset Geom.domain.Upper.Y    1.0
pfset Geom.domain.Upper.Z    1.0

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------

pfset Geom.Perm.Names                 "domain"

# Values in m/hour
pfset Geom.domain.Perm.Type            Constant
pfset Geom.domain.Perm.Value           100.

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
pfset Geom.domain.SpecificStorage.Value 0.0

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
pfset TimingInfo.StopTime        12.0
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
pfset Phase.RelPerm.GeomNames          domain

pfset Geom.domain.RelPerm.Alpha         2.0
pfset Geom.domain.RelPerm.N             2. 

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         domain

pfset Geom.domain.Saturation.Alpha        2.0
pfset Geom.domain.Saturation.N            2.
pfset Geom.domain.Saturation.SRes         0.2
pfset Geom.domain.Saturation.SSat         1.0

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
pfset Cycle.rainrec.rain.Length           20
pfset Cycle.rainrec.rec.Length            100
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
pfset Patch.z-upper.BCPressure.rain.Value	      -0.05
pfset Patch.z-upper.BCPressure.rec.Value	      0.00

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type "Constant"
pfset TopoSlopesX.GeomNames "domain"
pfset TopoSlopesX.Geom.domain.Value 0.5

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------


pfset TopoSlopesY.Type "Constant"
pfset TopoSlopesY.GeomNames "domain"
pfset TopoSlopesY.Geom.domain.Value 0.0

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

pfset Solver.TerrainFollowingGrid                        True

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

pfset Solver.Linear.Preconditioner                            PFMGOctree
pfset Solver.Linear.Preconditioner.PFMGOctree.MaxIter          1
pfset Solver.Linear.Preconditioner.PFMGOctree.MaxLevels       10
pfset Solver.Linear.Preconditioner.PFMGOctree.MaxPreRelax      1
pfset Solver.Linear.Preconditioner.PFMGOctree.MaxPostRelax     1
pfset Solver.Linear.Preconditioner.PFMGOctree.BoxSizePowerOf2  1

#pfset Solver.Linear.Preconditioner                            PFMG
#pfset Solver.Linear.Preconditioner.PCMatrixType               FullJacobian

pfset Solver.PrintSubsurf				True
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
pfwritedb terrainfollowing1
