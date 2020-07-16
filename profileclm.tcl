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

pfset ComputationalGrid.NX               100 
pfset ComputationalGrid.NY                 5
pfset ComputationalGrid.NZ               600

pfset ComputationalGrid.DX	         1.0
pfset ComputationalGrid.DY               0.2
pfset ComputationalGrid.DZ	         0.01

set nx [pfget ComputationalGrid.NX]
set dx [pfget ComputationalGrid.DX]
set ny [pfget ComputationalGrid.NY]
set dy [pfget ComputationalGrid.DY]
set nz [pfget ComputationalGrid.NZ]
set dz [pfget ComputationalGrid.DZ]

#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------
pfset GeomInput.Names                 "solidinput"

pfset GeomInput.solidinput.InputType  SolidFile
pfset GeomInput.solidinput.GeomNames  domain
pfset GeomInput.solidinput.FileName   geom.pfsol

pfset Geom.domain.Patches             "top bottom perimeter"

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------
pfset Geom.Perm.Names                 domain

# Values in m/hour

pfset Geom.domain.Perm.Type            Constant
pfset Geom.domain.Perm.Value           1.0

#-----------------------------------------------------------------------------
# Perm Tensors
#-----------------------------------------------------------------------------

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
pfset Geom.domain.SpecificStorage.Value 0.0001
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
pfset TimingInfo.BaseUnit               1.0
pfset TimingInfo.StartCount             0
pfset TimingInfo.StartTime              0.0
pfset TimingInfo.StopTime               24.0
pfset TimingInfo.DumpInterval           1.0
pfset TimeStep.Type                     Constant
pfset TimeStep.Value                    1.0

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames           domain

pfset Geom.domain.Porosity.Type          Constant
pfset Geom.domain.Porosity.Value         0.415

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------

pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------


pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          domain

pfset Geom.domain.RelPerm.Alpha         1.0
pfset Geom.domain.RelPerm.N             2. 

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         domain

pfset Geom.domain.Saturation.Alpha        1.0
pfset Geom.domain.Saturation.N            2.0
pfset Geom.domain.Saturation.SRes         0.11
pfset Geom.domain.Saturation.SSat         1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names constant
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1
 
#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames        "top bottom perimeter"

pfset Patch.perimeter.BCPressure.Type                      FluxConst
pfset Patch.perimeter.BCPressure.Cycle                     "constant"
pfset Patch.perimeter.BCPressure.alltime.Value             0.0

pfset Patch.bottom.BCPressure.Type                      FluxConst
pfset Patch.bottom.BCPressure.Cycle                     "constant"
pfset Patch.bottom.BCPressure.alltime.Value             0.0

pfset Patch.top.BCPressure.Type                       OverlandFlow
pfset Patch.top.BCPressure.Cycle                      "constant"
pfset Patch.top.BCPressure.alltime.Value              0.0

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type "Constant"
pfset TopoSlopesX.GeomNames "domain"
pfset TopoSlopesX.Geom.domain.Value 0.01

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------

pfset TopoSlopesY.Type "Constant"
pfset TopoSlopesY.GeomNames "domain"
pfset TopoSlopesY.Geom.domain.Value 0.01

#---------------------------------------------------------
# Mannings coefficient 
#---------------------------------------------------------

pfset Mannings.Type "Constant"
pfset Mannings.GeomNames "domain"
pfset Mannings.Geom.domain.Value 5.52e-6

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
pfset Solver.MaxIter                                     2500000

#pfset Solver.TerrainFollowingGrid                        True

pfset Solver.Nonlinear.MaxIter                           300
pfset Solver.Nonlinear.ResidualTol                       1.0e-8
pfset Solver.Nonlinear.EtaChoice                         Walker1 
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          0.001
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-16
pfset Solver.Nonlinear.StepTol				 1e-10
pfset Solver.Nonlinear.Globalization                     LineSearch
pfset Solver.Linear.KrylovDimension                      40
pfset Solver.Linear.MaxRestart                           2

pfset Solver.Linear.Preconditioner                       PFMGOctree

#pfset Solver.Linear.Preconditioner                       PFMG
#pfset Solver.Linear.Preconditioner.PCMatrixType          FullJacobian

#pfset Solver.Linear.Preconditioner.PCMatrixType     FullJacobian

#pfset Solver.Linear.Preconditioner                       MGSemi
#pfset Solver.Linear.Preconditioner.MGSemi.MaxIter        1
#pfset Solver.Linear.Preconditioner.MGSemi.MaxLevels      100

pfset Solver.Drop                                      1E-20
pfset Solver.AbsTol                                     1E-9

pfset Solver.LSM                                         CLM
pfset Solver.WriteSiloCLM                                True
pfset Solver.CLM.Print1dOut                              True
pfset Solver.BinaryOutDir                                ./binout 
pfset Solver.CLM.CLMDumpInterval                         8760

pfset Solver.CLM.MetForcing                              1D
pfset Solver.CLM.MetFileName                             atmforcing.txt
pfset Solver.CLM.MetFilePath                             ./
pfset Solver.CLM.IstepStart                              1

pfset Solver.PrintPressure                              True 
pfset Solver.PrintSaturation                            True 
pfset Solver.PrintSubsurfData                           True 
pfset Solver.PrintLSMSink                               True 
pfset Solver.PrintEvapTransSum                          True 
pfset Solver.PrintSlopes                                True 
pfset Solver.PrintMannings                              True 
pfset Solver.PrintMask                                  True 
pfset Solver.PrintOverlandSum                           False

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------

pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.Value                      -0.125

pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   top

#pfdist press.ini.pfb
#pfset ICPressure.Type                                   PFBFile
#pfset Geom.domain.ICPressure.FileName                   press.ini.pfb
#pfset ICPressure.GeomNames                              domain

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------

pfwritedb profileclm1

