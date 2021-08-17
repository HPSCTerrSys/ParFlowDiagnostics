# ParFlow Diagnostics

This Diagnositics class provides functions to calculated a global or local mass balance based on ParFlow output.

class Diagnostics  
def __init__(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split)
* Massk:            3D mask, tensor float (nz, ny, nx)
* Perm:             4D permeability, tensor float (3, nz, ny ,nx)
* Poro:             3D porosity, tensor float(nz, ny, nx)
* Sstorage:         3D specific storage, tensor float (nz, ny, nx)
* Ssat:             3D van Genuchten Ssat, tensor float (nz, ny, nx)
* Sres:             3D van Genuchten Sres, tensor float (nz, ny, nx)
* Nvg:              3D van Genuchten n, tensor float (nz, ny, nx)
* Alpha:            3D van genuchten alpha, tensor float (nz, ny, nx)
* Mannings:         quasi 2D Mannings n, tensor float (1, ny, nx)
* Slopex, Slopey:   quasi 2D slope in x- and y-direction, tensor float (1, ny, nx)
* Dx, Dy, Dz:       discretization in x-, y-, z-direction, float
* Nx, Ny, Nz:       number of grid points in x-, y-, z-direction, integer
* Terrainfollowing: flag (True/False) for terrain following grid, logical
* Split:            split axis for parallel computations

**Available Functions:**

* SubsurfaceStorage(self, Press, Satur) -> subsurface_storage  
Calculates water and compressible subsurface storage (L³) for each cell without distinguishing between active/inactive regions. Output is a 3D tensor.

* VolumetricMoisture(self, Satur) -> volumetric_moisture  
Calculates volumetric soil moisture (-) for each cell without distinguishing between active/inactive regions. Output is a 3D tensor.

* TopLayerPressure(self, Press) -> top_layer_pressure  
Extracts the pressure head (L) opf the cell at the surface taking into acount active/inactive regions. Output is a 2D tensor.

* SurfaceStorage(self,Toplayerpress) -> surface_storage  
Calculates water storage (L³) at the surface for each cell taking into account active/inactive regions. Output is a 2D tensor

* OverlandFlow(self,Toplayerpress) -> oflowx, oflowy  
Calculates overland flow (L²/T) in the x- and y-direction without distinguishing between active/inactive regions. Output is two 2D tensors. 

* NetLateralOverlandFlow(self, overland_flow_x, overland_flow_y) -> net_overland_flow  
Calculates net overland (L³/T) flow without distinguishing between active/inactive regions. Output is a 2D tensor.

* SubsurfaceFlow(self, Press, Krel) -> flowleft, flowright, flowfront, flowback, flowbottom, flowtop  
Calculates volume flux (L³/T) over all six cell faces taking into account active/inactive regions. Output is six 3D tensors.

* VanGenuchten(self,Press) -> satur, krel  
Calculated relative saturation (-) and relative conductivity (-) based on the van Genuchten function. Output is two 3D tensors.

**HeAT**
HeAT is a Python Array computation library, structured similar to numpy. Additional Features include GPU support and parallel processing capabilities for HPC systems by introducing a split axis. The Array is split to the mutliple processes along this axis. The Documentation as well as tutorials can be found at (https://github.com/helmholtz-analytics/heat).
On the HPC systems, working HeAT environments can be found at
'/p/project/cslts/local/juwels/HeAT/'
or
'/p/project/cslts/local/jureca/HeAT/'

To use the environment, source the most recent .ini file:
'source PythonHeAT_Stage2020_20200925_5fa39184_v0.5.0.ini'

**Usage:**  
The Diagnostics class is based on HeAT and requires methods to read ParFlow output (netCDF, PFB, Silo).
The Diagnostics class takes as input the split-axis for HeAT. Both splitting in x and y direction is supported (in a zyx coordinate system) as well as None for not splitting at all. In order to use mutliple nodes, the split axis must not be None. The job configuration allows basically arbitrary MPI and OpenMP parallelisation, however it is recommended to use 2 MPI-processes per Node and set the number of OpenMP threads accordingly, i.e. #number_of_CPU_cores_per_node / 2  (24 on JUWELS, 64 on JurecaDC).
To use the Diagnostics methods, clone the gitlab repository into your working directory:
'git clone https://icg4geo.icg.kfa-juelich.de/SoftwareTools/ana_parflow-diagnostics_pythonheat.git <WORKING_DIR>'
Then, you can use it in Python via
'from Diagnostics import Diagnostics'
Any operations on results obtain by applying the Diagnostics functions, such filtering of active/inactive regions, spatial/temporal averaging must be implemented by the user.
Example applications are provided in the accompanying test cases, which require to run ParFlow via the tcl scripts to produce output.

**Tests**
* All Methods are tested using both no split and splitting across multiple nodes
* All Methods iterating over the input array have been vectorized at a significant performance increase and been verified to yield the same result
* At the beginning of the development, several issues arose with HeAT (especially indexing), those are fixed as of HeAT version 1.1, do not use these Diagnostics with older versions of HeAT.

**Test cases:**  
* test1.tcl; test.py: Simply 3D box domain, infiltration at the top patch, no-flow over all other patches  
* slab.tcl, drainageslab.tcl, inflowslab.tcl, slab.pfsol, domain.pfsol; slap.py: Slab benchmark + variations including drainage at the bottom and inflow over the righ patch, variations can be analyzed with slab.py by changing the 'name'  
* profileclm.tcl, geom.pfsol, atmforcing.txt, drv_clmin.dat, drv_vegm.dat, drv_vegp.dat; profileclm.py: A quasi 2D profile including exchange with the landsurface based on coupled CLM  
* parkinglot.tcl; parkinglot.py: Simple parking lot runoff test case for overland flow only; extracts a hydrograph at an individual pixel  
* terrainfollowing.tcl, terrainfollowing.py: Terrain following grid test case based in quasi 2D cross-section
* cordex3km.py: Example of a comprehensive mass balance calculation; requires access to the shared drive

**Big-Data Benchmarking**
Generic Test for Big-Data capability and Scalability of all Methods provided in the Diagnostics Repository using the Global 1km Dataset by Stefan Kollet ( /p/scratch/cesmtst/kollet1/globaltest/ ).
Then using 4, 8, 16 nodes on juwels yielded the following performance-scaling plots:

These Scaling experiments are generated using the JUBE2 - Benchmark environment (available via module load JUBE) and the provided JUBE_diagnostics.xml file, which defines the Job configurations and runs the benchmark.py file. In order to run your own benchmarks, you can find more information here: https://icg4geo.icg.kfa-juelich.de/SoftwareTools/ana_BigDataAnalytics_PythonHeAT/-/tree/jube_benchmarking 
For visualization of jubes outputs, the jube_result.ipynb is available. It generates multiple plots showing the runtimes and speedups of the different sections. To run it, use the JSC Jupyter system (https://jupyter-jsc.fz-juelich.de/).

To test for the correctness of the parallel computations, a Reference Output was generated using only one MPI-process on JurecaDC ( /p/home/jusers/bourgart1/juwels/cesmtst_bourgart1/reference_output ) 
and then compared to calculate the difference. 
The difference introduced by applying parallelization is summarized by the follwing:
