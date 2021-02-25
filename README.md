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

**Usage:**  
The Diagnostics class is based on HeAT and requires methods to read ParFlow output (netCDF, PFB, Silo).
Any operations on results obtain by applying the Diagnostics functions, such filtering of active/inactive regions, spatial/temporal averaging must be implemented by the user.
Example applications are provided in the accompanying test cases, which require to run ParFlow via the tcl scripts to produce output.

Only use the Diagnostics in serial, ie using only one process. Parallelization is currently NOT supported and still in development.

**Test cases:**  
* test1.tcl; test.py: Simply 3D box domain, infiltration at the top patch, no-flow over all other patches  
* slab.tcl, drainageslab.tcl, inflowslab.tcl, slab.pfsol, domain.pfsol; slap.py: Slab benchmark + variations including drainage at the bottom and inflow over the righ patch, variations can be analyzed with slab.py by changing the 'name'  
* profileclm.tcl, geom.pfsol, atmforcing.txt, drv_clmin.dat, drv_vegm.dat, drv_vegp.dat; profileclm.py: A quasi 2D profile including exchange with the landsurface based on coupled CLM  
* parkinglot.tcl; parkinglot.py: Simple parking lot runoff test case for overland flow only; extracts a hydrograph at an individual pixel  
* terrainfollowing.tcl, terrainfollowing.py: Terrain following grid test case based in quasi 2D cross-section
* cordex3km.py: Example of a comprehensive mass balance calculation; requires access to the shared drive

