Authors: Stefan KOLLET (s.kollet@fz-juelich.de), Ben BOURGART b.bourgart@fz-juelich.de), Klaus GOERGEN (k.goergen@fz-juelich.de)

# ParFlow Diagnostics (PD)

ParFlow Diagnostics (PD) is a Python class that provides functions to calculate all variables needed for a global or local mass balance based on [ParFlow hydrologic model](https://www.parflow.org) output. An essential feature is that PD uses the [Helmholtz Analytics Toolkit (HeAT)](https://github.com/helmholtz-analytics/heat/) Python library. This makes PD big data-capable; with HeAT PD can be run in parallel on a single HPC node (or any multi-core machine, including notebooks) up many nodes of an HPC system using CPUs as well as GPUs.

For details on HeAT, see the respective information referenced on the [HeAT github repository](https://github.com/helmholtz-analytics/heat/). This is not a HeAT tutorial, but using PD one also learns the basics of using HeAT. **HeAT is a totally generic library, it is only used here, similarly to using numpy, to enable our PD calculations to be run in parallel. It is currently recommended to use the PD not in parallel and only use a single process.**

**All information and software tools to use PD are within this repository. This README.md file contains all information to get started with PD including simple cookbooks.**

**Important prerequisite: The PD assume you use HeAT >v1.1. And, as of Autumn 2021, the experimental HeAT version `experimental_heat` provided on JUWELS/JURECA. If you follow this README this is ensured.**

**The `parallel_tests` branch is the current git branch of PD.**

**All tests and scripts of Autumn 2021 refer to the application on CPUs. GPU tests are ongoing.**

This repository covers six related aspects of using and testing PD and not every user needs to cover them all; the information provided here addresses novel HPC users as well as experienced users who just want to start using PD with HeAT:

1. PD description
2. How to set up PD based on HeAT (relevant for everyone; for experienced ParFlow users, who already have ParFlow output, this is the most important information)
3. How to generate test input for PD and example scripts on using PD, which can be used as a starting point for own developments (for beginners this is especially interesting as this repository contains ParFlow test cases, i.e., the model configuration to generate such simulation results for testing; also experienced users find the necessary PD examples here)
4. Tests conducted with PD based on HeAT (this is for those who want deeper understanding in how it was made sure, that the HeAT implementation did not break PD); these test cases can be run by everybody to generate their own test data
5. Specific performance tests using large datasets giving proof that HeAT in fact makes PD big data capable (this is for developers, who might be expanding PD or use PD based on HeAT for own developments); these test cases cannnot be run by everybody, instead test data, i.e., ParFlow simulation results, are provided centrally for everybodies use
6. Correctness tests ensure the developments do not break the tools (this is for developers who might want to do their own developments with PD and HeAT)

ParFlow and PD and HeAT run literally on any Linux system. However the instructions below are meant for the primary HPC systems at the Jülich Supercomputing Centre (JSC); as of Autumn 2021 this is [JURECA](https://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JURECA/JURECA_node.html) and [JUWELS](https://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUWELS/JUWELS_node.html).

## PD description

### Calling PD

```python
class Diagnostics  
def __init__(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split)
```
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

### Available functions

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

### HeAT

HeAT is a Python Array computation library, structured similar to numpy. Additional Features include GPU support and parallel processing capabilities for HPC systems by introducing a split axis. The Array is split to the mutliple processes along this axis. The Documentation as well as tutorials can be found at [https://github.com/helmholtz-analytics/heat](https://github.com/helmholtz-analytics/heat).

## Getting-started, part 1: Get PD and prepare for usage

1. Login to any of the JURECA or JUWELS front nodes, e.g., `ssh -i <secret_ssh_key> <username>@jureca.fz-juelich.de` or `ssh -i <secret_ssh_key> <username>@juwels.fz-juelich.de`
2. To use the PD tools one needs to have a working HeAT implementation. as HeAT cannot be loaded at the moment at JSC using the `module` command from a system-wide installation, the PD maintainers provide a working HeAT environment as part of the HPSC TerrSys and SDLTS account on JSC systems. Working HeAT environments are installed and maintained for everybodies use under: `/p/project/cslts/local/juwels/HeAT/` (for JUWELS), and `/p/project/cslts/local/jureca/HeAT/` (for JURECA). In order to use the PD with all its functionalities as of Autumn 2021, one needs to use an experimental installation of HeAT. To load a virtual Python environment, which includes everything needed to use PD based on the latest experimental HeAT, just do this:
On JURECA:
`cd /p/project/cslts/local/jureca/HeAT`
`source experimental_HeAT.ini`
On JUWELS:
`cd /p/project/cslts/local/jureca/HeAT`
`source experimental_heat.ini`
Note the command line prompt changes,
You can copy these software environment initialisation files everwhere you like on the JSC machines and source them from there. As a normal user, use just use the software environments recommended here and PD should work.
3. Get yourselves a recent git version: `module load git`
4. Before using the PD methods, clone the gitlab repository into your working directory: `git clone https://icg4geo.icg.kfa-juelich.de/SoftwareTools/ana_parflow-diagnostics_pythonheat.git <WORKING_DIR>`. `cd <WORKING_DIR>`; then check all branches `git branch -a` but only get and checkout `parallel_tests` branch: `git checkout -b parallel_tests remotes/origin/parallel_tests`
5. Then, you would use PD in your Python script or in an interactive Python session via: `from Diagnostics import Diagnostics`. A number of example uses of PD based on HeAT are provided in this repository. Most of these usage examples comnsist of (i) the ParFlow configuration (tcl file) to *generate* ParFlow simulation output by running a short ParFlow test and (ii) to *analyse* this output using PD with HeAT.

**RECOMMENDATION: Go to the next section and check out the ParFlow test cases and the accomanying PD/HeAT Python analysis scripts.**

Important notes:
- There are multiple HeAT installations under `/p/project/cslts/local/{jureca,juwels}/HeAT`; each of these insdtallations comes with its own `ini`-file; these different installations feature different HeAT versions (see HeAT version tag and also HeAT git commit -- in case it is a not a release version) and have been installed at different dates using (potentially) different underlying software environments, i.e., toolchains (see the [JSC information](https://apps.fz-juelich.de/jsc/hps/jureca/software-modules.html) on the toolchains and software environments). Just follow advice in this README.md on which HeAT environment to load.
- The Diagnostics class is based on HeAT and requires methods to read ParFlow output (netCDF, PFB, Silo).
- The Diagnostics class takes as input the split-axis for HeAT. Both splitting in x and y direction is supported (in a zyx coordinate system) as well as None for not splitting at all. In order to use mutliple nodes, the split axis must not be None. The job configuration allows basically arbitrary MPI and OpenMP parallelisation, however it is recommended to use <!--- 2 MPI-processes per compute node and set the number of OpenMP threads accordingly, i.e. #number_of_CPU_cores_per_node / 2  (24 on JUWELS, 64 on JurecaDC). Tests have shown that, if one is using HeAT (and hence also PD tools) in a parallel setup on a single multi-core node or accross multiple nodes, the best performance is reached when running with one MPI task per CPU and mutpiple OpenMP threads on the individual CPU cores. --->**only one MPI-process with the PD.**
 The correct way to do this in the SBATCH script is by setting the according SLURM parameters:
```
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=<#number_of_CPU_cores_per_node>
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
```
More Details can be found at [writing-a-batch-script](https://apps.fz-juelich.de/jsc/hps/juwels/batchsystem.html#writing-a-batch-script).

- Any operations on results obtained by applying the diagnostics functions, such filtering of active/inactive regions, spatial/temporal averaging must be implemented by the user. I.e, PD is enabling users to start with the analyses quickly and provides therefore a common computational kernel, which should be used by all to aslo ensure comparable and reproducible results. Albeit embedding PD into your own workflows still needs to be done individually.
- Other usage examples of PD with HeAT are included also in [SLOTH](https://icg4geo.icg.kfa-juelich.de/SoftwareTools/SLOTH).

## Getting-started, part 2: Specific ParFlow test simulations and test applications of PD

**Important prerequisite: These test simulations require an existing ParFlow installation, see below.**

Most of the Python scripts (mentioned last in the list of files) contains system calls to run ParFlow; this means each HeAT analysis test can generate on-the-fly its own test data. BUT: because running ParFlow to generate test data requires a different software environment than doing the analysis with PD and HeAT, this functionality has been deactivated. Hence: (i) first run the short ParFlow tests, (ii) do the analysis. To help in doing this, we created some wrapper scripts to do (i) and (ii) combined, see below.

Overview of existing test cases, ParFlow runs and analysis:

* `test1.tcl`; `test.py`: Simply 3D box domain, infiltration at the top patch, no-flow over all other patches  
* `slab.tcl`, `drainageslab.tcl`, `inflowslab.tcl`, `slab.pfsol`, `domain.pfsol`; `slap.py`: Slab benchmark + variations including drainage at the bottom and inflow over the righ patch, variations can be analyzed with slab.py by changing the 'name'  
* `profileclm.tcl`, `geom.pfsol`, `atmforcing.txt`, `drv_clmin.dat`, `drv_vegm.dat`, `drv_vegp.dat`; `profileclm.py`: A quasi 2D profile including exchange with the landsurface based on coupled CLM
* `parkinglot.tcl`; `parkinglot.py`: Simple parking lot runoff test case for overland flow only; extracts a hydrograph at an individual pixel  
* `terrainfollowing.tcl`, `terrainfollowing.py`: Terrain following grid test case based in quasi 2D cross-section
* `cordex3km.py`: Special analysis code example of a comprehensive mass balance calculation for a pan-European 3km model domain; requires access to the shared drive where the model output is located
* `Driver.py`: Generic test, also demonstrating the use of PD

1. To run the tests, a working ParFlow installation including a software environment (i.e., loaded modules and paths) is needed. Instructions on how to build ParFlow can be found at [https://icg4geo.icg.kfa-juelich.de/ModelSystems/ParFlow_scripts](https://icg4geo.icg.kfa-juelich.de/ModelSystems/ParFlow_scripts). Because the ParFlow installation is seperate from the HeAT and the PD, the source code is separated as well. The ParFlow software environment initialisation is contained with the ParFlow built scripts.
2. The `tests.sh` script runs the ParFlow elements of the test cases from above and hence generates data as input to PD; `tests.sh` (see also the `tests.py` in order to run interactively on compute nodes) needs to be called from the ParFlow environment.
3. The Python files that run the HeAT based PD tool need to be called from the HeAT environment: `source /p/project/cslts/local/jureca/HeAT/experimental_HeAT.ini` and then run `sbatch heat_tests_serial.sh` or `sbatch heat_tests_parallel.sh`; please note: depending on which tests you do (the split axis in the Python codes called has to be set accordingly, basically "activating the parallelism" in HeAT!).

The ParFlow and the PD software environments can be switched by `sourcing` their respective `.ini` files, `source <inifile>`.

For simplicity: Just use two seperate sessions, one to run the ParFlow test simulations, the other one to use PD with HeAT.

## Advanced usage: PD standard tests

Note: Only relevant for users with HPC expertise.

During the development phase, a base version of PD was used (non-vectorized and serial) and (i) vectorized and (ii) efficiency-improved using HeAT. Estensive tests were done to ensure that both performance improvements of (i) and (ii) do not alter the results.

* All Methods have been tested using both "no split", "splitting on a single node", and "splitting across multiple nodes"; i.e. reference results of serial PD tests are compared as part of testing with parallel test results (see the correctness tests section below).
* All Methods iterating over an input array have been vectorized at a significant performance increase and been verified to yield the same result.
As of the merge on 11.03.2022, only the vectorized versions are provided in the `master` branch. If you need to recover the non-vectorized versions, please checkout an earlier commit from the master branch or the `parallel_tests` branch.
* At the beginning of the development, several issues arose with HeAT (especially indexing), those are fixed as of HeAT version 1.1. See the revision of HeAT in the [HeAT github](https://github.com/helmholtz-analytics/heat/) for details.

## Advanced usage: PD big-data tests and benchmarking

Note: Only relevant for users with HPC expertise.

The big data-capability and the scalability of all PD methods have been tested to make sure that PD can really scale to many compute nodes on the JURECA and JUWELS HPC systems and thereby allow for an efficient use of big data sets in the multi-TB-range.

- The test data set is the existing global 1km ParFlow test by Stefan KOLLET under `/p/scratch/cesmtst/kollet1/globaltest/`.
- The scaling tests were done using 4, 8, 16 compute nodes on JUWELS.
- The scaling plots are under `benchmarking/plots`.

Results clearly show HeAT's capability to efficiently handle datasets too large to fit onto a single node.

To repeat this yourself:
- All Files for these Experiments can be found in the `benchmarking` directory.
- These scaling experiments are generated using the JUBE2 - Benchmark environment (available via `module load JUBE` on JSC machines) and the provided `JUBE_diagnostics.xml` file, which defines the job configurations and runs the `benchmark.py` file.
- In order to run your own benchmarks, you can find more information here: [https://icg4geo.icg.kfa-juelich.de/SoftwareTools/ana_BigDataAnalytics_PythonHeAT/-/tree/jube_benchmarking](https://icg4geo.icg.kfa-juelich.de/SoftwareTools/ana_BigDataAnalytics_PythonHeAT/-/tree/jube_benchmarking)
- For visualization of JUBE's outputs, the `jube_result.ipynb` is available. It generates multiple plots showing the runtimes and speedups of the different sections. To run it, use the JSC Jupyter system [https://jupyter-jsc.fz-juelich.de/](https://jupyter-jsc.fz-juelich.de/).

## Correctness tests

Note: Only relevant for users with HPC expertise.

To test for the correctness of the parallel computations, a reference output was generated using only one MPI-process on JURECA-DC (`/p/home/jusers/bourgart1/juwels/cesmtst_bourgart1/reference_output`). The refernece was then compared with outputs based on a parallel computation. The difference introduced by applying parallelization is summarized by the following files: `diff.py`, `diff.out`, `diffs.png`.
