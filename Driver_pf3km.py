import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics
import IO as io
import glob
#import io_nc as nc
def _print(*args):
    """ Printing method for parallel environments. Only root prints """
    if ht.MPI_WORLD.rank is 0:
        print(*args, flush=True)

#Run ParFlow test case
#Export ParFlow install directory; needs to be changed by the user
#e.g. export PARFLOW_DIR=/home/s.kollet/restore/migrated/Programs/parflow/
#Oddly this can not be done with os.system

# This can be done with os.environ, e.g.:
# os.environ['PARFLOW_DIR'] = '/home/b.bourgart/parflow/install'

#cmd='tclsh test1.tcl'
#os.system(cmd)
#cmd='$PARFLOW_DIR/bin/parflow testOne'
#os.system(cmd)
name = 'cordexIHME'
path='/p/home/jusers/naz1/jureca/SCRATCH_jibg3103/parflow_3km/heat_analysis_1/output_27/'
hr_files = sorted(glob.glob(os.path.join(path, name +'.out.press.*.pfb')))
outpath='/p/home/jusers/naz1/jureca/SCRATCH_jibg3103/parflow_3km/heat_analysis_1/ana_parflow-diagnostics_pythonheat/output/'

split=None

dx = dy = 3000.0
dz = 2.0
#dz = 1.0
dzmult = [0.01,0.015,0.025,0.035,0.065,0.10,0.15,0.25,0.35,0.50,2.0,5.0,5.0,7.50,9.0] #List of dzScales
nx = 1592
ny = 1544
nz = 15
dt = ht.float64(1.0)

nt = len(hr_files)
_print(nt)
shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)



#Read static information
saturPF  = io.read_pfb(path + name + '.out.satur.00000.pfb',split)
sstorage = io.read_pfb(path + name + '.out.specific_storage.pfb',split)
mask     = io.read_pfb(path + name + '.out.mask.pfb',split)
poro     = io.read_pfb(path + name + '.out.porosity.pfb',split)
mannings = io.read_pfb(path + name + '.out.mannings.pfb',split)
#mannings = None
slopex   = io.read_pfb(path + 'ParFlow_MB3km_SLPX_x1592y1544.pfb',split)
#slopex   = None
slopey   = io.read_pfb(path + 'ParFlow_MB3km_SLPY_x1592y1544.pfb',split)
#slopey   = None
alpha = io.read_pfb(path + 'Alpha3D.pfb',split)
nvg   = io.read_pfb(path + 'Nvg3D.pfb',split)
sres  = io.read_pfb(path + 'Sres3D.pfb',split)


_print(saturPF.shape)
_print(sstorage.shape)
_print(mask.shape)
_print(poro.shape)
_print(mannings.shape)
_print(slopex.shape)
_print(slopey.shape)
_print(alpha.shape)
_print(nvg.shape)
_print(sres.shape)


split=2
#Set the mask to one in active and zero in inactive regions

mask  = ht.where(mask>0.0, 1.0, mask)

#Some other constant values
perm  = ht.full(shape3D,0.01,dtype=ht.float64,split=None)
ssat  = ht.full(shape3D,1.0,dtype=ht.float64,split=None)
#sres  = ht.full(shape3D,0.2,dtype=ht.float64,split=None)
#alpha = ht.full(shape3D,1.0,dtype=ht.float64,split=None)
#nvg   = ht.full(shape3D,2.0,dtype=ht.float64,split=None)
weekly_sm = ht.zeros(shape4D,split=split)
sub_storage = ht.zeros(shape4D,split=split)
lat_overlandflow = ht.zeros(shape3D,split=split)
#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, split)


for t, file in enumerate(hr_files):
	press = io.read_pfb(file, split=None)
	_print('processed file no:', t, '\tpath:', file)
	_print(press.shape)
	#if t>0:
	#   subsurface_storage_old = subsurface_storage

	#Calculate relative saturation and relative hydraulic conductivity
	satur,krel = diag.VanGenuchten(press)
	_print(satur.shape)

	#Returns an unmasked 3D field of subsurface storage, (L^3)
	subsurface_storage=diag.SubsurfaceStorage(press,satur)

	#Returns an unmasked 3D field of volumetric soil moisture, (L^3/L^3)
	volumetric_moisture = diag.VolumetricMoisture(satur)
	_print(volumetric_moisture.shape)
	weekly_sm[t] = volumetric_moisture

	#Calculate total subsurface storage, (L^3)
	subsurface_storage = subsurface_storage * mask
	#total_subsurface_storage = ht.sum(subsurface_storage)
	_print(subsurface_storage.shape)
	sub_storage[t] = subsurface_storage		

	#Return the top layer pressure 
	#top_layer = diag.TopLayerPressure(press) 
	#Calculate overland flow, (L^2/T)
	#overland_flow_x,overland_flow_y = diag.OverlandFlow(press,top_layer) 

	#Calculate net lateral overland flow
	#net_later_overland_flow = diag.NetLateralOverlandFlow(overland_flow_x,overland_flow_y)
	#_print(net_later_overland_flow.shape)
	#lat_overlandflow[t] = net_later_overland_flow
	

ht.save_netcdf(weekly_sm,outpath+name + '.out_soilmoisture_w01.nc','SM')
ht.save_netcdf(sub_storage,outpath+name + '.out_sub_storage_w01.nc','substor')
#ht.save_netcdf(lat_overlandflow,outpath+name + '.out_lateral_overlandflow_w01.nc','overflow')


    #Calculate subsurface flow in all 6 directions for each grid cell (L^3/T)
    #flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag.SubsurfaceFlow(press,krel)

    #Column balance
    #if t > 0:
    #  column_balance  = ht.zeros(shape2D,dtype=ht.float64,split=split)
      #Change in storage for each column
    #  column_balance  = (ht.sum(subsurface_storage_old - subsurface_storage,axis=0))
      #Add divergence of the flux for each column
    #  column_balance += dt * ht.sum(flowleft-flowright+flowfront-flowback+flowbottom-flowtop, axis=0)
      #Mass balance over full domain without flux at the top boundary
    #  _print('Time step:',t, ', Increment:',ht.sum(column_balance)/(nx*ny))

#_print()
#_print('At each time step, the increment should be equal to the flux at the top boundary: 0.0001 (L/T)')


#cmd='rm -r testOne*'
#os.system(cmd)
#my_dict = {'Sm': weekly_sm}
#io.write_nc4(my_dict, outpath+name + '.out_soilmoisture_01.nc') # not working error in ht.save_netcdf
