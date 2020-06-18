import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics 
import IO as io

#Run ParFlow test case
#Export ParFlow install directory; needs to be changed by the user
cmd='export PARFLOW_DIR=/home/s.kollet/restore/migrated/Programs/parflow/'
os.system(cmd)
cmd='tclsh slab.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow slab3.1m'
os.system(cmd)


split=None
name = 'slab3.1m'
#Read static information
saturPF  = io.read_pfb(name + '.out.satur.00000.pfb',split)
sstorage = io.read_pfb(name + '.out.specific_storage.pfb',split)
mask     = io.read_pfb(name + '.out.mask.pfb',split)
poro     = io.read_pfb(name + '.out.porosity.pfb',split)
perm     = io.read_pfb(name + '.out.perm_x.pfb',split)

dx = dy = 1. 
dz = 0.05
nx = 100
ny = 1
nz = 300
dzmult = ht.full(nz,1.0,split=None)

dt = .1
nt = 20

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)

#Set the mask to one in active and zero in inactive regions
mask  = ht.where(mask>0.0, 1.0, mask)
#Some other constant values
ssat  = ht.full(shape3D,1.0,split=None)
sres  = ht.full(shape3D,0.2,split=None)
alpha = ht.full(shape3D,1.0,split=None)
nvg   = ht.full(shape3D,2.0,split=None)
#The following fields are special; while they are 2D, we need to define them in 3D,
#because they are often read from pfb, which is always 3D (with 1 in z-direction)
mannings = ht.full((1,ny,nx),0.000001,split=None)
slopex   = ht.full((1,ny,nx),0.1,split=None)
slopey   = ht.full((1,ny,nx),0.0,split=None)

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, split)

for t in range (nt):
    print(name + '.out.press.'+('{:05d}'.format(t))+'.pfb')
    press    = io.read_pfb(name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split=split)

    #Calculate relative saturation and relative hydraulic conductivity
    satur,krel = diag.VanGenuchten(press)

    #Returns an unmasked 3D field of subsurface storage, (L^3)
    subsurface_storage=diag.SubsurfaceStorage(press,satur)

    #Calculate subsurface flow in all 6 directions for each grid cell (L^3/T)
    flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag.SubsurfaceFlow(press,krel)

    #Obtain pressure at the land surface
    top_layer_press = diag.TopLayerPressure(press)
    
    #Calculate overland flow
    flowx,flowy = diag.OverlandFlow(top_layer_press)

    #Calculate total storage (subsurface + surface); here we need to use the mask, because domain is based on pfsol
    total_storage_column = ht.sum(subsurface_storage * mask,axis=0) + ht.where(top_layer_press>0.0, top_layer_press, 0.0) * dx * dy

    #Column balance
    if t > 0:
      #Change in storage for each column
      dstorage_column = (total_storage_column_old - total_storage_column)
      print(dstorage_column)
      #Add divergence of the flux for each column
      #column_balance += dt * ht.sum(flowleft-flowright+flowfront-flowback+flowbottom-flowtop, axis=0)
      #Mass balance over full domain without flux at the top boundary
      print('Time step:',t, ', Increment:',ht.sum(dstorage_column))
    
    #Pass to old for next time iteration
    total_storage_column_old = total_storage_column  


print()
#print('At each time step, the increment should be equal to the flux at the top boundary: 0.0001 (L/T)')


#cmd='rm -r testOne*'
#os.system(cmd)
