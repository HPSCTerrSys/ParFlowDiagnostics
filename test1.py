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
cmd='tclsh test1.tcl'
os.system(cmd)
cmd='$PARFLOW_DIR/bin/parflow testOne'
os.system(cmd)


split=None
name = 'testOne'
#Read static information
saturPF  = io.read_pfb(name + '.out.satur.00000.pfb',split)
sstorage = io.read_pfb(name + '.out.specific_storage.pfb',split)
mask     = io.read_pfb(name + '.out.mask.pfb',split)
poro     = io.read_pfb(name + '.out.porosity.pfb',split)
mannings = None
slopex   = None
slopey   = None

dx = dy = dz = 1.0
nx = 10
ny = 10
nz = 8
dzmult = ht.full(nz,1.0,split=None)

dt = 1.0
nt = 10

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)

#Set the mask to one in active and zero in inactive regions
mask  = ht.where(mask>0.0, 1.0, mask)
#Some other constant values
perm  = ht.full(shape3D,0.01,split=None)
ssat  = ht.full(shape3D,1.0,split=None)
sres  = ht.full(shape3D,0.2,split=None)
alpha = ht.full(shape3D,1.0,split=None)
nvg   = ht.full(shape3D,2.0,split=None)

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, split)

for t in range (nt):
    print(name + '.out.press.'+('{:05d}'.format(t))+'.pfb')
    press    = io.read_pfb(name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split=split)

    if t>0:
        subsurface_storage_old = subsurface_storage

    #Calculate relative saturation and relative hydraulic conductivity
    satur,krel = diag.VanGenuchten(press)

    #Returns an unmasked 3D field of subsurface storage, (L^3)
    subsurface_storage=diag.SubsurfaceStorage(press,satur)

    #Calculate subsurface flow in all 6 directions for each grid cell (L^3/T)
    flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag.SubsurfaceFlow(press,krel)

    #Column balance
    if t > 0: 
      column_balance  = ht.zeros(shape2D,split=split)
      #Change in storage for each column
      column_balance  = (ht.sum(subsurface_storage_old - subsurface_storage,axis=0))
      #Add divergence of the flux for each column
      column_balance += dt * ht.sum(flowleft-flowright+flowfront-flowback+flowbottom-flowtop, axis=0)
      #Mass balance over full domain without flux at the top boundary
      print('Time step:',t, ', Increment:',ht.sum(column_balance)/(nx*ny))

print()
print('At each time step, the increment should be equal to the flux at the top boundary: 0.0001 (L/T)')


cmd='rm -r testOne*'
os.system(cmd)
