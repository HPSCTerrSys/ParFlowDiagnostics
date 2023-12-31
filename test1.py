import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics, printroot
print=printroot
import IO as io

#Run ParFlow test case
#Export ParFlow install directory; needs to be changed by the user
#e.g. export PARFLOW_DIR=/home/s.kollet/restore/migrated/Programs/parflow/
#Oddly this can not be done with os.system

# This can be done with os.environ, e.g.:
# os.environ['PARFLOW_DIR'] = '/home/b.bourgart/parflow/install'

# cmd='tclsh test1.tcl'
# os.system(cmd)
# cmd='$PARFLOW_DIR/bin/parflow testOne'
# os.system(cmd)


split=-1
name = 'testOne'
#Read static information
saturPF  = io.read_pfb(name + '.out.satur.00000.pfb', split=split)
sstorage = io.read_pfb(name + '.out.specific_storage.pfb', split=split)
permx    = io.read_pfb(name + '.out.perm_x.pfb', split=split)
permy    = io.read_pfb(name + '.out.perm_y.pfb', split=split)
permz    = io.read_pfb(name + '.out.perm_z.pfb', split=split)
mask     = io.read_pfb(name + '.out.mask.pfb', split=split)
poro     = io.read_pfb(name + '.out.porosity.pfb', split=split)
mannings = None
slopex   = None
slopey   = None

dx = dy = dz = ht.float64(1.0)
nx = 10
ny = 10
nz = 8
dzmult = ht.full(nz,1.0,dtype=ht.float64,split=None)

dt = ht.float64(1.0)
nt = 10

perm = ht.zeros((3,nz,ny,nx),split=split)
perm[0]=permz
perm[1]=permy
perm[2]=permx

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)

#Set the mask to one in active and zero in inactive regions
mask  = ht.where(mask>0.0, 1.0, mask)
#Some other constant values
ssat  = ht.full(shape3D,1.0,dtype=ht.float64,split=None)
sres  = ht.full(shape3D,0.2,dtype=ht.float64,split=None)
alpha = ht.full(shape3D,1.0,dtype=ht.float64,split=None)
nvg   = ht.full(shape3D,2.0,dtype=ht.float64,split=None)

terrainfollowing = False

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, terrainfollowing, split)
results=[]

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
      #Change in storage for cell
      cell_balance = subsurface_storage_old - subsurface_storage
      #Divergene over a cell
      cell_balance += dt * (flowleft-flowright+flowfront-flowback+flowbottom-flowtop)
      #Balance over columns
      column_balance = ht.sum(cell_balance,axis=0)
      #Mass balance over full domain without flux at the top boundary
      print('Time step:',t, ', Increment:',(ht.sum(cell_balance)/(nx*ny)).item())
      results.append(np.round(ht.sum(cell_balance).numpy()/(nx*ny), 6))

print()
ht.MPI_WORLD.Barrier()
print()
print('At each time step, the increment should be equal to the flux at the top boundary: 0.0001 (L/T)', results)


# cmd='rm -r testOne*'
# os.system(cmd)
