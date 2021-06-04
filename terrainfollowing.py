import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics, printroot
print=printroot
import IO as io

#Run ParFlow test case
#Export ParFlow install directory; needs to be changed by the user
# os.environ['PARFLOW_DIR'] = '/p/project/cslts/slts31/parflow_ben_original/run'
# cmd='tclsh terrainfollowing.tcl'
# os.system(cmd)
# cmd='$PARFLOW_DIR/bin/parflow terrainfollowing1'
# os.system(cmd)


#Define terrain following test case
name = 'terrainfollowing1'

split=-1
#Read static information
sstorage = io.read_pfb(name + '.out.specific_storage.pfb',split=split)
mask     = io.read_pfb(name + '.out.mask.pfb',split=split)
mask     = ht.where(mask>0.0,1.0,mask)
poro     = io.read_pfb(name + '.out.porosity.pfb',split=split)
mannings = io.read_pfb(name + '.out.mannings.pfb',split=split)
slopex   = io.read_pfb(name + '.out.slope_x.pfb',split=split)
slopey   = io.read_pfb(name + '.out.slope_y.pfb',split=split)
permx    = io.read_pfb(name + '.out.perm_x.pfb',split=split)
permy    = io.read_pfb(name + '.out.perm_y.pfb',split=split)
permz    = io.read_pfb(name + '.out.perm_y.pfb',split=split)

dx = dy = 1.
dz = 0.1
nx = 10
ny = 1
nz = 10
dzmult = ht.full(nz,1.0,split=None)
terrainfollowing = True

dt = 0.1
nt = 20

perm =ht.zeros((3,nz,ny,nx),split=split)
perm[0] = permz
perm[1] = permy
perm[2] = permx

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)
#Set the mask to one in active and zero in inactive regions
#mask  = ht.where(mask>0.0, 1.0, mask)
#Some other constant values
ssat  = ht.full(shape3D,1.0,split=split)
sres  = ht.full(shape3D,0.2,split=split)
alpha = ht.full(shape3D,2.0,split=split)
nvg   = ht.full(shape3D,2.0,split=split)
#The following fields are special; while they are 2D, we need to define them in 3D,
#because they are often read from pfb, which is always 3D (with 1 in z-direction)

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, terrainfollowing, Split=split)

for t in range (nt+1):
    print(name + '.out.press.'+('{:05d}'.format(t))+'.pfb')
    press = io.read_pfb(name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
    press = ht.where(mask==0.0,99999.0,press)

    print('press.shape=', press.lshape, 'press.split=',press.split, flush=True)
    #Calculate relative saturation and relative hydraulic conductivity
    satur,krel = diag.VanGenuchten(press)

    #Obtain pressure at the land surface
    top_layer_press = diag.TopLayerPressure(press)

    #Returns an unmasked 3D field of subsurface storage, (L^3)
    subsurface_storage=diag._SubsurfaceStorage(press,satur)

    #Returns an unmasked 2D field of surface storage, (L^3)
    surface_storage = diag.SurfaceStorage(top_layer_press)

    #Calculate subsurface flow in all 6 directions for each grid cell (L^3/T)
    flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag._SubsurfaceFlow(press,krel)
    print('subsurfaceflow', flush=True)

    #Calculate overland flow (L^3/T)
    oflowx,oflowy = diag.OverlandFlow(top_layer_press)
    print('overlandflow', flush=True)

    #Calculate net overland flow for each top layer cell (L/T)
    net_overland_flow = diag._NetLateralOverlandFlow(oflowx,oflowy)

    #Column balance
    if t > 0:
      #Change in subsurface storage for each cell
      dstorage_cell = old_subsurface_storage - subsurface_storage

      #Change in surface storage for each surface cell
      dsurface_storage_cell = old_surface_storage - surface_storage

      #Surface balanc
      balance_surface = dsurface_storage_cell - dt * net_overland_flow

      #Divergence of the flux for each cell
      divergence_cell = dt * (flowleft-flowright+flowfront-flowback+flowbottom-flowtop)

      #Balance for each cell
      balance_cell = dstorage_cell + divergence_cell

      #Change in storage for each column
      dstorage_column = ht.sum(dstorage_cell*mask,axis=0)

      #Divergence of the flux for each column
      divergence_column = ht.sum(divergence_cell*mask,axis=0)

      #Balance for each column
      balance_column  = ht.sum(balance_cell*mask,axis=0)
      balance_column += balance_surface

      #Discharge out of the domain at left column

      #Mass balance over full domain without flux at the top boundary
      print('Time step:',t, ', dstorage:',ht.sum(dstorage_column).item())
      print('Time step:',t, ', divergence:',ht.sum(divergence_column).item())
      print('Time step:',t, ', dsurface_storage:',ht.sum(dsurface_storage_cell).item())
      print('Time step:',t, ', netoverlandflow:',ht.sum(net_overland_flow).item())
      print('Time step:',t, ', surface_balance:',ht.sum(balance_surface).item())
      print('Time step:',t, ', total balance:',ht.sum(balance_column).item())


    #New becomes old in the ensuing time step
    old_subsurface_storage = subsurface_storage
    old_surface_storage = surface_storage


print()
print('In the first 2 hours, the balance increment due to rain is -0.5.')



#cmd='rm -r *slab1.*'
#os.system(cmd)
