import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics 
import IO as io

#Run ParFlow test case
#Export ParFlow install directory; needs to be changed by the user
# os.environ['PARFLOW_DIR'] = '/p/project/cslts/slts31/parflow_ben_original/run'
# cmd='tclsh parkinglot.tcl'
# os.system(cmd)
# cmd='$PARFLOW_DIR/bin/parflow parkinglot1'
# os.system(cmd)


#Define slab test case
name = 'parkinglot1'
#name = 'drainageslab1'
#name = 'inflowslab1'

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
dz = 1.0
nx = 10
ny = 10
nz = 1
dzmult = ht.full(nz,1.0,split=split)
terrainfollowing = False

dt = 1.0
nt = 10

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

#Initialize hydrograph information
hydrograph = ht.zeros(nt+1, split = None)
#Gauge is located at pixel (0,0)
gaugex = gaugey = 0  

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, terrainfollowing, split)

for t in range (nt+1):
    print(name + '.out.press.'+('{:05d}'.format(t))+'.pfb')
    press = io.read_pfb(name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
    press = ht.where(mask==0.0,99999.0,press)

    #Calculate relative saturation and relative hydraulic conductivity
    satur,krel = diag.VanGenuchten(press)

    #Obtain pressure at the land surface
    top_layer_press = diag.TopLayerPressure(press)
    
    #Returns an unmasked 3D field of subsurface storage, (L^3)
    subsurface_storage=diag.SubsurfaceStorage(press,satur)

    #Returns an unmasked 2D field of surface storage, (L^3)
    surface_storage = diag.SurfaceStorage(top_layer_press)

    #Calculate subsurface flow in all 6 directions for each grid cell (L^3/T)
    flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag.SubsurfaceFlow(press,krel)

    #Calculate overland flow (L^3/T)
    oflowx,oflowy = diag.OverlandFlow(top_layer_press)
    print(ht.MPI_WORLD.rank, 'oflowx calculated', 'gshape:',oflowx.shape, 'lshape:',oflowx.lshape, 'split',oflowx.split, flush=True)
    
    #Extract overland flow at the gauge and calculate absolute discharge (L^3/T)
    ht.MPI_WORLD.Barrier()
    for i in range(0,10): # broadcast from rank 1 to 0
        try:
            oflowx[i,i]
        except:
            print('failed', i, flush=True)
    ht.MPI_WORLD.Barrier()
    
    ht.MPI_WORLD.Barrier()
    print(oflowx[0,0])
    ht.MPI_WORLD.Barrier()
    
    hydrograph[t] = dy * ht.abs(oflowx[gaugey,gaugex]) + dx * ht.abs(oflowy[gaugey,gaugex])
    
    #Calculate net overland flow for each top layer cell (L/T)
    net_overland_flow = diag.NetLateralOverlandFlow(oflowx,oflowy)

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
      print(balance_column.shape, balance_column.split, balance_surface.shape, balance_surface.split)
      balance_column += balance_surface 

      #Mass balance over full domain without flux at the top boundary
      print('Time step:',t, ', dstorage:',ht.sum(dstorage_column))
      print('Time step:',t, ', divergence:',ht.sum(divergence_column))
      print('Time step:',t, ', dsurface_storage:',ht.sum(dsurface_storage_cell))
      print('Time step:',t, ', netoverlandflow:',ht.sum(net_overland_flow))
      print('Time step:',t, ', surface_balance:',ht.sum(balance_surface))
      print('Time step:',t, ', total balance:',ht.sum(balance_column))
      
    
    #New becomes old in the ensuing time step
    old_subsurface_storage = subsurface_storage
    old_surface_storage = surface_storage

print()
print('Discharge at gauge:',hydrograph)
