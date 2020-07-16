import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics
import IO as io

split=None
path = './'
name = 'profileclm1'
#Read static information
sstorage = io.read_pfb(path + name + '.out.specific_storage.pfb', split=split)
mask     = io.read_pfb(path + name + '.out.mask.pfb', split=split)
mask     = ht.where(mask>0.0, 1.0, mask)
permx    = io.read_pfb(path + name + '.out.perm_x.pfb', split=split)
permx    = ht.where(mask==1.0, permx, 0.0)
permy    = io.read_pfb(path + name + '.out.perm_y.pfb', split=split)
permy    = ht.where(mask==1.0, permy, 0.0)
permz    = io.read_pfb(path + name + '.out.perm_z.pfb', split=split)
permz    = ht.where(mask==1.0, permz, 0.0)
poro     = io.read_pfb(path + name + '.out.porosity.pfb', split=split)
mannings = io.read_pfb(path + name + '.out.mannings.pfb', split=split)
slopex   = io.read_pfb(path + name + '.out.slope_x.pfb', split=split)
slopey   = io.read_pfb(path + name + '.out.slope_y.pfb', split=split)

dx = 1.0
dy = 0.2
dz = 0.01
nx = 100
ny = 5
nz = 600
dt = ht.float64(1.0)
nt = 24

perm = ht.zeros((3,nz,ny,nx),split=split)
perm[0]=permz
perm[1]=permy
perm[2]=permx

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)

dzmult = ht.full(nz,1.0,split=split)
alpha  = ht.full(shape3D,1.0,split=split)
nvg    = ht.full(shape3D,2.0,split=split)
sres   = ht.full(shape3D,0.11,split=split)
ssat   = ht.full(shape3D,1.0,dtype=ht.float64,split=split)
#Set to False in case of profileclm and to True in case of tfgprofileclm
terrainfollowing = False 

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, terrainfollowing, split)

for t in range (nt):
    #Read ParFlow pressure file
    print(path + name + '.out.press.'+('{:05d}'.format(t))+'.pfb')
    press = io.read_pfb(path + name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
    press = ht.where(mask>0.0,press,99999.0)

    #Calculate relative saturation and relative hydraulic conductivity
    dummy,krel = diag.VanGenuchten(press)

    #Work with the ParFlow output for now, ...
    satur = io.read_pfb(path + name + '.out.satur.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
    satur = ht.where(mask==1.0,satur,0.0)

    #Obtain pressure at the land surface
    top_layer_press = diag.TopLayerPressure(press)

    #Returns an unmasked 3D field of subsurface storage, (L^3)
    subsurface_storage=diag.SubsurfaceStorage(press,satur)

    #Returns an unmasked 2D field of surface storage, (L^3)
    surface_storage = diag.SurfaceStorage(top_layer_press)

    #Calculate subsurface flow in all 6 directions for each grid cell (L^3/T)
    flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag.SubsurfaceFlow(press,krel)

    #Calculate overland flow (L^2/T)
    oflowx,oflowy = diag.OverlandFlow(top_layer_press)

    #Calculate net overland flow for each top layer cell (L^3/T)
    net_overland_flow = diag.NetLateralOverlandFlow(oflowx,oflowy)

    #Column balance
    if t > 0:
      #Read source/sink values coming from CLM
      sink = io.read_pfb(path + name + '.out.et.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
      sink = ht.where(mask==1.0,sink,0.0)
      for k in range(nz):
        sink[k,:,:] *= dz * dzmult[k] 

      #Read CLM sink/source files (mm/s)
      #qflx_tran_veg = io.read_pfb(path + name + '.out.qflx_tran_veg.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
      #qflx_infl     = io.read_pfb(path + name + '.out.qflx_infl.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
      #Convert to correct units, here m^3
      #qflx_tran_veg *= (-1.0) * dt * 3600.0 * dy * dx * 0.001
      #qflx_infl     *= dt * 3600.0 * dy * dx * 0.001

      #Source/sink (L^3)
      #sourcesink_ = (qflx_tran_veg + qflx_infl) * mask[nz-1,:,:]
      sourcesink = ht.sum(sink, axis=0) * dt * dy * dx 

      #Change in subsurface storage for each cell (L^3)
      dstorage_cell = old_subsurface_storage - subsurface_storage

      #Change in surface storage for each surface cell (L^3)
      dsurface_storage_cell = old_surface_storage - surface_storage

      #Surface balanc (L^3)
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
      balance_column += sourcesink

      #Mass balance over full domain without flux at the top boundary
      print('Time step:',t, ', dstorage:',ht.sum(dstorage_column))
      print('Time step:',t, ', divergence:',ht.sum(divergence_column))
      print('Time step:',t, ', dsurface_storage:',ht.sum(dsurface_storage_cell))
      print('Time step:',t, ', netoverlandflow:',ht.sum(net_overland_flow))
      print('Time step:',t, ', surface_balance:',ht.sum(balance_surface))
      print('Time step:',t, ', source/sink:',ht.sum(sourcesink))
      print('Time step:',t, ', total balance:',ht.sum(balance_column))

    #New becomes old in the ensuing time step
    old_subsurface_storage = subsurface_storage
    old_surface_storage = surface_storage
