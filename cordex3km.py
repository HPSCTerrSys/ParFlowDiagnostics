import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics
import IO as io

split=None
path = '/p/scratch/cjibg31/jibg3103/parflow_3km/heat_analysis_1/output_27/'
name = 'cordexIHME'
#Read static information
sstorage = io.read_pfb(path + name + '.out.specific_storage.pfb', split=split)
permx    = io.read_pfb(path + name + '.out.perm_x.pfb', split=split)
permy    = io.read_pfb(path + name + '.out.perm_y.pfb', split=split)
permz    = io.read_pfb(path + name + '.out.perm_z.pfb', split=split)
mask     = io.read_pfb(path + name + '.out.mask.pfb', split=split)
mask     = ht.where(mask>0.0, 1.0, mask)
poro     = io.read_pfb(path + name + '.out.porosity.pfb', split=split)
mannings = io.read_pfb(path + name + '.out.mannings.pfb', split=split)
slopex   = io.read_pfb(path + 'ParFlow_MB3km_SLPX_x1592y1544.pfb', split=split)
slopey   = io.read_pfb(path + 'ParFlow_MB3km_SLPY_x1592y1544.pfb', split=split)
#alpha    = io.read_pfb(path + 'Alpha3D.pfb',split)
#alpha    = ht.where(mask>0.0,alpha,1.0)
#nvg      = io.read_pfb(path + 'Nvg3D.pfb',split)
#nvg      = ht.where(mask>0.0,nvg,2.0)
#sres     = io.read_pfb(path + 'Sres3D.pfb',split)
#sres     = ht.where(mask>0.0,sres,0.1)

dx = dy = 3000.0
dz = 2.0
dzmult = [0.01,0.015,0.025,0.035,0.065,0.10,0.15,0.25,0.35,0.50,2.0,5.0,5.0,7.50,9.0] #List of dzScales
nx = 1592
ny = 1544
nz = 15
dt = ht.float64(1.0)
nt = 10

perm = ht.zeros((3,nz,ny,nx),split=split)
perm[0]=permz
perm[1]=permy
perm[2]=permx

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)

ssat     = ht.full(shape3D,1.0,dtype=ht.float64,split=None)
terrainfollowing = True

# Generate van Genuchten fields ################################################################################################################################################
path = '/p/scratch/cjibg31/jibg3103/parflow_3km/heat_analysis_1/output_27/'
Indi3D    = io.read_pfb(path + 'ParFlow_SOIL_INDICATOR3_x1592y1544z15.pfb',None)
Alpha3D = ht.full(shape3D,2.0,split=None)
Nvg3D   = ht.full(shape3D,3.0,split=None)
Sres3D  = ht.full(shape3D,0.1,split=None)
Geom3D  = ht.zeros(shape3D,split=None)

IndicatorInput = [   1,   2,   3,   4,   5,   6,   9999,      18,      19,      20,      21,      23,      24,      25,      26,      28,      29,   30,   31,   32,   33,   40]
GeomInput      = ['F1','F2','F3','F4','F5','F6','water',    'W1',    'W2',    'W3',    'W4',    'W6',    'W7',    'W8',    'W9',   'W11',   'W12','W13','W14','W15','W16','B40']
Alpha          = [ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,    2.0,3.548134,3.467369,2.691535,0.501187,1.122018,2.089296,0.831764,1.584893,1.621810,1.513561,  2.0,  2.0,  2.0,  2.0,  2.0]
Nvg            = [ 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,    3.0,3.162278,    2.01,    2.01,    2.01,    2.01,    2.01,    2.01,    2.01,    2.01,    2.01,  3.0,  3.0,  3.0,  3.0,  3.0]
Sres           = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,    3.0,   0.076,  0.0628, 0.05037,0.074032,0.076441,0.082031,0.093361,0.084361,0.125384,0.106704,  0.1,  0.1,  0.1,  0.1,  0.1]

for k in range(len(IndicatorInput)):
    alpha    = ht.where(Indi3D == IndicatorInput[k], Alpha[k], Alpha3D)
    nvg      = ht.where(Indi3D == IndicatorInput[k], Nvg[k],   Nvg3D)
    sres     = ht.where(Indi3D == IndicatorInput[k], Sres[k],  Sres3D)
################################################################################################################################################################################

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, terrainfollowing, split)

for t in range (nt):
    #Read ParFlow pressure file
    print(path + name + '.out.press.'+('{:05d}'.format(t))+'.pfb')
    press = io.read_pfb(path + name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
    press = ht.where(mask>0.0,press,99999.0)

    #Calculate relative saturation and relative hydraulic conductivity
    satur_,krel = diag.VanGenuchten(press)

    #Work with the ParFlow output for now, ...
    satur = io.read_pfb(path + name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
    satur = ht.where(mask>0.0,press,0.0)
    #..., because something is wrong with the generated van Genuchten fields
    print(ht.sum((satur-satur_)*mask))

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

    #Calculate net overland flow for each top layer cell (L/T)
    net_overland_flow = diag.NetLateralOverlandFlow(oflowx,oflowy)

    #Column balance
    if t > 0:
      #Read CLM sink/source files (mm/s)
      qflx_tran_veg = io.read_pfb(path + name + '.out.qflx_tran_veg.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
      qflx_infl     = io.read_pfb(path + name + '.out.qflx_infl.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
      #Convert to correct units, here m^3/h
      qflx_tran_veg = qflx_tran_veg * 3600.0 * 3000.0 * 3000.0 * 0.001
      qflx_infl     = qflx_infl     * 3600.0 * 3000.0 * 3000.0 * 0.001

      #Source/sink
      sourcesink = (qflx_tran_veg + qflx_infl) * mask[nz-1,:,:]

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
