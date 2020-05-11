import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics 
import IO as io

split=None
path = '/p/scratch/cesmtst/esmtst00/cordex3km/'
#Read static information
saturPF  = io.read_pfb(path + 'cordexIHME.out.satur.00000.pfb',split)
sstorage = io.read_pfb(path + 'cordexIHME.out.specific_storage.pfb',split)
mask     = io.read_pfb(path + 'cordexIHME.out.mask.pfb',split)
poro     = io.read_pfb(path + 'cordexIHME.out.porosity.pfb',split)
mannings = io.read_pfb(path + 'cordexIHME.out.mannings.pfb',split)
slopex   = io.read_pfb(path + 'ParFlow_MB3km_SLPX_x1592y1544.pfb',split)
slopey   = io.read_pfb(path + 'ParFlow_MB3km_SLPY_x1592y1544.pfb',split)
mannings = io.read_pfb(path + 'cordexIHME.out.mannings.pfb',split)

dx = dy = 3000.0
dz = 2.0
dzmult = [0.01,0.015,0.025,0.035,0.065,0.10,0.15,0.25,0.35,0.50,2.0,5.0,5.0,7.50,9.0] #List of dzScales
nx = 1592
ny = 1544
nz = 15
nt = 2

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)

#Set the mask to one in active and zero in inactive regions

mask  = ht.where(mask>0.0, 1.0, mask)
#Default values
perm  = ht.full(shape3D,1.0,split=None)
ssat  = ht.full(shape3D,1.0,split=None)
sres  = ht.full(shape3D,0.0,split=None)
alpha = ht.full(shape3D,1.0,split=None)
nvg   = ht.full(shape3D,2.0,split=None)

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, split)

for t in range (nt):
    print('cordexIHME.out.press.'+('{:05d}'.format(t))+'.pfb')
    press    = io.read_pfb(path + 'cordexIHME.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split)

    #Calculate relative saturation and relative hydraulic conductivity
    satur,krel = diag.VanGenuchten(press)

    #Returns an unmasked 3D field of subsurface storage, (L/L^2)
    subsurface_storage=diag.SubsurfaceStorage(press,satur)

    #Returns an unmasked 3D field of volumetric soil moisture, (L^3/L^3)
    volumetric_moisture = diag.VolumetricMoisture(satur)

    #Calculate total subsurface storage, (L)
    subsurface_storage = subsurface_storage * mask
    total_subsurface_storage = ht.sum(subsurface_storage)

    #Generate a vector with the k-indices of the top layer
    top_layer = diag.Toplayer()
    top_press = ht.zeros(shape2D,split=split)
    for j in range(ny):
        for i in range(nx):
            k = top_layer[j,i].astype(int)
            print(k,j,i)
            top_press[j,i] = press[k,j,i]
    #top_press[0,0:ny-1,0:nx-1]  = press[top_vector[0:ny-1,0:nx-1],0:ny-1,0:nx-1]
    #top_press = [press[k,:,:] for k in top_vector[:,:]]
    #print(top_vector.shape)
    #press[[top_vector[:,:]],:,:]

    #Calculate overland flow, (L^2/T)
    overland_flow_x,overland_flow_y = diag.OverlandFlow(press,top_layer) 

    #Calculate net lateral overland flow
    net_later_overland_flow = diag.NetLateralOverlandFlow(overland_flow_x,overland_flow_y)

    #Calculate subsurface flow in all 6 directions for each grid cell
    flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag.SubsurfaceFlow(press,krel)

    #Column balance
    if t > 0: 
      colume_balance  = ht.zero(shape2D,split)
      column_balance  = ht.sum(satur_old - satur,axis=0)
      column_balance += ht.sum(flowleft+flowright+flowfront+flowback+flowbottom+flowtop, axis=0)
      column_balance += net_lateral_overland_flow

