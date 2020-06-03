import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics 
import IO as io

split=None
#path = '/p/scratch/cesmtst/esmtst00/cordex3km/'
path = '/home/s.kollet/restore/migrated/Programs/parflow/test/'
#name = cordexIHME
name = 'default_richards'
#Read static information
saturPF  = io.read_pfb(path + name + '.out.satur.00000.pfb',split)
sstorage = io.read_pfb(path + name + '.out.specific_storage.pfb',split)
mask     = io.read_pfb(path + name + '.out.mask.pfb',split)
poro     = io.read_pfb(path + name + '.out.porosity.pfb',split)
#mannings = io.read_pfb(path + name + '.out.mannings.pfb',split)
mannings = None
#slopex   = io.read_pfb(path + 'ParFlow_MB3km_SLPX_x1592y1544.pfb',split)
slopex   = None
#slopey   = io.read_pfb(path + 'ParFlow_MB3km_SLPY_x1592y1544.pfb',split)
slopey   = None

#dx = dy = 3000.0
dx = 8.8888888888888893
dy = 10.666666666666666
#dz = 2.0
dz = 1.0
#dzmult = [0.01,0.015,0.025,0.035,0.065,0.10,0.15,0.25,0.35,0.50,2.0,5.0,5.0,7.50,9.0] #List of dzScales
#nx = 1592
#ny = 1544
#nz = 15
nx = 10
ny = 10
nz = 8
dzmult = ht.full(nz,1.0,split=None)

dt = 0.001
nt = 6

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
shape4D=(nt, nz, ny, nx)

#Set the mask to one in active and zero in inactive regions

mask  = ht.where(mask>0.0, 1.0, mask)
#Default values
perm  = ht.full(shape3D,1.0,split=None)
ssat  = ht.full(shape3D,1.0,split=None)
sres  = ht.full(shape3D,0.2,split=None)
alpha = ht.full(shape3D,0.005,split=None)
nvg   = ht.full(shape3D,2.0,split=None)

#Initialize Diagnostics class
#Diagnostics(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Split):
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, split)

for t in range (nt):
    print(name + '.out.press.'+('{:05d}'.format(t))+'.pfb')
    press    = io.read_pfb(path + name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb',split=split)

    if t>0:
        subsurface_storage_old = subsurface_storage

    #Calculate relative saturation and relative hydraulic conductivity
    satur,krel = diag.VanGenuchten(press)

    #Returns an unmasked 3D field of subsurface storage, (L^3)
    subsurface_storage=diag.SubsurfaceStorage(press,satur)

    #Returns an unmasked 3D field of volumetric soil moisture, (L^3/L^3)
    volumetric_moisture = diag.VolumetricMoisture(satur)

    #Calculate total subsurface storage, (L^3)
    subsurface_storage = subsurface_storage * mask
    total_subsurface_storage = ht.sum(subsurface_storage)

    #Calculate overland flow, (L^2/T)
    #overland_flow_x,overland_flow_y = diag.OverlandFlow(press,top_layer) 

    #Calculate net lateral overland flow
    #net_later_overland_flow = diag.NetLateralOverlandFlow(overland_flow_x,overland_flow_y)

    #Calculate subsurface flow in all 6 directions for each grid cell (L^3/T)
    flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag.SubsurfaceFlow(press,krel)

    #Column balance
    if t > 0: 
      print(t)
      column_balance  = ht.zeros(shape2D,split=split)
      column_balance  = (ht.sum(subsurface_storage_old - subsurface_storage,axis=0))
      print('Storage:', column_balance, subsurface_storage)
      column_balance += dt * ht.sum(flowleft+flowright+flowfront+flowback+flowbottom+flowtop, axis=0)
      print('Total:',column_balance)
