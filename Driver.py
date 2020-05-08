import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics 
import IO as io

Split=None
path = '/p/scratch/cesmtst/esmtst00/cordex3km/'
press    = io.read_pfb(path + 'cordexIHME.out.press.00000.pfb',Split)
satur    = io.read_pfb(path + 'cordexIHME.out.satur.00000.pfb',Split)
sstorage = io.read_pfb(path + 'cordexIHME.out.specific_storage.pfb',Split)
mask     = io.read_pfb(path + 'cordexIHME.out.mask.pfb',Split)
poro     = io.read_pfb(path + 'cordexIHME.out.porosity.pfb',Split)
mannings = io.read_pfb(path + 'cordexIHME.out.mannings.pfb',Split)
slopex   = io.read_pfb(path + 'ParFlow_MB3km_SLPX_x1592y1544.pfb',Split)
slopey   = io.read_pfb(path + 'ParFlow_MB3km_SLPY_x1592y1544.pfb',Split)
mannings = io.read_pfb(path + 'cordexIHME.out.mannings.pfb',Split)

Dx = Dy = 3000.0
Dz = 2.0
#List of dzScales
Dzmult = [0.01,0.015,0.025,0.035,0.065,0.10,0.15,0.25,0.35,0.50,2.0,5.0,5.0,7.50,9.0]
Nx = 1592
Ny = 1544
Nz = 15
Nt = 1

shape2D=(Ny, Nx)
shape3D=(Nz, Ny, Nx)
shape4D=(Nt, Nz, Ny, Nx)

#Set the mask to one in active and zero in inactive regions
mask = ht.where(mask>0.0, 1.0, mask)

#Initialize Diagnostics class
diag = Diagnostics(press, satur, mask, poro, sstorage, mannings, slopex, slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, None)

#Returns an unmasked 3D field of subsurface storage, (L/L^2)
subsurface_storage=diag.SubsurfaceStorage()

#Returns an unmasked 3D field of volumetric soil moisture, (L^3/L^3)
volumetric_moisture = diag.VolumetricMoisture()

#Calculate total subsurface storage, (L)
subsurface_storage = subsurface_storage * mask
total_subsurface_storage = ht.sum(subsurface_storage)

#Generate a vector with the k-indices of the top layer
top_vector=diag.Toplayer()

#Calculate overland flow, (L^2/T)
overland_flow_x,overland_flow_y = diag.OverlandFlow(top_vector)

#Calculate net lateral overland flow
net_later_overland_flow = diag.NetLateralOverlandFlow(overland_flow_x,overland_flow_y)
