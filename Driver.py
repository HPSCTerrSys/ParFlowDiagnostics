import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics 
import IO as io

path = '/p/scratch/cesmtst/esmtst00/cordex3km/'
press    = io.read_pfb(path + 'cordexIHME.out.press.00000.pfb',None)
satur    = io.read_pfb(path + 'cordexIHME.out.satur.00000.pfb',None)
sstorage = io.read_pfb(path + 'cordexIHME.out.specific_storage.pfb',None)
mask     = io.read_pfb(path + 'cordexIHME.out.mask.pfb',None)
poro     = io.read_pfb(path + 'cordexIHME.out.porosity.pfb',None)
mannings = io.read_pfb(path + 'cordexIHME.out.mannings.pfb',None)
slopex   = io.read_pfb(path + 'ParFlow_MB3km_SLPX_x1592y1544.pfb',None)
slopey   = io.read_pfb(path + 'ParFlow_MB3km_SLPY_x1592y1544.pfb',None)
mannings = io.read_pfb(path + 'cordexIHME.out.mannings.pfb',None)

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
diag = Diagnostics(press, satur, mask, poro, sstorage, mannings, slopex, slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz)

#Returns an unmasked 3D field of subsurface storage, (L/L^2)
subsurface_storage=diag.SubsurfaceStorage()

#Returns an unmasked 3D field of volumetric soil moisture, (L^3/L^3)
volumetric_moisture = diag.VolumetricMoisture()

#Calculate total subsurface storage, (L)
subsurface_storage = subsurface_storage * mask
total_subsurface_storage = ht.sum(subsurface_storage)

#Generate a 2D field with the k-indices of top layer
diag.Toplayer()

#Calculate overland flow, (L^2/T)
overland_flow_x,overland_flow_y = diag.OverlandFlow()

#Calculate net lateral overland flow for each grid cell; should that be a function?
Nix = ht.zeros(shape2D, split=None)

##Calc flow east
##ParFlow:ke_[io] = pfmax(qx_[io], 0.0) - pfmax(-qx_[io + 1], 0.0);
flow_east = ht.zeros(shape2D, split=None)
for i in range (Nx-1):
    flow_east[:,i]  = ht.maximum(overland_flow_x[:,i], Nix[:,i]) 
    flow_east[:,i] -= ht.maximum((-1)*overland_flow_x[:,i+1], Nix[:,i+1])

##Calc flow west
##ParFlow:kw_[io] = pfmax(qx_[io - 1], 0.0) - pfmax(-qx_[io], 0.0);
flow_west = ht.zeros(shape2D, split=None)
for i in range (1,Nx):
    flow_west[:,i]  = ht.maximum(overland_flow_x[:,i-1], Nix[:,i-1])
    flow_west[:,i] -= ht.maximum((-1)*overland_flow_x[:,i], Nix[:,i])

##Calc flow north
##ParFlow:kn_[io] = pfmax(qy_[io], 0.0) - pfmax(-qy_[io + sy_p], 0.0);
flow_north = ht.zeros(shape2D, split=None)
for j in range (Ny-1):
    flow_north[j,:]  = ht.maximum(overland_flow_y[j,:], Nix[j,:])
    flow_north[j,:] -= ht.maximum((-1)*overland_flow_y[j+1,:], Nix[j+1,:])

##Calc flow south
##ParFlow:ks_[io] = pfmax(qy_[io - sy_p], 0.0) - pfmax(-qy_[io], 0.0);
flow_south = ht.zeros(shape2D, split=None)
for i in range (1,Ny):
    flow_south[j,:]  = ht.maximum(overland_flow_x[j-1,:], Nix[j-1,:])
    flow_south[j,:] -= ht.maximum((-1)*overland_flow_x[j,:], Nix[j,:])

##Calc net lateral overland flow for each grid cell, (L/T)
##ParFlow: ((ke_[io] - kw_[io]) / dx + (kn_[io] - ks_[io]) / dy
net_lateral_overlandflow = (flow_east - flow_west)/Dx + (flow_north - flow_south)/Dy 
