import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics 
import IO as io

path = '/p/scratch/cesmtst/esmtst00/cordex3km/cordexIHME'
press    = io.read_pfb(path + '.out.press.00000.pfb',None)
satur    = io.read_pfb(path + '.out.satur.00000.pfb',None)
sstorage = io.read_pfb(path + '.out.specific_storage.pfb',None)
mask     = io.read_pfb(path + '.out.mask.pfb',None)
poro     = io.read_pfb(path + '.out.porosity.pfb',None)
Dx = Dy = 3000.0
Dz = 2.0
#List of dzScales
Dzmult = [0.01,0.015,0.025,0.035,0.065,0.10,0.15,0.25,0.35,0.50,2.0,5.0,5.0,7.50,9.0]
Nx = 1592
Ny = 1544
Nz = 15

mask = ht.where(mask>0.0, 1.0, mask)

#Initialize Diagnostics class
diag = Diagnostics(press, satur, mask, poro, sstorage, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz)

#Returns an unmasked 3D field of subsurface storage, (L/L^2)
subsurface_storage=diag.SubsurfaceStorage()

#Returns an unmasked 3D field of volumetric soil moisture, (L^3/L^3)
volumetric_moisture = diag.VolumetricMoisture()

#Calculate total subsurface storage, (L)
subsurface_storage = subsurface_storage * mask
total_subsurface_storage = ht.sum(subsurface_storage)

#Generate a 2D field with the k-indices of top layer
diag.Toplayer()

