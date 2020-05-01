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
Nx = 1592
Ny = 1544
Nz = 15

diag = Diagnostics(press, satur, mask, poro, sstorage, Dx, Dy, Dz, Nx, Ny, Nz)

total_subsurface_storage=diag.TotalSubsurfaceStorage()
print(total_subsurface_storage)

volumetric_moisture = diag.VolumetricMoisture()
