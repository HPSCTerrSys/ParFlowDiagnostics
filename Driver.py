import heat as ht
import numpy as np
import sys
import os
import Diagnostics as diag
import IO as io

path = '/p/scratch/cesmtst/esmtst00/cordex3km/'
self.press    = io.readpfb(path + '.out.press.00001.pfb')
self.mask     = io.readpfb(path + '.out.press.00001.pfb')
self.sstorage = io.readpfb(path + '.out.press.00001.pfb')
self.satur    = io.readpfb(path + '.out.press.00001.pfb')
self.poro     = io.readpfb(path + '.out.press.00001.pfb')
self.Dx = self.Dy = 3000.0
self.Dz = 1
self.Nx = 
self.Ny =
Self.Nz = 15

diag.TotalSubsurfaceStorage()

