import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics 
import IO as io

path = '/p/scratch/cesmtst/esmtst00/cordex3km/'
Indi3D    = io.read_pfb(path + 'ParFlow_SOIL_INDICATOR3_x1592y1544z15.pfb',None)
Dx = Dy = 3000.0
Dz = 2.0
#List of dzScales
Dzmult = [0.01,0.015,0.025,0.035,0.065,0.10,0.15,0.25,0.35,0.50,2.0,5.0,5.0,7.50,9.0]
Nx = 1592
Ny = 1544
Nz = 15
shape3D = (Ny,Nx)
Alpha3D = ht.full(shape3D,2.0,split=None)
Nvg3D   = ht.full(shape3D,3.0,split=None)
Sres3D  = ht.full(shape3D,0.1,split=None)
Geom3D  = ht.zeros(shape3D,split=None)

IndicatorInput = [ 1,   2,  3,  4,  5,  6,  9999,      18,      19,      20,      21,      23,      24,      25,      26,      28,      29, 30, 31, 32, 33, 40]
GeomInput      = [ F1, F2, F3, F4, F5, F6, water,      W1,      W2,      W3,      W4,      W6,      W7,      W8,      W9,     W11,     W12,W13,W14,W15,W16,B40]
Alpha          = [2.0,2.0,2.0,2.0,2.0,2.0,   2.0,3.548134,3.467369,2.691535,0.501187,1.122018,2.089296,0.831764,1.584893,1.621810,1.513561,2.0,2.0,2.0,2.0,2.0]
Nvg            = [3.0,3.0,3.0,3.0,3.0,3.0,   3.0,3.162278,    2.01,    2.01,    2.01,    2.01,    2.01,    2.01,    2.01,    2.01,    2.01,3.0,3.0,3.0,3.0,3.0]
Sres           = [0.1,0.1,0.1,0.1,0.1,0.1,   3.0,   0.076,  0.0628, 0.05037,0.074032,0.076441,0.082031,0.093361,0.084361,0.125384,0.106704,0.1,0.1,0.1,0.1,0.1]

for k in range(len(IndicatorInput)):
    Alpha3D  = ht.where(Indi3D == IndicatorInput[k], Alpha[k], Alpha3D)
    Nvg3D    = ht.where(Indi3D == IndicatorInput[k], Nvg[k],   Nvg3D)
    Sres3D   = ht.where(Indi3D == IndicatorInput[k], Sres[k],  Sres3D)


io.write_pfb(path + 'Alpha3D.pfb',None)
io.write_pfb(path + 'Nvg3D.pfb',None)
io.write_pfb(path + 'Sres3D.pfb',None)
