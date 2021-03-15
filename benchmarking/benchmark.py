#TIME_INIT#
import os
import numpy as np
import heat as ht
from Diagnostics import Diagnostics, printroot
import IO

path = '/p/scratch/cesmtst/kollet1/global/'
name = 'Global1km'
split = #SPLIT#

#TIME_BEGIN#
mask = IO.read_pfb(path + 'Global1km.out.mask.pfb', split=split)
press = IO.read_pfb(path + 'Global1km.out.press.00000.pfb', split=split)
satur = IO.read_pfb(path + 'Global1km.out.satur.00000.pfb', split=split)
#TIME_END#

# Use multiple Layers
layers = 3
mask = mask * ht.ones((layers,1,1))
press = press * (1. + 0.1 * press.std() * ht.random.randn(layers,1,1)).clip(a_min=0., a_max=None)
satur = satur * (1. + 0.1 * satur.std() * ht.random.randn(layers,1,1)).clip(a_min=0., a_max=None)

nz, ny, nx = mask.shape
shape2D=(ny, nx)
shape3D=(nz, ny, nx)
printroot('shape:', shape3D)
dz, dy, dx = 0.1, 1000., 1000.
dzmult = ht.full(nz,1.0)

#Set the mask to one in active and zero in inactive regions
mask  = ht.where(mask>0.0, 1.0, mask)
#Default values
perm  = ht.full((3,*shape3D),1.0,split=split)  # perm needs to be 4D ?!
ssat  = ht.full(shape3D,1.0,split=split)
sres  = ht.full(shape3D,0.2,split=split)
alpha = ht.full(shape3D,0.005,split=split)
nvg   = ht.full(shape3D,2.0,split=split)
slopex   = IO.read_pfb(path + 'global_1km_XSLOPE_MERIT_sea_corr.pfb',split=split)
mannings = ht.ones_like(slopex)
slopey   = IO.read_pfb(path + 'global_1km_YSLOPE_MERIT_sea_corr.pfb',split=split)
#missing data
poro = ht.full(shape3D,1.0,split=split)
sstorage = ht.full(shape3D,1.0,split=split)
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, None, Split=split)

#TIME_BEGIN#
a = diag.SubsurfaceStorage(press, satur)
#TIME_END#
#TIME_BEGIN#
b = diag._SubsurfaceStorage(press, satur)
#TIME_END#
printroot('SubsurfaceStorage', ht.allclose(a, b))


#TIME_BEGIN#
a = diag.TopLayerPressure(press)
#TIME_END#
#TIME_BEGIN#
b = diag._TopLayerPressure(press)
#TIME_END#
printroot('TopLayerPressure', ht.allclose(a, b), flush=True)
toplayerpress = a


#TIME_BEGIN#
surf_storage = diag.SurfaceStorage(toplayerpress)
#TIME_END#
printroot('surf_storage', flush=True)

#TIME_BEGIN#
vol_moisture = diag.VolumetricMoisture(satur)
#TIME_END#
printroot('vol_moisture', flush=True)

#TIME_BEGIN#
overland_flow_x, overland_flow_y = diag.OverlandFlow(toplayerpress)
#TIME_END#
printroot('overland_flow', overland_flow_x.shape, overland_flow_x.split, overland_flow_y.shape, overland_flow_y.split, flush=True)

#TIME_BEGIN#
# a = diag.NetLateralOverlandFlow(overland_flow_x, overland_flow_y)
#TIME_END#
#TIME_BEGIN#
b = diag._NetLateralOverlandFlow(overland_flow_x, overland_flow_y)
#TIME_END#
# # printroot('NetLateralOverlandFlow', ht.allclose(a, b))
printroot('NetLateralOverlandFlow', flush=True)

#TIME_BEGIN#
_, krel = diag.VanGenuchten(press)
#TIME_END#
printroot('vanGenuchten', flush=True)

#TIME_BEGIN#
# a = diag.SubsurfaceFlow(press, krel)
#TIME_END#
#TIME_BEGIN#
# b = diag._SubsurfaceFlow(press, krel)
#TIME_END#
# printroot('SubsurfaceFlow', ht.allclose(a, b), flush=True)

#TIME_FINAL#
