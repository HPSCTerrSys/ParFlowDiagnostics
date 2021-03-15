import heat as ht
import numpy as np
import sys
import os
from Diagnostics import Diagnostics
import IO as io

split=-1
name = 'testOne'
#Read static information
saturPF  = io.read_pfb(name + '.out.satur.00000.pfb', split=split, comm=ht.MPI_SELF)
sstorage = io.read_pfb(name + '.out.specific_storage.pfb', split=split, comm=ht.MPI_SELF)
permx    = io.read_pfb(name + '.out.perm_x.pfb', split=split, comm=ht.MPI_SELF)
permy    = io.read_pfb(name + '.out.perm_y.pfb', split=split, comm=ht.MPI_SELF)
permz    = io.read_pfb(name + '.out.perm_z.pfb', split=split, comm=ht.MPI_SELF)
mask     = io.read_pfb(name + '.out.mask.pfb', split=split, comm=ht.MPI_SELF)
poro     = io.read_pfb(name + '.out.porosity.pfb', split=split, comm=ht.MPI_SELF)

saturPF = ht.array(saturPF.larray, split=split)
sstorage = ht.array(sstorage.larray, split=split)
permx = ht.array(permx.larray, split=split)
permy = ht.array(permy.larray, split=split)
permz = ht.array(permz.larray, split=split)
mask = ht.array(mask.larray, split=split)
poro = ht.array(poro.larray, split=split)


saturPF2  = io.read_pfb_mpi(name + '.out.satur.00000.pfb', split=split)
sstorage2 = io.read_pfb_mpi(name + '.out.specific_storage.pfb', split=split)
permx2    = io.read_pfb_mpi(name + '.out.perm_x.pfb', split=split)
permy2    = io.read_pfb_mpi(name + '.out.perm_y.pfb', split=split)
permz2    = io.read_pfb_mpi(name + '.out.perm_z.pfb', split=split)
mask2     = io.read_pfb_mpi(name + '.out.mask.pfb', split=split)
poro2     = io.read_pfb_mpi(name + '.out.porosity.pfb', split=split)

print(ht.allclose(saturPF, saturPF2))
print(ht.allclose(sstorage, sstorage2))
print(ht.allclose(permx, permx2))
print(ht.allclose(permy, permy2))
print(ht.allclose(permz, permz2))
print(ht.allclose(mask, mask2))
print(ht.allclose(poro, poro2))

name = 'terrainfollowing1'
#Read static information
sstorage = io.read_pfb(name + '.out.specific_storage.pfb', split=split, comm=ht.MPI_SELF)
permx    = io.read_pfb(name + '.out.perm_x.pfb', split=split, comm=ht.MPI_SELF)
permy    = io.read_pfb(name + '.out.perm_y.pfb', split=split, comm=ht.MPI_SELF)
permz    = io.read_pfb(name + '.out.perm_z.pfb', split=split, comm=ht.MPI_SELF)
mask     = io.read_pfb(name + '.out.mask.pfb', split=split, comm=ht.MPI_SELF)
poro     = io.read_pfb(name + '.out.porosity.pfb', split=split, comm=ht.MPI_SELF)

sstorage = ht.array(sstorage.larray, split=split)
permx = ht.array(permx.larray, split=split)
permy = ht.array(permy.larray, split=split)
permz = ht.array(permz.larray, split=split)
mask = ht.array(mask.larray, split=split)
poro = ht.array(poro.larray, split=split)

sstorage2 = io.read_pfb_mpi(name + '.out.specific_storage.pfb', split=split)
permx2    = io.read_pfb_mpi(name + '.out.perm_x.pfb', split=split)
permy2    = io.read_pfb_mpi(name + '.out.perm_y.pfb', split=split)
permz2    = io.read_pfb_mpi(name + '.out.perm_z.pfb', split=split)
mask2     = io.read_pfb_mpi(name + '.out.mask.pfb', split=split)
poro2     = io.read_pfb_mpi(name + '.out.porosity.pfb', split=split)

print(ht.allclose(sstorage, sstorage2))
print(ht.allclose(permx, permx2))
print(ht.allclose(permy, permy2))
print(ht.allclose(permz, permz2))
print(ht.allclose(mask, mask2))
print(ht.allclose(poro, poro2))

name = 'slab1'
#Read static information
sstorage = io.read_pfb(name + '.out.specific_storage.pfb', split=split, comm=ht.MPI_SELF)
permx    = io.read_pfb(name + '.out.perm_x.pfb', split=split, comm=ht.MPI_SELF)
permy    = io.read_pfb(name + '.out.perm_y.pfb', split=split, comm=ht.MPI_SELF)
permz    = io.read_pfb(name + '.out.perm_z.pfb', split=split, comm=ht.MPI_SELF)
mask     = io.read_pfb(name + '.out.mask.pfb', split=split, comm=ht.MPI_SELF)
poro     = io.read_pfb(name + '.out.porosity.pfb', split=split, comm=ht.MPI_SELF)

sstorage = ht.array(sstorage.larray, split=split)
permx = ht.array(permx.larray, split=split)
permy = ht.array(permy.larray, split=split)
permz = ht.array(permz.larray, split=split)
mask = ht.array(mask.larray, split=split)
poro = ht.array(poro.larray, split=split)

sstorage2 = io.read_pfb_mpi(name + '.out.specific_storage.pfb', split=split)
permx2    = io.read_pfb_mpi(name + '.out.perm_x.pfb', split=split)
permy2    = io.read_pfb_mpi(name + '.out.perm_y.pfb', split=split)
permz2    = io.read_pfb_mpi(name + '.out.perm_z.pfb', split=split)
mask2     = io.read_pfb_mpi(name + '.out.mask.pfb', split=split)
poro2     = io.read_pfb_mpi(name + '.out.porosity.pfb', split=split)

print(ht.allclose(sstorage, sstorage2))
print(ht.allclose(permx, permx2))
print(ht.allclose(permy, permy2))
print(ht.allclose(permz, permz2))
print(ht.allclose(mask, mask2))
print(ht.allclose(poro, poro2))

