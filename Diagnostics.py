# <editor-fold imports
import heat as ht
import numpy as np
import sys
import os
from . import IO

def printroot(*args, **kwargs):
    """ Only Root Process prints """
    if ht.MPI_WORLD.rank == 0:
        print(*args, **kwargs)
# </editor-fold>


class Diagnostics():  # Make this a subclass of ht.DNDarray?
    def __init__(self, Press, Satur, Mask, Permx, Permy, Permz,
                 Alpha, N, Sres, Poro, Sstorage, Mannings,
                 Slopex, Slopey, Sourcesink, Dx, Dy, Dz, Nx, Ny, Nz, split=None):
        self.Press      = ht.array(Press,      split=split)
        self.Satur      = ht.array(Satur,      split=split)
        self.Mask       = ht.array(Mask,       split=split)
        self.Permx      = ht.array(Permx,      split=split)
        self.Permy      = ht.array(Permy,      split=split)
        self.Permz      = ht.array(Permz,      split=split)
        self.Alpha      = ht.array(Alpha,      split=split)
        self.N          = ht.array(N,          split=split)
        self.Sres       = ht.array(Sres,       split=split)
        self.Poro       = ht.array(Poro,       split=split)
        self.Sstorage   = ht.array(Sstorage,   split=split)
        self.Mannings   = ht.array(Mannings,   split=split)
        self.Slopex     = ht.array(Slopex,     split=split)
        self.Slopey     = ht.array(Slopey,     split=split)
        self.Sourcesink = ht.array(Sourcesink, split=split)
        self.Dx         = float(Dx)
        self.Dy         = float(Dy)
        self.Dz         = float(Dz)
        self.Nx         = int(Nx)
        self.Ny         = int(Ny)
        self.Nz         = int(Nz)
        self.split      = split

    def TotalSubsurfaceStorage(self):
        shape3D = (self.Nx, self.Ny, self.Nz)
        subsurface_storage = ht.zeros(shape3D, split=self.split)
        subsurface_storage = self.Satur * self.Poro * self.Dx * self.Dy * self.Dz * self.Mask
        subsurface_storage += self. Press * self.Sstorage * self.Satur * self.Dx * self.Dy * self.Dz * self.Mask
        total_subsurface_storage = np.sum(subsurface_storage)
        return {total_subsurface_storage}


if __name__ == '__main__':
    pass
