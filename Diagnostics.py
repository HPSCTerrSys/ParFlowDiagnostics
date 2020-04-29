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
                 Slopex, Slopey, Sourcesink, Dx, Dy, Dz, Nx, Ny, Nz):
        self.Press = Press
        self.Satur = Satur
        self.Mask = Mask
        self.Permx = Permx
        self.Permy = Permy
        self.Permz = Permz
        self.Alpha = Alpha
        self.N = N
        self.Sres = Sres
        self.Poro = Poro
        self.Sstorage = Sstorage
        self.Mannings = Mannings
        self.Slopex = Slopex
        self.Slopey = Slopey
        self.Sourcesink = Sourcesink
        self.Dx = Dx
        self.Dy = Dy
        self.Dz = Dz
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz

    def TotalSubsurfaceStorage(self):
        shape3D = (self.Nx, self.Ny, self.Nz)
        subsurface_storage = ht.zeros(shape3D, split=self.split)
        subsurface_storage = self.Satur * self.Poro * self.Dx * self.Dy * self.Dz * self.Mask
        subsurface_storage += self. Press * self.Sstorage * self.Satur * self.Dx * self.Dy * self.Dz * self.Mask
        total_subsurface_storage = np.sum(subsurface_storage)
        return {total_subsurface_storage}


if __name__ == '__main__':
    pass
