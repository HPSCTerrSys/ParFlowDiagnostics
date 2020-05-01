# <editor-fold imports
import heat as ht
import numpy as np
import sys
import os
#from . import IO

def printroot(*args, **kwargs):
    """ Only Root Process prints """
    if ht.MPI_WORLD.rank == 0:
        print(*args, **kwargs)
# </editor-fold>


class Diagnostics:  # Make this a subclass of ht.DNDarray?
    def __init__(self, Press, Satur, Mask, Poro, Sstorage, Dx, Dy, Dz, Nx, Ny, Nz):
        self.Press      = Press
        self.Satur      = Satur
        self.Mask       = Mask
        self.Poro       = Poro
        self.Sstorage   = Sstorage
        self.Dx         = Dx
        self.Dy         = Dy
        self.Dz         = Dz
        self.Nx         = Nx
        self.Ny         = Ny
        self.Nz         = Nz

    def TotalSubsurfaceStorage(self):
        shape3D = (self.Nx, self.Ny, self.Nz)
        subsurface_storage = ht.zeros(shape3D, split=None)
        subsurface_storage = self.Satur * self.Poro * self.Dx * self.Dy * self.Dz 
        subsurface_storage += self. Press * self.Sstorage * self.Satur * self.Dx * self.Dy * self.Dz 
        total_subsurface_storage = ht.sum(subsurface_storage)
        return(total_subsurface_storage)

    def VolumetricMoisture(self):
        shape3D = (self.Nx, self.Ny, self.Nz)
        volumetric_moisture = ht.zeros(shape3D, split=None)
        volumetric_moisture = self.Satur * self.Poro 
        return(volumetric_moisture)

if __name__ == '__main__':
    pass
