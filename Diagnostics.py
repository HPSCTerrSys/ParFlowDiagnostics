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
    def __init__(self, Press, Satur, Mask, Poro, Sstorage, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz):
        self.Press      = Press
        self.Satur      = Satur
        self.Mask       = Mask
        self.Poro       = Poro
        self.Sstorage   = Sstorage
        self.Mannings   = Mannings
        self.Slopex     = Slopex
        self.Slopey     = Slopey
        self.Dx         = Dx
        self.Dy         = Dy
        self.Dz         = Dz
        self.Dzmult     = Dzmult
        self.Nx         = Nx
        self.Ny         = Ny
        self.Nz         = Nz
        self.Top        = None

    def SubsurfaceStorage(self):
        shape3D = (self.Nz, self.Ny, self.Nx)
        subsurface_storage = ht.zeros(shape3D, split=None)
        for k in range(self.Nz):
            subsurface_storage = self.Satur * self.Poro *  self.Dz * self.Dzmult[k]
            subsurface_storage += self. Press * self.Sstorage * self.Satur * self.Dz * self.Dzmult[k]
        return(subsurface_storage)

    def VolumetricMoisture(self):
        shape3D = (self.Nz, self.Ny, self.Nx)
        volumetric_moisture = ht.zeros(shape3D, split=None)
        volumetric_moisture = self.Satur * self.Poro
        return(volumetric_moisture)

    def Toplayer(self):
        shape2D = (self.Ny, self.Nx)
        self.Top = ht.zeros(shape2D,split=None)
        check = ht.full(shape2D,-1.0,split=None)
        for k in reversed(range(self.Nz)):
            self.Top = ht.where((self.Mask[k,:,:]>0.0) & (check<0.0) , k, self.Top)
            check = ht.where(self.Top>0.0, 0.0, check)

    def OverlandFlow(self):
        shape2D = (self.Ny,self.Nx)
        dirx = ht.full(shape2D,-1.0,split=None)
        diry = ht.full(shape2D,-1.0,split=None)

        dirx = ht.where(self.Slopex>0.0, 1.0, dirx)
        diry = ht.where(self.Slopey>0.0, 1.0, diry)

        #We need to pick the pressure of the top cell/layer
        press = ht.zeros(shape2D,split=None)
        check = ht.full(shape2D,-1.0,split=None)
        for k in reversed(range(self.Nz)):
            press[:,:] = ht.where((self.Mask[k,:,:]>0.0) & (check<0.0) , self.Press[k,:,:], press)
            check[:,:] = ht.where(press!=0.0, 0.0, check)


        #We need only the positive pressure values and set the rest to zero, which results in zero overland flow
        press = ht.where(press>0.0, press, 0.0)

        flowx=ht.zeros(shape2D,split=None)
        flowy=ht.zeros(shape2D,split=None)
        flowx[:,:] = dirx * (ht.absolute(self.Slopex[0,:,:]))**(1/2)/self.Mannings[0,:,:] * press**(5/3)
        flowy[:,:] = diry * (ht.absolute(self.Slopey[0,:,:]))**(1/2)/self.Mannings[0,:,:] * press**(5/3)

        return(flowx, flowy)

if __name__ == '__main__':
    pass
