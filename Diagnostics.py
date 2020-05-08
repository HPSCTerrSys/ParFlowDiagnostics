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
    def __init__(self, Press, Satur, Mask, Poro, Sstorage, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Split):
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
        self.Split      = Split

    def SubsurfaceStorage(self):
        shape3D = (self.Nz, self.Ny, self.Nx)
        subsurface_storage = ht.zeros(shape3D, split=self.Split)
        for k in range(self.Nz):
            subsurface_storage = self.Satur * self.Poro *  self.Dz * self.Dzmult[k]
            subsurface_storage += self. Press * self.Sstorage * self.Satur * self.Dz * self.Dzmult[k]
        return(subsurface_storage)

    def VolumetricMoisture(self):
        shape3D = (self.Nz, self.Ny, self.Nx)
        volumetric_moisture = ht.zeros(shape3D, split=self.Split)
        volumetric_moisture = self.Satur * self.Poro
        return(volumetric_moisture)

    def Toplayer(self):
        shape1D   = (self.Ny*self.Nx)
        Topvector = ht.zeros(shape1D,split=self.Split)
        shape2D   = (self.Ny, self.Nx)
        Toplayer  = ht.zeros(shape2D,split=self.Split)
        check     = ht.full(shape2D,-1.0,split=self.Split)
        for k in reversed(range(self.Nz)):
            Toplayer[:,:] = ht.where((self.Mask[k,:,:]>0.0) & (check<0.0) , k, Toplayer)
            check = ht.where(Toplayer>0.0, 0.0, check)
        Topvector = Toplayer.flatten()
        return(Topvector)

    def OverlandFlow(self,top_vector):
        shape2D = (self.Ny,self.Nx)
        dirx = ht.full(shape2D,-1.0,split=self.Split)
        diry = ht.full(shape2D,-1.0,split=self.Split)

        dirx = ht.where(self.Slopex>0.0, 1.0, dirx)
        diry = ht.where(self.Slopey>0.0, 1.0, diry)

        #We need to pick the pressure of the top cell/layer
        #Don't know how to do that with the top_vector therefore...
        press = ht.zeros(shape2D,split=self.Split)
        check = ht.full(shape2D,-1.0,split=self.Split)
        for k in reversed(range(self.Nz)):
            press[:,:] = ht.where((self.Mask[k,:,:]>0.0) & (check<0.0) , self.Press[k,:,:], press)
            check[:,:] = ht.where(press!=0.0, 0.0, check)

        #We need only the positive pressure values and set the rest to zero, which results in zero overland flow
        press = ht.where(press>0.0, press, 0.0)

        flowx=ht.zeros(shape2D,split=self.Split)
        flowy=ht.zeros(shape2D,split=self.Split)
        flowx[:,:] = dirx * (ht.absolute(self.Slopex[0,:,:]))**(1/2)/self.Mannings[0,:,:] * press**(5/3)
        flowy[:,:] = diry * (ht.absolute(self.Slopey[0,:,:]))**(1/2)/self.Mannings[0,:,:] * press**(5/3)

        return(flowx, flowy)

    def NetLateralOverlandFlow(self,overland_flow_x,overland_flow_y):
        shape2D = (self.Ny, self.Nx)
        Nix = ht.zeros(shape2D, split=self.Split)

        #Calc flow east
        #ParFlow:ke_[io] = pfmax(qx_[io], 0.0) - pfmax(-qx_[io + 1], 0.0);
        flow_east = ht.zeros(shape2D, split=self.Split)
        for i in range (self.Nx-1):
            flow_east[:,i]  = ht.maximum(overland_flow_x[:,i], Nix[:,i])
            flow_east[:,i] -= ht.maximum((-1)*overland_flow_x[:,i+1], Nix[:,i+1])

        #Calc flow west
        #ParFlow:kw_[io] = pfmax(qx_[io - 1], 0.0) - pfmax(-qx_[io], 0.0);
        flow_west = ht.zeros(shape2D, split=self.Split)
        for i in range (1,self.Nx):
            flow_west[:,i]  = ht.maximum(overland_flow_x[:,i-1], Nix[:,i-1])
            flow_west[:,i] -= ht.maximum((-1)*overland_flow_x[:,i], Nix[:,i])

        #Calc flow north
        #ParFlow:kn_[io] = pfmax(qy_[io], 0.0) - pfmax(-qy_[io + sy_p], 0.0);
        flow_north = ht.zeros(shape2D, split=self.Split)
        for j in range (self.Ny-1):
            flow_north[j,:]  = ht.maximum(overland_flow_y[j,:], Nix[j,:])
            flow_north[j,:] -= ht.maximum((-1)*overland_flow_y[j+1,:], Nix[j+1,:])

        #Calc flow south
        #ParFlow:ks_[io] = pfmax(qy_[io - sy_p], 0.0) - pfmax(-qy_[io], 0.0);
        flow_south = ht.zeros(shape2D, split=self.Split)
        for i in range (1,self.Ny):
            flow_south[j,:]  = ht.maximum(overland_flow_x[j-1,:], Nix[j-1,:])
            flow_south[j,:] -= ht.maximum((-1)*overland_flow_x[j,:], Nix[j,:])

        #Calc net lateral overland flow for each grid cell, (L/T)
        #ParFlow: ((ke_[io] - kw_[io]) / dx + (kn_[io] - ks_[io]) / dy
        net_lateral_overlandflow = (flow_east - flow_west)/self.Dx + (flow_north - flow_south)/self.Dy

        return(net_lateral_overlandflow)

if __name__ == '__main__':
    pass
