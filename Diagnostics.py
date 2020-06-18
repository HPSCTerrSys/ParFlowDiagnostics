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
    def __init__(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Split):
        self.Mask       = Mask
        self.Perm       = Perm
        self.Poro       = Poro
        self.Sstorage   = Sstorage
        self.Ssat       = Ssat
        self.Sres       = Sres
        self.Nvg        = Nvg
        self.Alpha      = Alpha
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

    def SubsurfaceStorage(self,Press,Satur):
        shape3D = (self.Nz, self.Ny, self.Nx)
        subsurface_storage = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        for k in range(self.Nz):
            subsurface_storage = Satur * self.Poro *  self.Dx * self.Dy * self.Dz * self.Dzmult[k]
            #print(self.Poro)
            subsurface_storage += Press * self.Sstorage * Satur * self.Dx * self.Dy * self.Dz * self.Dzmult[k]
        return(subsurface_storage)

    def VolumetricMoisture(self,Satur):
        shape3D = (self.Nz, self.Ny, self.Nx)
        volumetric_moisture = ht.zeros(shape3D, split=self.Split)
        volumetric_moisture = Satur * self.Poro
        return(volumetric_moisture)

    def TopLayerPressure(self,Press):
        shape2D        = (self.Ny, self.Nx)
        Toplayerpress  = ht.full(shape2D,99999.0,split=self.Split)
        check     = ht.full(shape2D,-1,split=self.Split)
        for k in reversed(range(self.Nz)):
            Toplayerpress[:,:] = ht.where((self.Mask[k,:,:]>0.0) & (check<0), Press[k,:,:], Toplayerpress)
            check = ht.where(Toplayerpress!=99999.0, 0, check)
        return(Toplayerpress)

    def OverlandFlow(self,Toplayerpress):
        shape2D = (self.Ny,self.Nx)
        dirx = ht.full(shape2D,-1.0,split=self.Split)
        diry = ht.full(shape2D,-1.0,split=self.Split)

        dirx = ht.where(self.Slopex>0.0, 1.0, dirx)
        diry = ht.where(self.Slopey>0.0, 1.0, diry)

        #We need only the positive pressure values and set the rest to zero, which results in zero overland flow
        Toplayerpress = ht.where(Toplayerpress>0.0, Toplayerpress, 0.0)

        flowx=ht.zeros(shape2D,split=self.Split)
        flowy=ht.zeros(shape2D,split=self.Split)
        flowx[:,:] = dirx * (ht.absolute(self.Slopex[0,:,:]))**(1/2)/self.Mannings[0,:,:] * Toplayerpress**(5/3)
        flowy[:,:] = diry * (ht.absolute(self.Slopey[0,:,:]))**(1/2)/self.Mannings[0,:,:] * Toplayerpress**(5/3)

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

    def SubsurfaceFlow(self,Press,Krel):
        shape3D     = (self.Nz,self.Ny,self.Nx)
        Kmean=ht.full(shape3D,1.0,dtype=ht.float64,split=self.Split)
        
        #Calculate the flux across the right and left face
        grad      = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        flowright = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        flowleft  = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        for i in range(self.Nx-1):
            #Right and left face
            grad[:,:,i] = (Press[:,:,i+1] - Press[:,:,i])/self.Dx
            Kmean[:,:,i]= 2./(1.0/self.Perm[:,:,i] + 1.0/self.Perm[:,:,i+1])
            flowright[:,:,i] = ( ht.float(-1.0) *  Kmean[:,:,i] * ht.where(grad[:,:,i]>0.0,Krel[:,:,i+1],Krel[:,:,i]) * grad[:,:,i] )
            flowright[:,:,i] = flowright[:,:,i] * self.Mask[:,:,i]
            flowleft[:,:,i+1] = flowright[:,:,i]
        for k in range(self.Nz):
            flowright[k,:,:] = self.Dy * self.Dz * self.Dzmult[k] * flowright[k,:,:]
            flowleft[k,:,:]  = self.Dy * self.Dz * self.Dzmult[k] * flowleft[k,:,:]

        flowback = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        flowfront  = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        for j in range(self.Ny-1):
            #Back and front face
            grad[:,j,:] = (Press[:,j+1,:] - Press[:,j,:])/self.Dy 
            Kmean[:,j,:]= 2./(1.0/self.Perm[:,j,:] + 1.0/self.Perm[:,j+1,:])
            flowback[:,j,:] = ( ht.float64(-1.0) *  Kmean[:,j,:] * ht.where(grad[:,j,:]>0.0,Krel[:,j+1,:],Krel[:,j,:]) * grad[:,j,:] )
            flowback[:,j,:] = flowback[:,j,:] * self.Mask[:,j,:]
            flowfront[:,j+1,:] = flowback[:,j,:]
        for k in range(self.Nz):
            flowback[k,:,:] = self.Dx * self.Dz * self.Dzmult[k] * flowback[k,:,:]
            flowfront[k,:,:]  = self.Dx * self.Dz * self.Dzmult[k] * flowfront[k,:,:]

        flowtop = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        flowbottom  = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        for k in range(self.Nz-1):
            #Top and bottom face
            grad[k,:,:] = (Press[k+1,:,:] - Press[k,:,:])/(self.Dz * (self.Dzmult[k]/2.0 + self.Dzmult[k+1]/2.0)) + 1.0
            Kmean[k,:,:]= ( (self.Dz * (self.Dzmult[k]+self.Dzmult[k+1])) / 
                            (self.Dz*self.Dzmult[k]/self.Perm[k,:,:] + self.Dz*self.Dzmult[k+1]/self.Perm[k+1,:,:]) )
            flowtop[k,:,:] = ( ht.float64(-1.0) *  Kmean[k,:,:] * ht.where(grad[k,:,:]>0.0,Krel[k+1,:,:],Krel[k,:,:]) * grad[k,:,:] )
            flowtop[k,:,:] = flowtop[k,:,:] * self.Mask[k,:,:]
            flowbottom[k+1,:,:] = flowtop[k,:,:]
        flowtop = self.Dx * self.Dy * flowtop
        flowbottom = self.Dx * self.Dy * flowbottom

        return(flowleft,flowright,flowfront,flowback,flowbottom,flowtop)

    def VanGenuchten(self,Press):
        #ParFlow:
        #alpha = alphas[ir];
        #n = ns[ir];
        #m = 1.0e0 - (1.0e0 / n);
        #s_res = s_ress[ir];
        #s_dif = s_difs[ir];

        #if (ppdat[ipp] >= 0.0)
        #   psdat[ips] = s_dif + s_res;
        #   else
        #   {
        #     head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);
        #     psdat[ips] = s_dif / pow(1.0 + pow((alpha * head), n), m)
        #                  + s_res;
        m = 1.0 - 1.0/self.Nvg
        Satur = ht.where(Press<0.0, (self.Ssat - self.Sres)/((1.0+ (self.Alpha*ht.absolute(Press))**self.Nvg)**m) + self.Sres, 1.0)

        #ParFlow:
        #opahn = 1.0 + pow(alpha * head, n);
        #ahnm1 = pow(alpha * head, n - 1);
        #prdat[ipr] = pow(1.0 - ahnm1 / (pow(opahn, m)), 2)
        #             / pow(opahn, (m / 2));
        opahn = 1.0 + (self.Alpha * ht.absolute(Press))**self.Nvg
        ahnm1 = (self.Alpha * ht.absolute(Press))**(self.Nvg-1)
        Krel  = ht.where(Press<0.0, (1.0-ahnm1 / (opahn)**m)**2.0 / opahn**(m/2.0), 1.0)
        return(Satur,Krel)

if __name__ == '__main__':
    pass
