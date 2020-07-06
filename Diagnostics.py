# <editor-fold imports
import heat as ht
import numpy as np
import sys
import os
import math

def printroot(*args, **kwargs):
    """ Only Root Process prints """
    if ht.MPI_WORLD.rank == 0:
        print(*args, **kwargs)
# </editor-fold>


class Diagnostics:  # Make this a subclass of ht.DNDarray?
    def __init__(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split):
        self.Mask             = Mask
        self.Perm             = Perm
        self.Poro             = Poro
        self.Sstorage         = Sstorage
        self.Ssat             = Ssat
        self.Sres             = Sres
        self.Nvg              = Nvg
        self.Alpha            = Alpha
        self.Mannings         = Mannings
        self.Slopex           = Slopex
        self.Slopey           = Slopey
        self.Dx               = Dx
        self.Dy               = Dy
        self.Dz               = Dz
        self.Dzmult           = Dzmult
        self.Nx               = Nx
        self.Ny               = Ny
        self.Nz               = Nz
        self.Terrainfollowing = Terrainfollowing
        self.Split            = Split

    def SubsurfaceStorage(self, Press, Satur):
        shape3D = (self.Nz, self.Ny, self.Nx)
        subsurface_storage = ht.zeros(shape3D, dtype=ht.float64, split=self.Split)
        for k in range(self.Nz):
            subsurface_storage[k,:,:] = Satur[k,:,:] * self.Poro[k,:,:] * self.Dx * self.Dy * self.Dz * self.Dzmult[k] 
            subsurface_storage[k,:,:] += Press[k,:,:] * self.Sstorage[k,:,:] * Satur[k,:,:] * self.Dx * self.Dy * self.Dz * self.Dzmult[k]
        return(subsurface_storage)

    def VolumetricMoisture(self, Satur):
        volumetric_moisture = Satur * self.Poro
        return(volumetric_moisture)

    def TopLayerPressure(self, Press):
        shape2D = (self.Ny, self.Nx)
        fill_val = 99999.0
        Toplayerpress = ht.full(shape2D, fill_val, split=self.Split)
        check = ht.full(shape2D, -1, split=self.Split)
        for k in reversed(range(self.Nz)):
            Toplayerpress[:,:] = ht.where((self.Mask[k,:,:]>0.0) & (check[:,:]<0), Press[k,:,:], Toplayerpress[:,:])
            #Check also contains the the layer index k of the top layer
            check[:,:] = ht.where((self.Mask[k,:,:]>0.0) & (check[:,:]<0), k, check[:,:])
        return(Toplayerpress)

    def SurfaceStorage(self,Toplayerpress):
        shape2D = (self.Ny,self.Nx)
        Surfacestorage = self.Dx * self.Dy * ht.where(Toplayerpress>0.0,Toplayerpress,0.0)
        return(Surfacestorage)

    def OverlandFlow(self,Toplayerpress):
        shape2D = (self.Ny, self.Nx)
        flowx= ht.zeros(shape2D, split=self.Split)
        flowy= ht.zeros(shape2D, split=self.Split)
        dirx = ht.where(self.Slopex > 0.0, -1.0, 1.0)  
        diry = ht.where(self.Slopey > 0.0, -1.0, 1.0)

        Ponding = ht.where(Toplayerpress>0,Toplayerpress,0.0)
        #We need only the positive pressure values and set the rest to zero, which results in zero overland flow
        flowx[:,:] = self.Dy * dirx * (ht.absolute(self.Slopex[0,:,:]))**(1./2.)/self.Mannings[0,:,:] * Ponding**(5./3.)
        flowy[:,:] = self.Dx * diry * (ht.absolute(self.Slopey[0,:,:]))**(1./2.)/self.Mannings[0,:,:] * Ponding**(5./3.)
        return(flowx, flowy)

    def NetLateralOverlandFlow(self, overland_flow_x, overland_flow_y):
        shape2D = (self.Ny, self.Nx)
        Nix = ht.zeros(shape2D, split=self.Split)

        #Calc flow east
        #ParFlow:ke_[io] = pfmax(qx_[io], 0.0) - pfmax(-qx_[io + 1], 0.0);
        flow_east = ht.zeros(shape2D, split=self.Split)
        flow_east[:,self.Nx-1] = ht.where(overland_flow_x[:,self.Nx-1]>0.0,overland_flow_x[:,self.Nx-1],Nix[:,self.Nx-1])
        for i in range (self.Nx-1):
            flow_east[:,i]  = ht.maximum(overland_flow_x[:,i], Nix[:,i])
            flow_east[:,i] -= ht.maximum((-1.0)*overland_flow_x[:,i+1], Nix[:,i+1])

        #Calc flow west
        #ParFlow:kw_[io] = pfmax(qx_[io - 1], 0.0) - pfmax(-qx_[io], 0.0);
        flow_west = ht.zeros(shape2D, split=self.Split)
        flow_west[:,0] = ht.where(overland_flow_x[:,0]<0.0, overland_flow_x[:,0],Nix[:,0]) 
        for i in range (1,self.Nx):
            flow_west[:,i]  = ht.maximum(overland_flow_x[:,i-1], Nix[:,i-1])
            flow_west[:,i] -= ht.maximum((-1.0)*overland_flow_x[:,i], Nix[:,i])

        #Calc flow north
        #ParFlow:kn_[io] = pfmax(qy_[io], 0.0) - pfmax(-qy_[io + sy_p], 0.0);
        flow_north = ht.zeros(shape2D, split=self.Split)
        flow_north[self.Ny-1,:] = ht.where(overland_flow_x[self.Ny-1,:]>0.0, overland_flow_y[self.Ny-1,:],Nix[self.Ny-1,:])
        for j in range (self.Ny-1):
            flow_north[j,:]  = ht.maximum(overland_flow_y[j,:], Nix[j,:])
            flow_north[j,:] -= ht.maximum((-1.0)*overland_flow_y[j+1,:], Nix[j+1,:])

        #Calc flow south
        #ParFlow:ks_[io] = pfmax(qy_[io - sy_p], 0.0) - pfmax(-qy_[io], 0.0);
        flow_south = ht.zeros(shape2D, split=self.Split)
        flow_south[0,:] = ht.where(overland_flow_y[0,:]<0, overland_flow_x[0,:], Nix[0,:])
        for i in range (1,self.Ny):
            flow_south[j,:]  = ht.maximum(overland_flow_x[j-1,:], Nix[j-1,:])
            flow_south[j,:] -= ht.maximum((-1.0)*overland_flow_x[j,:], Nix[j,:])

        #Calc net lateral overland flow for each grid cell, (L/T)
        #ParFlow: ((ke_[io] - kw_[io]) / dx + (kn_[io] - ks_[io]) / dy
        net_lateral_overlandflow = (flow_east - flow_west)/self.Dx + (flow_north - flow_south)/self.Dy

        return(net_lateral_overlandflow)

    def SubsurfaceFlow(self, Press, Krel):
        #ParFlow
        #x_dir_g = Mean(gravity * sin(atan(x_ssl_dat[io])), gravity * sin(atan(x_ssl_dat[io + 1])));
        #x_dir_g_c = Mean(gravity * cos(atan(x_ssl_dat[io])), gravity * cos(atan(x_ssl_dat[io + 1])));
        #y_dir_g = Mean(gravity * sin(atan(y_ssl_dat[io])), gravity * sin(atan(y_ssl_dat[io + sy_p])));
        #y_dir_g_c = Mean(gravity * cos(atan(y_ssl_dat[io])), gravity * cos(atan(y_ssl_dat[io + sy_p])));
        shape3D   = (self.Nz,self.Ny,self.Nx)
        x_dir_g   = ht.zeros(shape3D, split=self.Split)
        x_dir_g_c = ht.zeros(shape3D, split=self.Split)
        y_dir_g   = ht.zeros(shape3D, split=self.Split)
        y_dir_g_c = ht.zeros(shape3D, split=self.Split)
        if self.Terrainfollowing:
            print('Terrain following')
            for k in range (self.Nz):
                #Expression for sin(atan(x)), because atan(x) does not exist in heat
                x_dir_g[k,:,:-1]   = ( self.Slopex[0,:,:-1]/ht.sqrt(self.Slopex[0,:,:-1]**2.0+1.0) + self.Slopex[0,:,1:]/ht.sqrt(self.Slopex[0,:,1:]**2.0+1.0) )/2.0
                #Expression for cos(atan(x)), because atan(x) does not exist in heat
                x_dir_g_c[k,:,:-1] = ( 1./ht.sqrt(self.Slopex[0,:,:-1] + 1.0) + 1./ht.sqrt(self.Slopex[:,:,1:] + 1.0) )/2.0 
                #Expression for sin(atan(x)), because atan(x) does not exist in heat
                y_dir_g[k,:-1,:]   = ( self.Slopey[0,:-1,:]/ht.sqrt(self.Slopey[0,:-1,:]**2.0+1.0) + self.Slopey[0,1:,:]/ht.sqrt(self.Slopey[0,1:,:]**2.0+1.0) )/2.0
                #Expression for cos(atan(x)), because atan(x) does not exist in heat
                y_dir_g_c[k,:-1,:] = ( 1./ht.sqrt(self.Slopey[0,:-1,:] + 1.0) + 1./ht.sqrt(self.Slopey[:,1:,:] + 1.0) )/2.0 

        Dzmult3D = ht.array(self.Dzmult, dtype=ht.float64).expand_dims(axis=-1).expand_dims(axis=-1)
        inv_perm = 1.0 / self.Perm

        #Calculate the flux across the right and left face
        # Left and Right
        flowright = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        flowleft  = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)

        #Note, in the inactive cells, Perm is zero, thus, 1/Perm results in inf, which then results in Kmean = 0!
        Kmean = 2. / (inv_perm[2,:, :, :-1] + inv_perm[2,:, :, 1:])
        #diff = pp[ip] - pp[ip + 1];
        #updir = (diff / dx) * x_dir_g_c - x_dir_g;
        grad = ht.diff(Press, axis=2)/self.Dx
        grad = grad * x_dir_g_c[:,:,:-1] - x_dir_g[:,:,:-1]

        flowright[:, :, :-1] = -1. * Kmean * grad * ht.where(grad > 0.0, Krel[:, :, 1:], Krel[:, :, :-1]) 
        flowleft[:,:,1:] = flowright[:,:,:-1]

        flowright *= self.Dy * self.Dz * Dzmult3D
        flowleft  *= self.Dy * self.Dz * Dzmult3D  # save this by setting flowleft after multiplication

        # Front and Back
        flowback = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        flowfront = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)

        #Note, in the inactive cells, Perm is zero, thus, 1/Perm results in inf, which then results in Kmean = 0!
        Kmean = 2. / (inv_perm[1,:, :-1, :] + inv_perm[1,:, 1:, :])
        grad = ht.diff(Press, axis=1)/self.Dy 
        grad = grad * y_dir_g_c[:,:-1,:] - y_dir_g[:,:-1,:]

        flowback[:, :-1, :] = -1. * Kmean * grad * ht.where(grad > 0.0, Krel[:, 1:, :], Krel[:, :-1, :]) 
        flowfront[:, 1:, :] = flowback[:, :-1, :]

        flowback *= self.Dx * self.Dz * Dzmult3D
        flowfront  *= self.Dx * self.Dz * Dzmult3D

        #  Top and Bottom
        flowtop = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)
        flowbottom = ht.zeros(shape3D,dtype=ht.float64,split=self.Split)

        #Note, in the inactive cells, Perm is zero, thus, 1/Perm results in inf, which then results in Kmean = 0!
        Kmean = ( (Dzmult3D[:-1] + Dzmult3D[1:]) /
                (Dzmult3D[:-1]/self.Perm[0,:-1,:,:] + Dzmult3D[1:]/self.Perm[0,1:,:,:]) )

        grad = ht.float64(1.) + ht.diff(Press, axis=0) * 2. / (self.Dz * (Dzmult3D[:-1] + Dzmult3D[1:]))

        #Application of mask checks if node k is active
        flowtop[:-1, :, :] = ht.float64(-1.) * Kmean * grad * ht.where(grad > 0.0, Krel[1:, :, :], Krel[:-1, :, :]) 
        flowbottom[1:, :, :] = flowtop[:-1, :, :]

        flowtop *= self.Dx * self.Dy
        flowbottom *= self.Dx * self.Dy

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
        m = ht.float64(1.0) - ht.float64(1.0)/self.Nvg
        Satur = ht.where(Press<0.0, (self.Ssat - self.Sres)/((ht.float64(1.0)+ (self.Alpha*ht.absolute(Press))**self.Nvg)**m) + self.Sres, ht.float64(1.0))

        #ParFlow:
        #opahn = 1.0 + pow(alpha * head, n);
        #ahnm1 = pow(alpha * head, n - 1);
        #prdat[ipr] = pow(1.0 - ahnm1 / (pow(opahn, m)), 2)
        #             / pow(opahn, (m / 2));
        opahn = ht.float64(1.0) + (self.Alpha * ht.abs(Press))**self.Nvg
        ahnm1 = (self.Alpha * ht.abs(Press))**(self.Nvg-ht.float64(1.))
        Krel  = ht.where(Press<ht.float64(0.0), (ht.float64(1.0)-ahnm1 / (opahn)**m)**ht.float64(2.0) / opahn**(m/ht.float64(2.0)), ht.float64(1.0))
        return(Satur,Krel)

if __name__ == '__main__':
    pass
