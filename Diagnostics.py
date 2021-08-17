# <editor-fold imports
import heat as ht
import numpy as np
import sys
import os
import math

def printroot(*args, **kwargs):
    """ Only Root Process prints """
    kwargs['flush'] = True
    if ht.MPI_WORLD.rank == 0:
        print(*args, **kwargs)
# </editor-fold>


class Diagnostics:  # Make this a subclass of ht.DNDarray?
    def __init__(self, Mask, Perm, Poro, Sstorage, Ssat, Sres, Nvg, Alpha, Mannings, Slopex, Slopey, Dx, Dy, Dz, Dzmult, Nx, Ny, Nz, Terrainfollowing, Split=np.nan, Split3D=np.nan, Split2D=np.nan):
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
        if Split is None or not np.isnan(Split):
            self.Split3D = Split
            self.Split2D = Split
            self.Split = Split
        elif Split3D is None or not np.isnan(Split3D):
            self.Split3D = Split3D
            if Split2D is None or not np.isnan(Split2D):
                self.Split2D = Split2D
            else:
                self.Split2D = self.Split3D
                self.Split = None
        else:
            self.Split3D = None
            self.Split2D = None
            self.Split = None
        #self.Split3D          = Split3D if Split is None else Split
        #self.Split2D          = Split2D if Split2D is not None else self.Split3D

    def SubsurfaceStorage(self, Press, Satur):
        shape3D = (self.Nz, self.Ny, self.Nx)
        subsurface_storage = ht.zeros(shape3D, dtype=ht.float64, split=self.Split3D)
        for k in range(self.Nz):
            subsurface_storage[k,:,:] = Satur[k,:,:] * self.Poro[k,:,:] * self.Dx * self.Dy * self.Dz * self.Dzmult[k]
            subsurface_storage[k,:,:] += Press[k,:,:] * self.Sstorage[k,:,:] * Satur[k,:,:] * self.Dx * self.Dy * self.Dz * self.Dzmult[k]
        return(subsurface_storage)

    def VolumetricMoisture(self, Satur):
        volumetric_moisture = ht.mul(Satur,self.Poro)
        return(volumetric_moisture)

    def TopLayerPressure(self, Press):
        shape2D = (self.Ny, self.Nx)
        fill_val = 99999.0
        #Toplayerpress = ht.full(shape2D, fill_val, dtype=ht.float64, split=self.Split2D)
#         Toplayerpress = ht.full_like(Press[0], fill_val)
        check = ht.full(shape2D, -1, split=self.Mask[0].split)
        Toplayerpress = ht.full_like(check, fill_val, dtype=ht.float64)
        for k in reversed(range(self.Nz)):
            mask = self.Mask[k]      #ht.asarray(self.Mask[k].larray, is_split=self.Split2D)
            press = Press[k]         #ht.asarray(Press[k].larray, is_split=self.Split2D)
            cond1 = (mask>0).astype(ht.bool)
            cond2 = (check<0).astype(ht.bool)
            #print(cond1.lnumel, cond1.dtype, cond2.lnumel, cond2.dtype, flush=True)
            Toplayerpress[:,:] = ht.where(cond1 & cond2, press, Toplayerpress[:,:])
            #Check also contains the the layer index k of the top layer
            check[:,:] = ht.where((mask>0.0) & (check[:,:]<0), k, check[:,:])
        return(Toplayerpress)

    def SurfaceStorage(self,Toplayerpress):
        shape2D = (self.Ny,self.Nx)
        Surfacestorage = self.Dx * self.Dy * ht.where(Toplayerpress>0.0,Toplayerpress,0.0)
        return(Surfacestorage)

    def OverlandFlow(self,Toplayerpress):
        shape2D = (self.Ny, self.Nx)
        flowx= ht.zeros(shape2D, split=self.Split2D)
        flowy= ht.zeros(shape2D, split=self.Split2D)
        dirx = ht.where(self.Slopex[0] > 0.0, -1.0, 1.0)
        diry = ht.where(self.Slopey[0] > 0.0, -1.0, 1.0)

        Ponding = ht.where(Toplayerpress>0,Toplayerpress,0.0)
        #We need only the positive pressure values and set the rest to zero, which results in zero overland flow
        #x_slope_mannings = ht.asarray(((ht.absolute(self.Slopex[0,:,:]))**(1./2.)/self.Mannings[0,:,:]).larray, is_split=self.Split2D)
        x_slope_mannings = (ht.absolute(self.Slopex[0,:,:]))**(1./2.)/self.Mannings[0,:,:]
        #y_slope_mannings = ht.asarray(((ht.absolute(self.Slopey[0,:,:]))**(1./2.)/self.Mannings[0,:,:]).larray, is_split=self.Split2D)
        y_slope_mannings = (ht.absolute(self.Slopey[0,:,:]))**(1./2.)/self.Mannings[0,:,:]
        flowx[:,:] = dirx * x_slope_mannings * Ponding**(5./3.)
        flowy[:,:] = diry * y_slope_mannings * Ponding**(5./3.)
        return(flowx, flowy)

    def NetLateralOverlandFlow(self, overland_flow_x, overland_flow_y):
        shape2D = (self.Ny, self.Nx)
        #Nix = ht.zeros(shape2D, split=self.Split2D)

        #Calc flow east
        #ParFlow:ke_[io] = pfmax(qx_[io], 0.0) - pfmax(-qx_[io + 1], 0.0);
        flow_east = ht.zeros(shape2D, split=0)
        overland_flow_x = ht.resplit(overland_flow_x, 0) if overland_flow_x.split is not None and overland_flow_x.split != 0 else overland_flow_x

        flow_east[:,self.Nx-1] = ht.where(overland_flow_x[:,self.Nx-1]>0.0, overland_flow_x[:,self.Nx-1], 0.0)
        for i in range (self.Nx-1):
            #printroot('.', end='')
            flow_east[:,i]  = ht.maximum(overland_flow_x[:,i], 0.0)
            flow_east[:,i] -= ht.maximum((-1.0)*overland_flow_x[:,i+1], 0.0)
        flow_east.resplit_(self.Split2D)
        printroot('flow east')

        #Calc flow west
        #ParFlow:kw_[io] = pfmax(qx_[io - 1], 0.0) - pfmax(-qx_[io], 0.0);
        flow_west = ht.zeros(shape2D, split=0)
        flow_west[:,0] = ht.where(overland_flow_x[:,0]<0.0, overland_flow_x[:,0], 0.0)
        for i in range (1,self.Nx):
            #printroot('.', end='')
            flow_west[:,i]  = ht.maximum(overland_flow_x[:,i-1], 0.0)
            flow_west[:,i] -= ht.maximum((-1.0)*overland_flow_x[:,i], 0.0)
        flow_west.resplit_(self.Split2D)
        printroot('flow west')

        #Calc flow north
        #ParFlow:kn_[io] = pfmax(qy_[io], 0.0) - pfmax(-qy_[io + sy_p], 0.0);
        overland_flow_y = ht.resplit(overland_flow_y, 1) if overland_flow_y.split is not None and overland_flow_y.split != 1 else overland_flow_y
        flow_north = ht.zeros(shape2D, split=1)
        flow_north[self.Ny-1,:] = ht.where(overland_flow_y[self.Ny-1,:]>0.0, overland_flow_y[self.Ny-1,:],0.0)
        for j in range(self.Ny-1):
            #printroot('.', end='')
            flow_north[j,:]  = ht.maximum(overland_flow_y[j,:], 0.0)
            flow_north[j,:] -= ht.maximum((-1.0)*overland_flow_y[j+1,:], 0.0)
        flow_north.resplit_(self.Split2D)
        printroot('flow_north')

        #Calc flow south
        #ParFlow:ks_[io] = pfmax(qy_[io - sy_p], 0.0) - pfmax(-qy_[io], 0.0);
        flow_south = ht.zeros(shape2D, split=1)
        flow_south[0,:] = ht.where(overland_flow_y[0,:]<0.0, overland_flow_y[0,:], 0.0)
        for j in range (1,self.Ny):
            #printroot('.', end='')
            flow_south[j,:]  = ht.maximum(overland_flow_y[j-1,:], 0.0)
            flow_south[j,:] -= ht.maximum((-1.0)*overland_flow_y[j,:], 0.0)
        flow_south.resplit_(self.Split2D)
        printroot('flow south')

        #Calc net lateral overland flow for each grid cell, (L^3/T)
        #ParFlow: ((ke_[io] - kw_[io]) / dx + (kn_[io] - ks_[io]) / dy
        net_lateral_overlandflow = ht.zeros(shape2D, split=self.Split2D)
        net_lateral_overlandflow = self.Dy * (flow_east - flow_west) + self.Dx * (flow_north - flow_south)

        return(net_lateral_overlandflow)

    def SubsurfaceFlow(self, Press, Krel):
        #ParFlow
        #x_dir_g = Mean(gravity * sin(atan(x_ssl_dat[io])), gravity * sin(atan(x_ssl_dat[io + 1])));
        #x_dir_g_c = Mean(gravity * cos(atan(x_ssl_dat[io])), gravity * cos(atan(x_ssl_dat[io + 1])));
        #y_dir_g = Mean(gravity * sin(atan(y_ssl_dat[io])), gravity * sin(atan(y_ssl_dat[io + sy_p])));
        #y_dir_g_c = Mean(gravity * cos(atan(y_ssl_dat[io])), gravity * cos(atan(y_ssl_dat[io + sy_p])));
        shape3D   = (self.Nz,self.Ny,self.Nx)
        x_dir_g   = ht.zeros(shape3D, split=self.Split3D)
        x_dir_g_c = ht.full(shape3D, 1.0, split=self.Split3D)
        y_dir_g   = ht.zeros(shape3D, split=self.Split3D)
        y_dir_g_c = ht.full(shape3D, 1.0, split=self.Split3D)

        if self.Terrainfollowing:
            printroot('Terrain following')
            for k in range (self.Nz):
                x_dir_g[k,:,:-1]   = ( ht.arctan(self.Slopex[0,:,:-1]) + ht.arctan(self.Slopex[0,:,1:]) )/2.0
                x_dir_g_c[k,:,:-1] = ( ht.arctan(self.Slopex[0,:,:-1]) + ht.arctan(self.Slopex[0,:,1:]) )/2.0
                y_dir_g[k,:-1,:]   = ( ht.arctan(self.Slopey[0,:-1,:]) + ht.arctan(self.Slopey[0,1:,:]) )/2.0
                y_dir_g_c[k,:-1,:] = ( ht.arctan(self.Slopey[0,:-1,:]) + ht.arctan(self.Slopey[0,1:,:]) )/2.0

        Dzmult3D = ht.array(self.Dzmult, dtype=ht.float64).expand_dims(axis=-1).expand_dims(axis=-1)
        inv_perm = 1.0 / self.Perm

        #Calculate the flux across the right and left face
        # Left and Right
        flowright = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)
        flowleft  = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)

        #Note, in the inactive cells, Perm is zero, thus, 1/Perm results in inf, which then results in Kmean = 0!
        Kmean = 2. / (inv_perm[2,:, :, :-1] + inv_perm[2,:, :, 1:])
        #diff = pp[ip] - pp[ip + 1];
        #updir = (diff / dx) * x_dir_g_c - x_dir_g;
        grad = ht.diff(Press, axis=2)/self.Dx
        # + sign, because we later multiply by (-1.0)
        grad = grad * x_dir_g_c[:,:,:-1] + x_dir_g[:,:,:-1]

        flowright[:, :, :-1] = -1. * Kmean * grad * ht.where(grad > 0.0, Krel[:, :, 1:], Krel[:, :, :-1])
        flowleft[:,:,1:] = flowright[:,:,:-1]

        flowright *= self.Dy * self.Dz * Dzmult3D
        flowleft  *= self.Dy * self.Dz * Dzmult3D  # save this by setting flowleft after multiplication

        # Front and Back
        flowback = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)
        flowfront = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)

        #Note, in the inactive cells, Perm is zero, thus, 1/Perm results in inf, which then results in Kmean = 0!
        Kmean = 2. / (inv_perm[1,:, :-1, :] + inv_perm[1,:, 1:, :])
        #print('before diff Press.shape=', Press.gshape, 'Press.lshape', Press.lshape, flush=True)
        grad = ht.diff(Press, axis=1)/self.Dy
        # + sign, because we later multiply by (-1.0)
        grad = grad * y_dir_g_c[:,:-1,:] + y_dir_g[:,:-1,:]

        flowback[:, :-1, :] = -1. * Kmean * grad * ht.where(grad > 0.0, Krel[:, 1:, :], Krel[:, :-1, :])
        flowfront[:, 1:, :] = flowback[:, :-1, :]

        #print('before mul', (self.Dx * self.Dz * Dzmult3D).lshape, flowback.lshape, flush=True)
        flowback = flowback * self.Dx * self.Dz * Dzmult3D
        flowfront  = flowfront * self.Dx * self.Dz * Dzmult3D
        #print('after flows', flush=True)

        #  Top and Bottom
        flowtop = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)
        flowbottom = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)

        #Note, in the inactive cells, Perm is zero, thus, 1/Perm results in inf, which then results in Kmean = 0!
        Kmean = ( (Dzmult3D[:-1] + Dzmult3D[1:]) /
                (Dzmult3D[:-1]/self.Perm[0,:-1,:,:] + Dzmult3D[1:]/self.Perm[0,1:,:,:]) )

        grad = 1. + ht.diff(Press, axis=0) * 2. / (self.Dz * (Dzmult3D[:-1] + Dzmult3D[1:]))

        #Application of mask checks if node k is active
        flowtop[:-1, :, :] = -1. * Kmean * grad * ht.where(grad > 0.0, Krel[1:, :, :], Krel[:-1, :, :])
        flowbottom[1:, :, :] = flowtop[:-1, :, :]

        flowtop *= self.Dx * self.Dy
        flowbottom *= self.Dx * self.Dy
        #print('subsurfaceflow finished', flush=True)

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
        opahn = ht.float64(1.0) + (self.Alpha * ht.abs(Press))**self.Nvg
        ahnm1 = (self.Alpha * ht.abs(Press))**(self.Nvg-ht.float64(1.))
        Krel  = ht.where(Press<ht.float64(0.0), (ht.float64(1.0)-ahnm1 / (opahn)**m)**ht.float64(2.0) / opahn**(m/ht.float64(2.0)), ht.float64(1.0))
        return(Satur,Krel)

    def _SubsurfaceStorage(self, Press, Satur):
        shape3D = (self.Nz, self.Ny, self.Nx)
        subsurface_storage = ht.zeros(shape3D, dtype=ht.float64, split=self.Split3D)
        # for k in range(self.Nz):
        #     subsurface_storage[k,:,:] = Satur[k,:,:] * self.Poro[k,:,:] * self.Dx * self.Dy * self.Dz * self.Dzmult[k]
        #     subsurface_storage[k,:,:] += Press[k,:,:] * self.Sstorage[k,:,:] * Satur[k,:,:] * self.Dx * self.Dy * self.Dz * self.Dzmult[k]
        subsurface_storage[:] = (self.Poro + Press * self.Sstorage)
        subsurface_storage[:] *= Satur
        subsurface_storage[:] *= self.Dx * self.Dy * self.Dz
        subsurface_storage[:] *= ht.array(self.Dzmult, dtype=ht.float64).expand_dims(axis=-1).expand_dims(axis=-1)
        return(subsurface_storage)

    def _TopLayerPressure(self, Press, fill_val=99999.0):
        layers = (self.Mask > 0) * ht.arange(1, 1+self.Nz, dtype=ht.long).expand_dims(-1).expand_dims(-1)
        toplayer = layers.max(0) - 1
        #toplayer = ht.array(toplayer.larray, copy=False, is_split=self.Split2D)
        # toplayer contains the index of the highest layer and -1 if there is no highest layer
        y, x = np.indices(toplayer.larray.shape, sparse=True)  # sparse=True is important, otherwise x, y are unsplit(numpy) and of shape2D -> memory
        # do these need to be converted to heat tensors? -> No

        Toplayerpress = ht.array(Press.larray[toplayer.larray, y, x], copy=False, is_split=self.Split2D)
        Toplayerpress.larray[toplayer.larray < 0] = fill_val  # is this guaranteed to be balanced?
#         Toplayerpress = ht.where(toplayer.larray < 0, fill_val, Toplayerpress)
#         Toplayerpress[ht.nonzero(toplayer < 0)] = fill_val

        # alternative:
        # Toplayerpress = ht.full(shape2D, fill_val, split=self.Split)
        # Toplayerpress[toplayer >= 0] = Press[toplayer, y, x][toplayer >= 0]


        # check = ht.full(shape2D, -1, split=self.Split)
        # for k in reversed(range(self.Nz)):
        #     Toplayerpress[:,:] = ht.where((self.Mask[k,:,:]>0.0) & (check[:,:]<0), Press[k,:,:], Toplayerpress[:,:])
        #     #Check also contains the the layer index k of the top layer
        #     check[:,:] = ht.where((self.Mask[k,:,:]>0.0) & (check[:,:]<0), k, check[:,:])
        return Toplayerpress

    def _NetLateralOverlandFlow(self, overland_flow_x, overland_flow_y):
        shape2D = (self.Ny, self.Nx)
        Nix = ht.zeros(shape2D, split=self.Split2D)

        #Calc flow east
        #ParFlow:ke_[io] = pfmax(qx_[io], 0.0) - pfmax(-qx_[io + 1], 0.0);
        flow_east = ht.zeros_like(overland_flow_x)
        flow_east[:, :-1]  = ht.clip(overland_flow_x[:, :-1], min=0, max=None)
        flow_east[:, :-1] -= ht.clip(-1 * overland_flow_x[:, 1:], min=0, max=None)
        index_last = self.Nx-1 #slice(self.Nx-1, self.Nx)
        flow_east[:, index_last] = ht.where(overland_flow_x[:,index_last]>0.0, overland_flow_x[:,index_last], Nix[:,index_last])
        printroot('flow_east', flush=True)

        # for i in range (self.Nx-1):
        #     flow_east[:,i]  = ht.maximum(overland_flow_x[:,i], Nix[:,i])
        #     flow_east[:,i] -= ht.maximum((-1.0)*overland_flow_x[:,i+1], Nix[:,i+1])

        #Calc flow west
        #ParFlow:kw_[io] = pfmax(qx_[io - 1], 0.0) - pfmax(-qx_[io], 0.0);
        #flow_west = ht.zeros(shape2D, split=self.Split)
        flow_west = ht.zeros_like(overland_flow_x)
        flow_west[:, 1:]  = ht.clip(overland_flow_x[:, :-1], min=0, max=None)
        flow_west[:, 1:] -= ht.clip(-1 * overland_flow_x[:, 1:], min=0, max=None)
        tmp = ht.where(overland_flow_x[:,0]<0.0, overland_flow_x[:,0], Nix[:,0])
        #print(flow_west.gshape, 'flow_west:',flow_west[:, 0].gshape, flow_west[:, 0].lshape, flow_west[:, 0].split, 'tmp:', tmp.gshape, tmp.lshape, tmp.split, flush=True)
        flow_west[:, 0] = tmp #ht.where(overland_flow_x<0.0, overland_flow_x, Nix)[:,0:1]
        printroot('flow_west', flush=True)
        # for i in range (1,self.Nx):
        #     flow_west[:,i]  = ht.maximum(overland_flow_x[:,i-1], Nix[:,i-1])
        #     flow_west[:,i] -= ht.maximum((-1.0)*overland_flow_x[:,i], Nix[:,i])

        #Calc flow north
        #ParFlow:kn_[io] = pfmax(qy_[io], 0.0) - pfmax(-qy_[io + sy_p], 0.0);
        flow_north = ht.zeros(shape2D, split=self.Split2D)
        flow_north[:-1]  = ht.clip(overland_flow_y[:-1], min=0, max=None)
        flow_north[:-1] -= ht.clip(-1 * overland_flow_y[1:], min=0, max=None)
        index_last = self.Ny-1
        flow_north[index_last] = ht.where(overland_flow_y[index_last,:]>0.0, overland_flow_y[index_last,:], Nix[index_last,:])
        printroot('flow_north', flush=True)
        # for j in range (self.Ny-1):
        #     flow_north[j,:]  = ht.maximum(overland_flow_y[j,:], Nix[j,:])
        #     flow_north[j,:] -= ht.maximum((-1.0)*overland_flow_y[j+1,:], Nix[j+1,:])

        #Calc flow south
        #ParFlow:ks_[io] = pfmax(qy_[io - sy_p], 0.0) - pfmax(-qy_[io], 0.0);
        flow_south = ht.zeros(shape2D, split=self.Split2D)
        flow_south[1:]  = ht.clip(overland_flow_y[:-1], min=0, max=None)
        flow_south[1:] -= ht.clip(-1 * overland_flow_y[1:], min=0, max=None)
        flow_south[0] = ht.where(overland_flow_y[0]<0, overland_flow_y[0], Nix[0])
        #flow_south[0] = ht.clip(overland_flow_y[0], min=None, max=0)
        printroot('flow_south', flush=True)
        # for j in range (1,self.Ny):
        #     flow_south[j,:]  = ht.maximum(overland_flow_y[j-1,:], Nix[j-1,:])
        #     flow_south[j,:] -= ht.maximum((-1.0)*overland_flow_y[j,:], Nix[j,:])

        #Calc net lateral overland flow for each grid cell, (L^3/T)
        #ParFlow: ((ke_[io] - kw_[io]) / dx + (kn_[io] - ks_[io]) / dy
        # net_lateral_overlandflow = ht.zeros(shape2D, split=self.Split)
        net_lateral_overlandflow = self.Dy * (flow_east - flow_west) + self.Dx * (flow_north - flow_south)

        return(net_lateral_overlandflow)

    def _SubsurfaceFlow(self, Press, Krel):
        #ParFlow
        #x_dir_g = Mean(gravity * sin(atan(x_ssl_dat[io])), gravity * sin(atan(x_ssl_dat[io + 1])));
        #x_dir_g_c = Mean(gravity * cos(atan(x_ssl_dat[io])), gravity * cos(atan(x_ssl_dat[io + 1])));
        #y_dir_g = Mean(gravity * sin(atan(y_ssl_dat[io])), gravity * sin(atan(y_ssl_dat[io + sy_p])));
        #y_dir_g_c = Mean(gravity * cos(atan(y_ssl_dat[io])), gravity * cos(atan(y_ssl_dat[io + sy_p])));
        shape3D   = (self.Nz,self.Ny,self.Nx)
        x_dir_g   = ht.zeros(shape3D, split=self.Split3D)
        x_dir_g_c = ht.ones(shape3D, split=self.Split3D)
        y_dir_g   = ht.zeros(shape3D, split=self.Split3D)
        y_dir_g_c = ht.ones(shape3D, split=self.Split3D)
        if self.Terrainfollowing:
            printroot('Terrain following', flush=True)
            x_dir_g[:,:,:-1]   = (( ht.arctan(self.Slopex[0,:,:-1]) + ht.arctan(self.Slopex[0,:,1:]) )/2.0).expand_dims(0)
            x_dir_g_c[:,:,:-1] = (( ht.arctan(self.Slopex[0,:,:-1]) + ht.arctan(self.Slopex[0,:,1:]) )/2.0).expand_dims(0)
            y_dir_g[:,:-1,:]   = (( ht.arctan(self.Slopey[0,:-1,:]) + ht.arctan(self.Slopey[0,1:,:]) )/2.0).expand_dims(0)
            y_dir_g_c[:,:-1,:] = (( ht.arctan(self.Slopey[0,:-1,:]) + ht.arctan(self.Slopey[0,1:,:]) )/2.0).expand_dims(0)

        Dzmult3D = ht.array(self.Dzmult, dtype=ht.float64).expand_dims(axis=-1).expand_dims(axis=-1)
        inv_perm = 1.0 / self.Perm

        #Calculate the flux across the right and left face
        # Left and Right
        flowright = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)
        flowleft  = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)

        #Note, in the inactive cells, Perm is zero, thus, 1/Perm results in inf, which then results in Kmean = 0!
        Kmean = 2. / (inv_perm[2,:, :, :-1] + inv_perm[2,:, :, 1:])
        #diff = pp[ip] - pp[ip + 1];
        #updir = (diff / dx) * x_dir_g_c - x_dir_g;
        grad = ht.diff(Press, axis=2)/self.Dx
        # + sign, because we later multiply by (-1.0)
        grad = grad * x_dir_g_c[:,:,:-1] + x_dir_g[:,:,:-1]

        flowright[:, :, :-1] = -1. * Kmean * grad * ht.where(grad > 0.0, Krel[:, :, 1:], Krel[:, :, :-1])
        flowleft[:,:,1:] = flowright[:,:,:-1]

        flowright *= self.Dy * self.Dz * Dzmult3D
        flowleft  *= self.Dy * self.Dz * Dzmult3D  # save this by setting flowleft after multiplication

        # Front and Back
        flowback = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)
        flowfront = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)

        #Note, in the inactive cells, Perm is zero, thus, 1/Perm results in inf, which then results in Kmean = 0!
        Kmean = 2. / (inv_perm[1,:, :-1, :] + inv_perm[1,:, 1:, :])
        grad = ht.diff(Press, axis=1)/self.Dy
        # + sign, because we later multiply by (-1.0)
        grad = grad * y_dir_g_c[:,:-1,:] + y_dir_g[:,:-1,:]

        flowback[:, :-1, :] = -1. * Kmean * grad * ht.where(grad > 0.0, Krel[:, 1:, :], Krel[:, :-1, :])
        flowfront[:, 1:, :] = flowback[:, :-1, :]

        flowback *= self.Dx * self.Dz * Dzmult3D
        flowfront  *= self.Dx * self.Dz * Dzmult3D

        #  Top and Bottom
        flowtop = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)
        flowbottom = ht.zeros(shape3D,dtype=ht.float64,split=self.Split3D)

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

if __name__ == '__main__':
    pass
