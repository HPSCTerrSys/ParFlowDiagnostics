from Diagnostics import *

class SoilMoisture(Diagnostic):
    def calculate(self):
        #debug = 'wplotpor'
        dic = {
            # satfile: 'time',
            self.lmfile: 'LAND_MASK',
            self.porfile: 'PORO'
        }
        loaded = load(dic, self.ncDir, self.split)
        lm, por = loaded[self.lmfile]['LAND_MASK'], loaded[self.porfile]['PORO']
        # load saturation at all timesteps
        dic = {
            self.satfile + str(i + 1).zfill(2) + '.nc': 'SATUR' for i in range(self.timesteps)
        }
        sat_all_times = loadMerge(
            dic, self.ncDir, mergeAxis=0, keepDims=False, split=self.split)
        # this broadcast should work if mergeAxis==0
        array = ht.multiply(sat_all_times, por)
        # sum up times from all timesteps, is this correct?
        dic = {
            self.satfile + str(i + 1).zfill(2) + '.nc': 'time' for i in range(self.timesteps)
        }
        times = loadMerge(dic, self.ncDir, mergeAxis=0,
                          keepDims=False, split=self.split).sum(0)

        arraysurf = array[:, 13:15, :, :]
        arrayroot = array[:, 4:15, :, :]
        array = None
        arraysurf = ht.mean(arraysurf, axis=1)
        arrayroot = ht.mean(arrayroot, axis=1)
        printroot(arrayroot.shape)

        """#  where does por_nan come from??
        if (debug == 'plotpor') and ht.MPI_WORLD.rank == 0:
            por_nan2d = por_nan[(gz - 1), :, :]
            por_nan2d = np.reshape(por_nan2d, shape2D_notime)
            ax.scatter(x2d[skip2d], y2d[skip2d],
                       por_nan2d[skip2d], color='y', s=0.5)
            ax.set_xlabel('X Axis')
            ax.set_ylabel('Y Axis')
            ax.set_zlabel('Porosity')
            ax.set_zlim(0, 1.0)
            plt.show()
        #"""

        self.data = {
            'SM_Surf': ('m3/m3', arraysurf),
            'SM_Root': ('m3/m3', arrayroot)
        }
        return self.data
