from Diagnostics import *


class Velocities(Diagnostic):
    def calculate(self, velxfile, velyfile, velzfile):
        # Calculate_LateralGroundwaterFlow_viaVelocity
        #debug = 'wplotpor'
        dic = {
            self.lmfile: 'LAND_MASK',
            self.porfile: 'PORO',
            velxfile: ['VEL_X', 'time'],
            velyfile: 'VEL_Y',
            velzfile: 'VEL_Z',
        }
        loaded = load(dic, self.ncDir, self.split)
        lm = loaded[self.lmfile]['LAND_MASK']
        por = loaded[self.porfile]['PORO']
        time = loaded[velxfile]['time']
        arrvelx = loaded[velxfile]['VEL_X']
        arrvely = loaded[velyfile]['VEL_Y']
        arrvelz = loaded[velzfile]['VEL_Z']

        # We are calculating overland flow at EACH grid cell on the surface
        lm2D = lm[self.gz - 1, :, :]
        mask = lm < 1
        por[ht.nonzero(mask)] = ht.nan
        times = ht.arange(len(time), dtype='float32')
        del time
        tstart = ht.arange(0, self.timesteps - 1, 96)
        tend = ht.empty_like(tstart)
        tend[:-1] = tstart[1:]
        tend[-1] = self.timesteps
        printroot(tstart)
        printroot(tend)
        for i in range(len(tstart)):
            timesteps = tend[i] - tstart[i]
            printroot("timesteps %d".format(timesteps))

            velx = arrvelx[tstart[i]:tend[i], :, :, 1:self.gx + 1]
            vely = arrvely[tstart[i]:tend[i], :, :, 1:self.gy + 1]
            velz = arrvelz[tstart[i]:tend[i], :, :, 1:self.gz + 1]

            lm4D = ht.empty((timesteps, *lm.shape),
                            split=None if lm.split is None else lm.split + 1)
            for i in range(timesteps):
                lm4D[i] = lm
            mask = lm4D < 1

            # make the velocity arrays:
            lateralvel_mag = ht.sqrt(
                por * ht.power(velx, 2) + por * ht.power(vely, 2))
            lateralvel_mag[ht.nonzero(mask)] = ht.nan
            verticalvel_mag = ht.sqrt(por * ht.power(velz, 2))
            verticalvel_mag[ht.nonzero(mask)] = ht.nan
            vel_mag = ht.sqrt(por * ht.power(velx, 2) + por *
                              ht.power(vely, 2) + por * ht.power(velz, 2))
            vel_mag[ht.nonzero(mask)] = ht.nan
            vel_x = por * velx
            vel_x[ht.nonzero(mask)] = ht.nan
            vel_y = por * vely
            vel_y[ht.nonzero(mask)] = ht.nan
            vel_z = por * velz
            vel_z[ht.nonzero(mask)] = ht.nan
            self.data = {
                'LATERALFLOW_MAG': ('m/hour', lateralvel_mag),
                'VERTICALFLOW_MAG': ('m/hour', verticalvel_mag),
                'VEL_MAG': ('m/hour', vel_mag),
                'VEL_X': ('m/hour', vel_x),
                'VEL_Y': ('m/hour', vel_y),
                'VEL_Z': ('m/hour', vel_z)
            }
            yield self.data
