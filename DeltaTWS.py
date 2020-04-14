from Diagnostics import *

class DeltaTWS(Diagnostic):
    def calculate(self, storagesfile):
        """
        NOTES: is self.timesteps the number of months or timesteps?
        TODO:
        """
        dic = {
            self.lmfile: 'LAND_MASK',
            storagesfile: ['SURF_TWS', 'SUB_SURF_TWS'],
        }
        loaded = load(dic, self.ncDir, self.split)
        lm = loaded[self.lmfile]['LAND_MASK']
        surf_tws = loaded[storagesfile]['SURF_TWS']
        subsurf_tws = loaded[storagesfile]['SUB_SURF_TWS']

        lm2D = lm[self.gz - 1, :, :]
        mask = lm < 1
        mask2D = lm2D < 1
        del(lm, lm2D)
        tstart = []
        tend = []
        times = []  # what do we need this for?

        summ = 0
        # iterate through months (is timesteps the number of months or the number of timesteps??)
        for i in range(self.timesteps):
            num_hours = 24 * getDaysInMonth(self.year, i)
            tstart.append(0)
            times.append(summ)
            tend.append(num_hours - 1)
            summ += num_hours

        deltatws = ht.zeros(
            (self.timesteps, self.gy, self.gx), split=self.split)
        for i in range(len(tstart)):
            timestep_s = tstart[i]
            timestep_e = tend[i]
            #month = i + 1
            #month = str(month).zfill(2)

            tws_start = subsurf_tws[timestep_s]
            tws_end = subsurf_tws[timestep_e]
            surf_tws_start = surf_tws[timestep_s]
            surf_tws_end = surf_tws[timestep_e]

            beg_tws = ht.sum(tws_start, axis=0) + surf_tws_start
            end_tws = ht.sum(tws_end, axis=0) + surf_tws_end
            diff = end_tws - beg_tws
            diff[ht.nonzero(mask2D)] = ht.nan
            deltatws[i] = diff
        self.data = {'MONTLY_DELTA_TWS': ('mm', deltatws)}
        return self.data
