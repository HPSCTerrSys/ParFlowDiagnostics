from Diagnostics import *

class Storage(Diagnostic):

    def calculate(self):
        """
        TODO: - Move loading and saving out of this Method
            - make static with all data as input and provide non-static wrapper
        ignore g-coordinates for now as heat does not support reshape and they
        should equal the dndarrays shape.
        """

        # We are assuming monthly time scales (monmean)
        # TODO enforce monthly means by using aggregate
        # load vars
        dic = {
            self.satfile: ['SATUR', 'time'],
            self.pressfile: ['PSI'],  # 'time', 'PSI'],
            self.specstorfile: 'SPEC_STO_GEOL1',  # 'SPEC_STO',
            self.lmfile: 'LAND_MASK',
            self.porfile: 'PORO'
        }
        loaded = load(dic, self.ncDir, self.split)
        sat_all = loaded[self.satfile]['SATUR']
        times = loaded[self.satfile]['time']
        press_all = loaded[self.pressfile]['PSI']
        specsto = loaded[self.specstorfile]['SPEC_STO_GEOL1']
        # somehow, squeeze gives an mpiError
        lm = loaded[self.lmfile]['LAND_MASK'][0]
        por = loaded[self.porfile]['PORO'][0]

        printroot('loading finished', flush=True)
        times = ht.arange(len(times), dtype=ht.float32,
                          split=None)  # self.split)
        # lm = np.reshape(lm, shape3D_notime)
        printroot("compare shapes:", lm.shape, (self.gz, self.gy, self.gx))
        lm2D = lm[self.gz - 1]
        mask = lm < 1
        mask2D = lm2D < 1

        del(lm, lm2D)

        # por = np.reshape(por, shape3D_notime)
        # specsto = np.reshape(specsto, (gz, gy, gx))
        """
        tstart = []
        tend = []
        for i in range(0, tottimesteps - 1, 96):
            # go in incremements of 96
            tstart.append(i)
            if (i < (tottimesteps - 96)):
                tend.append(i + 96)
            else:
                tend.append(tottimesteps)
        """
        tstart = ht.arange(0, self.timesteps - 1, 96)
        tend = ht.empty_like(tstart)
        tend[:-1] = tstart[1:]
        tend[-1] = self.timesteps

        #printroot(tstart)
        #printroot(tend)
        for i in range(len(tstart)):
            timesteps = (tend[i] - tstart[i]).item()
            printroot("timesteps %d" % timesteps, flush=True)
            shape2D = (timesteps, self.gy, self.gx)
            shape2D_notime = (self.gy, self.gx)
            shape3D = (timesteps, self.gz, self.gy, self.gx)
            shape3D_notime = (self.gz, self.gy, self.gx)

            # Some kind of factory for generating these?? No?
            sat_stor_3D = ht.zeros(shape3D, split=self.split)
            unsat_stor_3D = ht.zeros(shape3D, split=self.split)
            tot_stor_3D = ht.zeros(shape3D, split=self.split)
            capi_stor_3D = ht.zeros(shape3D, split=self.split)
            tws_3D = ht.zeros(shape3D, split=self.split)
            pond_stor_2D = ht.zeros(shape2D, split=self.split)
            surf_stor_sat_2D = ht.zeros(shape2D, split=self.split)
            surf_stor_unsat_2D = ht.zeros(
                shape2D, split=self.split)
            surf_stor_2D = ht.zeros(shape2D, split=self.split)
            tws_2D = ht.zeros(shape2D, split=self.split)
            sat_stor_3D_array = ht.zeros(shape3D_notime, split=self.split)
            unsat_stor_3D_array = ht.zeros(shape3D_notime, split=self.split)
            tot_stor_3D_array = ht.zeros(shape3D_notime, split=self.split)
            capi_stor_3D_array = ht.zeros(shape3D_notime, split=self.split)
            tws_3D_array = ht.zeros(shape3D_notime, split=self.split)
            pond_stor_2D_array = ht.zeros(shape2D_notime, split=self.split)
            surf_stor_sat_2D_array = ht.zeros(shape2D_notime, split=self.split)
            surf_stor_unsat_2D_array = ht.zeros(
                shape2D_notime, split=self.split)
            surf_stor_2D_array = ht.zeros(shape2D_notime, split=self.split)
            tws_2D_array = ht.zeros(shape2D_notime, split=self.split)

            printroot('memory allocated', flush=True)
            # surf_ss = np.reshape(specsto[gz - 1, :, :], shape2D_notime)
            surf_ss = specsto[self.gz - 1, :, :]
            # surf_poro = np.reshape(por[gz - 1, :, :], shape2D_notime)
            surf_poro = por[self.gz - 1, :, :]
            for j in range(timesteps):  # idx = j
                press = press_all[j, :, :, :]
                # press = np.reshape(press, shape3D_notime)
                sat = sat_all[j, :, :, :]
                # sat = np.reshape(sat, shape3D_notime)
                # surface where gz=15 is the surface
                # surf_sat = np.reshape(sat[gz - 1, :, :], shape2D_notime)
                surf_sat = sat[self.gz - 1, :, :]
                # surf_press = np.reshape(press[gz - 1, :, :], shape2D_notime)
                surf_press = press[self.gz - 1, :, :]
                # get surface indicies where press>=0:
                # ind_array = np.ma.masked_where(surf_press>=0,surf_press)
                mask_array = surf_press >= 0.0
                mask_array.astype(ht.int, False)
                # print mask_array
                # make an index array of 1s and 0s /p/project/cjibg31/jibg3101/SharedData/all_scripts/analysis_pfl_clm_va
                # mask_array = np.zeros((gy,gz))
                # mask_array[ind_array.mask] = 1
                surf_sat_sat = surf_sat * mask_array
                surf_press_sat = surf_press * mask_array
                inverse = 1 - mask_array
                surf_sat_unsat = surf_sat * inverse
                surf_press_unsat = surf_press * inverse
                # print np.shape(surf_sat_sat)
                pond_stor_2D_array = surf_press_sat * self.dx * self.dy
                surf_stor_sat_2D_array = surf_sat_sat * surf_poro * self.dx * self.dy + \
                    surf_ss * surf_poro * surf_press_sat * surf_sat_sat * self.dx * self.dy
                surf_stor_unsat_2D_array = surf_sat_unsat * surf_poro * self.dx * self.dy + \
                    surf_ss * surf_poro * surf_press_unsat * surf_sat_unsat * self.dx * self.dy
                # take the whole field no masks
                surf_stor_2D_array = surf_sat * surf_poro * self.dx * self.dy + \
                    surf_ss * surf_poro * surf_press * surf_sat * self.dx * self.dy
                tws_2D_array = 1000.0 * \
                    (surf_stor_2D_array + pond_stor_2D_array) / (self.dx * self.dy)

                del(surf_sat_sat, surf_press_sat, surf_sat_unsat, surf_press_unsat,
                    surf_sat, surf_press, mask_array)

                # sub surface
                for k in range(self.gz):
                    printroot('subsurface storage level:', k, flush=True)
                    sat_levk = sat[k, :, :]
                    press_levk = press[k, :, :]
                    por_levk = por[k, :, :]
                    ss_levk = specsto[k, :, :]
                    mask_array_2a = press_levk >= 0
                    mask_array_2a.astype(ht.int, False)
                    # print mask_array_2a
                    mask_array_2b = (press_levk < 0) & (sat_levk > 0.99)
                    mask_array_2b.astype(ht.int, False)
                    sat_levk_m1 = sat_levk * mask_array_2a
                    sat_levk_invm1 = sat_levk * (1 - mask_array_2a)
                    sat_levk_m2 = sat_levk * mask_array_2b
                    press_levk_m1 = press_levk * mask_array_2a
                    press_levk_invm1 = press_levk * (1 - mask_array_2a)
                    press_levk_m2 = press_levk * mask_array_2b
                    # NOTE: Does this indexing require the 3D variables to have
                    # a different split-Axis?
                    if(k != (self.gz - 1)):
                        sat_stor_3D_array[k, :, :] = sat_levk_m1 * por_levk * self.dx * self.dy * self.dz[k] + \
                            por_levk * press_levk_m1 * ss_levk * \
                            sat_levk_m1 * self.dx * self.dy * self.dz[k]
                    else:
                        unsat_stor_3D_array[k, :, :] = sat_levk_invm1 * por_levk * self.dx * self.dy * self.dz[k] + \
                            por_levk * press_levk_invm1 * ss_levk * \
                            sat_levk_invm1 * self.dx * self.dy * self.dz[k]
                    # take the whole field no masks
                    tot_stor_3D_array[k, :, :] = sat_levk * por_levk * self.dx * self.dy * \
                        self.dz[k] + por_levk * press_levk * \
                        ss_levk * sat_levk * self.dx * self.dy * self.dz[k]
                    tws_3D_array[k, :, :] = 1000.0 * \
                        tot_stor_3D_array[k, :, :] / (self.dx * self.dy)
                    capi_stor_3D_array[k, :, :] = sat_levk_m2 * por_levk * self.dx * self.dy * self.dz[k] + \
                        por_levk * press_levk_m2 * ss_levk * \
                        sat_levk_m2 * self.dx * self.dy * self.dz[k]

                    del(sat_levk, press_levk, por_levk, ss_levk, sat_levk_m1,
                        press_levk_m1, sat_levk_m2, press_levk_m2, sat_levk_invm1,
                        press_levk_invm1, mask_array_2a, mask_array_2b)

                # mask array here
                # sat = None
                # press = None
                printroot('masking timestep:', j, flush=True)
                mask_ind = ht.nonzero(mask)
                mask_ind_2D = ht.nonzero(mask2D)
                sat_stor_3D_array[mask_ind] = ht.nan
                sat_stor_3D[j, :, :, :] = sat_stor_3D_array
                unsat_stor_3D_array[mask_ind] = ht.nan
                unsat_stor_3D[j, :, :, :] = unsat_stor_3D_array
                tot_stor_3D_array[mask_ind] = ht.nan
                tot_stor_3D[j, :, :, :] = tot_stor_3D_array
                capi_stor_3D_array[mask_ind] = ht.nan
                capi_stor_3D[j, :, :, :] = capi_stor_3D_array
                tws_3D_array[mask_ind] = ht.nan
                tws_3D[j, :, :, :] = tws_3D_array
                surf_stor_2D_array[mask_ind_2D] = ht.nan
                surf_stor_2D[j, :, :] = surf_stor_2D_array
                surf_stor_sat_2D_array[mask_ind_2D] = ht.nan
                surf_stor_sat_2D[j, :, :] = surf_stor_sat_2D_array
                pond_stor_2D_array[mask_ind_2D] = ht.nan
                pond_stor_2D[j, :, :] = pond_stor_2D_array
                surf_stor_unsat_2D_array[mask_ind_2D] = ht.nan
                surf_stor_unsat_2D[j, :, :] = surf_stor_unsat_2D_array
                tws_2D_array[mask_ind_2D] = ht.nan
                tws_2D[j, :, :] = tws_2D_array

            self.data = {
                'SUB_SURF_SAT_STO': ('m3', sat_stor_3D),
                'SUB_SURF_UNSAT_STO': ('m3', unsat_stor_3D),
                'SUB_SURF_CAPI_STO': ('m3', capi_stor_3D),
                'SUB_SURF_TOT_STO': ('m3', tot_stor_3D),
                'SUB_SURF_TWS': ('mm', tws_3D),
                'SURF_SAT_STO': ('m3', surf_stor_sat_2D),
                'SURF_TOT_STO': ('m3', surf_stor_2D),
                'POND_STO': ('m3', pond_stor_2D),
                'SURF_UNSAT_STO': ('m3', surf_stor_unsat_2D),
                'SURF_TWS': ('mm', tws_2D)
            }
            yield self.data


if __name__ == '__main__':
    storage = Storage(*sys.argv[1:], timeAxis=0, split=-1)
    for timestep in storage.calculate():
        storage.writeToNC(filename='storages.nc', slices=[slice(None, -1, -1)])
