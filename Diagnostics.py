# <editor-fold imports
import heat as ht
import numpy as np
from netCDF4 import Dataset
import calendar
import sys
import os

HOURLY = "hourly"
DAILY = "daily"
WEEKLY = "weekly"
MONTHLY = "monthly"
YEARLY = "yearly"
TIMES = [HOURLY, DAILY, WEEKLY, MONTHLY, YEARLY]


def printroot(*args, **kwargs):
    """ Only Root Process prints """
    if ht.MPI_WORLD.rank == 0:
        print(*args, **kwargs)
# </editor-fold>
# static functions first


def withoutNan(dndarray):
    """ not needed. Maybe balance resulting array. """
    return dndarray[ht.nonzero(dndarray == dndarray)]


def getDaysInMonth(year, month):
    month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    leap_month_days = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    monn = int(month) - 1
    days = month_days[monn]
    if (calendar.isleap(year)):
        days = leap_month_days[monn]
    return days


def load(dict, dir='.', split=None):
    """
    dict: files -> variables of that file;
         multiple variables can be given in a list
    TODO: - Allow different splitAxes
          - Add support for same Variable for every file
          - Introduce names, using filenames is not feasible for users
          - Combine names with Option to flatten output:
                    name -> dndarray
            (might use named tensors for it (torch 1.4))
    return:     nested dict: file -> variable -> dndarray
                access data using: data[file][variable]
    """
    data = {}
    for file, variables in dict.items():
        if isinstance(variables, str):  # allow single variable per file
            variables = [variables]
        if file not in data.keys():  # add new key
            data[file] = {}
        for var in variables:  # load all variables
            #assert type(var) is str, "variable name must be a string"
            path = file  # keep directory out of dict key
            if not file.startswith('/') or dir != '.':
                path = dir + '/' + file
            try:
                data[file][var] = ht.load_netcdf(
                    path, split=split, variable=var)
            except ValueError:
                data[file][var] = ht.load_netcdf(path, split=0, variable=var)
    return data


def loadMerge(dict, dir='.', mergeAxis=0, keepDims=True, split=None):
    """
    Loads multiple Files into **one** DNDarray.
    e.g. load monhly files into a yearly DNDarray.
    The Datasets are concatenated along the mergeAxis.
    If keepDims is False, a new dimension is added along which the files are
    concatenated. This dimension is added at the position of the mergeAxis.

    dict: files -> variables of that file;
         multiple variables can be given in a list
    TODO: Add support for same Variable for every file
    returns: DNDarray
    """
    # load
    if mergeAxis < split and not keepDims:
        """ If keepDims is False and and mergeAxis < split, the splitAxis of
        the result will be split +1 (Because a new Axis is added "in front of"
        the splitAxis). Therefore this correction."""
        dic = load(dict, dir, split - 1)
    else:
        dic = load(dict, dir, split)

    shape = None
    for vars in dic.values():
        for data in vars.values():
            if shape is None:  # initialize shape
                shape = list(data.shape)
                if not keepDims:
                    shape.insert(mergeAxis, 1)
            else:
                if keepDims:
                    shape[mergeAxis] += data.shape[mergeAxis]
                else:
                    shape[mergeAxis] += 1

    index = [slice(None) for _ in shape]
    ret = ht.empty(shape, split=split)

    # merge
    offset = 0
    for vars in dic.values():
        for data in vars.values():
            if keepDims:
                index[mergeAxis] = slice(
                    offset, offset + data.shape[mergeAxis])
                offset += data.shape[mergeAxis]
            else:
                index[mergeAxis] = slice(offset, offset + 1)
                offset += 1
            ret[tuple(index)] = data

    """ # ht.concatenate only supports two dndarrays
    dndarrays = tuple(
        data for variables in dic.values() for data in variables.values()
        )
    return ht.concatenate(dndarrays, mergeAxis)
    """
    return ret


def writeVarsToNCFile(dict, filename, missval,
                      lati, loni, delta,
                      ntimes, gz, nlats, nlons,  # use shape for this?
                      tstart=0, slices=[slice(None)], year=None, month=None):
    """
    dic : {varName : (unit, dndarray)}
    For 3D-dndarrays, the shape must always be (time, level, lat, lon)
    For 2D-dndarrays, the shape must always be (time, lat, lon)

    TODO:   How to implement time? Different timesteps can be written at
                different times (try as currently is, otherwise first write
                the whole timeaxis).
            Add metadata (units etc)
    Issues: HeAT automatically generates Dimensions for every variable and
    Dimensions cannot be changed later. Thereby every variable creates their own
    dimensions.
    An option to pass the dimension_names to ht.save_netcdf could work, but needs
    to ensure that every dimension exists (or is only available in 'r+' mode and
    lets user handle errors)
    -- solved locally in env: PythonHeAT_2019-12-19_1d1c514
    """
    dims3D = ('time', 'level', 'lon', 'lat')
    dims2D = ('time', 'lon', 'lat')
    if not os.path.isfile(filename):  # generate construct of file
        # does this append?
        time_out = tstart + ht.arange(ntimes, dtype=ht.float32, split=0)
        lats_out = lati + delta * ht.arange(nlats, dtype=ht.float32, split=0)
        lons_out = loni + delta * ht.arange(nlons, dtype=ht.float32, split=0)
        lev_out = 1.0 + ht.arange(gz, dtype=ht.float32, split=0)
        ht.save_netcdf(time_out, filename, 'time', 'w',
                       ['time'], is_unlimited=True)
        # Does this override existing variables?
        ht.save_netcdf(lats_out, filename, 'lat', 'r+', ['lat'])
        # If so, distinguish between writing and appending
        ht.save_netcdf(lons_out, filename, 'lon', 'r+', ['lon'])
        ht.save_netcdf(lev_out, filename, 'level', 'r+', ['level'])

    for varName, [unit, data] in dict.items():
        printroot('writing', varName, 'is balanced:',
                  data.is_balanced(), flush=True)
        data.balance_()
        if len(data.shape) == 4:  # 4 because time is included
            printroot('check shapes:', data.shape,
                      (ntimes, gz, nlons, nlats), flush=True)
            dim_names = dims3D
        else:
            dim_names = dims2D
        ht.save_netcdf(data, filename, varName, 'r+', dim_names, file_slices=slices,
                       fill_value=missval, zlib=False, least_significant_digit=6)  # zlib should be True
    return


class Diagnostics():  # Make this a subclass of ht.DNDarray?
    """
    Diagnostics is the class for processing ParFlow-netCDF4 output.
    TODO:   Parallelize for-loops by iterating through lloc[:] -> make sure all
                data is on local process (would this be useful?)
            Make methods return a dict (same as loading method) rather than
                a tuple -> Consistent to the writeVarsToNCFile and load method.
            Original scripts work by processing one timestep at a time and
                appending the output. storage and velocities currently 'yield'
                every timestep. Add a wrapper to either store every timestep or
                write every timestep.
    Notes:  How about loading all data in init?
            Should every Method automatically write the output?
            Are distributed masks useful in HeAT? Or is it better to not split
                them (load landmask on all procs)? Here, masks are used through
                ht.nonzero(arr), which has split=0 if arr.split!=None.Therefore
                every proc receives different indices -> Does every proc
                receive its local indices? If every proc receives all indices,
                do they automatically choose their local ones for processing?
                (probably yes)
            Which split is recommended? Maybe different for different methods
                storage iterates through 0,1
                deltaTWS through 0
                velocities through 0 and 3
            At which points is it smart to rebalance the data? e.g. landmasks
                lead to imbalances.
    """
    # non-static functions

    def __init__(self, year, month, ncDir, lmfile, satfile, pressfile,
                 specstorfile, timesteps, porfile, lati, loni, delta, dx, dy,
                 dz_list, gx, gy, gz, outDir=None, data=None, timeResolution=None,
                 timeAxis=None, split=None,
                 missval=-340282346638528859811704183484516925440):
        """
        Is it adequate to use missval instead of nan, because heat lacks support
        of nan-safe methods? (are nan-safe methods needed?)
        Does specifying outDir make sense here?
        -> for convienience: if given, write output
        Use Hard-Coded variables for ParFlow output Variables instead
        Every Diagnostic has:
            - Land-mask
            - Saturation
            - porosity
            - Time: year, month, timesteps
            - Lat; Lon
            - inDirectory, outDirectory
            - delta <- ?
            - d and g cooridnates <- ?

        files = [satfile, pressfile, specstorfile, lmfile, porfile]
        variables = [['SATUR'],['time', 'PSI'],['SPEC_STO'],['LAND_MASK'],['PORO']]
        """
        self.data = None  # store all loaded variables in data
        if data is not None:  # This doesnt make sense, provide files:variables dict instead
            self.data = ht.array(data, split=split)

        self.timeResolution = None
        if timeResolution in TIMES:
            self.timeResolution = timeResolution

        self.timeAxis = timeAxis
        self.split = split
        self.missval = missval
        self.year = year
        self.month = month
        self.outDir = outDir
        self.timesteps = int(timesteps)
        self.ncDir = ncDir
        self.porfile = porfile
        self.lmfile = lmfile
        self.satfile = satfile
        self.pressfile = pressfile
        self.specstorfile = specstorfile
        self.lati = float(lati)
        self.loni = float(loni)
        self.delta = float(delta)
        # input list here
        dz_list = dz_list.strip('[]')
        self.dz = [float(c) for c in dz_list.split(',')]
        self.dx = float(dx)
        self.dy = float(dy)
        self.gx = int(gx)
        self.gy = int(gy)
        self.gz = int(gz)

    def writeToNC(self, dict, filename, tstart=0, slices=[slice(None)]):
        writeVarsToNCFile(dict, filename, self.missval,
                          self.lati, self.loni, self.delta,
                          self.timesteps, self.gz, self.gx, self.gy,
                          tstart, slices)

    def aggregate(self, func, axis=None, time_resolution=None):
        """ NOT USABLE AT THIS POINT
        reduces the temporal resolution.
        should support different methods, e.g. summarize using mean or first
        value of every new timestep.
        TODO: integrate time_resolution, better implementation of timeAxis
        """
        assert self.timeResolution is not None, "no current temporal resolution set"
        assert time_resolution in TIMES, "no usable temporal resolution provided"

        if axis is None:
            axis = self.timeAxis
        self.data = func(self.data, axis)
        return self

    def storage(self):
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

        printroot(tstart)
        printroot(tend)
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

            yield {
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

    def soilMoisture(self):
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

        return {
            'SM_Surf': ('m3/m3', arraysurf),
            'SM_Root': ('m3/m3', arrayroot)
        }

    def deltaTWS(self, storagesfile):
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
        return {'MONTLY_DELTA_TWS': ('mm', deltatws)}

    def velocities(self, velxfile, velyfile, velzfile):
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
            yield {
                'LATERALFLOW_MAG': ('m/hour', lateralvel_mag),
                'VERTICALFLOW_MAG': ('m/hour', verticalvel_mag),
                'VEL_MAG': ('m/hour', vel_mag),
                'VEL_X': ('m/hour', vel_x),
                'VEL_Y': ('m/hour', vel_y),
                'VEL_Z': ('m/hour', vel_z)
            }


if __name__ == '__main__':
    # for storage
    data = Diagnostics(year=sys.argv[1],
                       month=sys.argv[2],
                       ncDir=sys.argv[3],
                       lmfile=sys.argv[4],
                       satfile=sys.argv[5],
                       pressfile=sys.argv[6],
                       specstorfile=sys.argv[7],
                       timesteps=sys.argv[8],
                       porfile=sys.argv[9],
                       lati=sys.argv[10],
                       loni=sys.argv[11],
                       delta=sys.argv[12],
                       dx=sys.argv[13],
                       dy=sys.argv[14],
                       dz_list=sys.argv[15],
                       gx=sys.argv[16],
                       gy=sys.argv[17],
                       gz=sys.argv[18],
                       outDir=sys.argv[19],
                       timeAxis=0, split=-1)

    # calc soil moisture at a grid cell/point
    for timestep in data.storage():
        data.writeToNC(timestep, 'storages.nc')
        printroot('writing finished')
        #raise SystemExit()
