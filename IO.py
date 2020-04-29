import heat as ht
import numpy as np
from struct import pack, unpack


def read_pfb(filename, split=None):
    with open(filename, "rb") as f:
        # read meta informations of datafile
        meta_inf = np.fromfile(f, dtype='>f8', count=3)
        x1 = meta_inf[0]
        y1 = meta_inf[1]
        z1 = meta_inf[2]

        meta_inf = np.fromfile(f, dtype='>i4', count=3)
        nx = meta_inf[0]
        ny = meta_inf[1]
        nz = meta_inf[2]
        nn = nx * ny * nz

        meta_inf = np.fromfile(f, dtype='>f8', count=3)
        dx = meta_inf[0]
        dy = meta_inf[1]
        dz = meta_inf[2]

        meta_inf = np.fromfile(f, dtype='>i4', count=1)
        nsubgrid = meta_inf[0]

        data = np.ndarray(shape=(nz, ny, nx), dtype='>f8')

        for s in range(nsubgrid):
            meta_inf = np.fromfile(f, dtype='>i4', count=9)
            ix = meta_inf[0]
            iy = meta_inf[1]
            iz = meta_inf[2]
            # print("---{0} Start Index (X,Y,Z):".format(s+1), ix, iy, iz)

            nx = meta_inf[3]
            ny = meta_inf[4]
            nz = meta_inf[5]
            nn = nx * ny * nz
            # print("---{0} Dimensions (X,Y,Z):".format(s+1), nx, ny, nz)

            rx = meta_inf[6]
            ry = meta_inf[7]
            rz = meta_inf[8]
            # print("---{0} Offsets (X,Y,Z):".format(s+1), rx, ry, rz)

            tmp_data = np.fromfile(
                f, dtype='>f8', count=nn).reshape((nz, ny, nx))

            data[iz:iz + nz, iy:iy + ny, ix:ix + nx] = tmp_data

    return ht.array(data, split=split)


def read_nc4(dict, dir='.', split=None):
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
    return_dict = {}
    for file, variables in dict.items():
        if isinstance(variables, str):  # allow single variable per file
            variables = [variables]
        if file not in return_dict.keys():  # add new key
            return_dict[file] = {}
        for var in variables:  # load all variables
            # assert type(var) is str, "variable name must be a string"
            path = file  # keep directory out of dict key
            if not file.startswith('/') or dir != '.':
                path = dir + '/' + file
            try:
                return_dict[file][var] = ht.load_netcdf(
                    path, split=split, variable=var)
            except ValueError:
                return_dict[file][var] = ht.load_netcdf(
                    path, split=0, variable=var
                    )
    return return_dict
