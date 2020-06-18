import heat as ht
import numpy as np
from struct import pack, unpack
import os


def write_packed(f, fmt, val):
    f.write(pack(fmt, val))


def create_pfb(filename, dndarray, delta=(1, 1, 1), subgrids=(1, 1, 1)):
    # dndarray = dndarray.astype(ht.float32)
    var = dndarray.numpy()
    if dndarray.comm.rank == 0:
        print(var.shape)
        nz, ny, nx = var.shape
        dz, dy, dx = delta
        sz, sy, sx = subgrids

        filepfb = open(filename, "wb")

        # Write start indices of global domain in x, y, z direction
        write_packed(filepfb, ">d", 0)
        write_packed(filepfb, ">d", 0)
        write_packed(filepfb, ">d", 0)

        # Write number of global gridpoints in x, y, z direction
        write_packed(filepfb, ">i", nx)
        write_packed(filepfb, ">i", ny)
        write_packed(filepfb, ">i", nz)

        # Write delta x, delta y and delta z
        write_packed(filepfb, ">d", dx)
        write_packed(filepfb, ">d", dy)
        write_packed(filepfb, ">d", dz)

        nSubGrid = np.prod(subgrids)

        nnx = int(nx / sx)
        nny = int(ny / sy)
        nnz = int(nz / sz)

        # Write the subgrid grid ID
        write_packed(filepfb, ">i", nSubGrid)

        for iz in np.arange(sz) * nnz:
            for iy in np.arange(sy) * nny:
                for ix in np.arange(sx) * nnx:
                    # print(ix,iy,iz, nnx,nny,nnz)
                    # Write start indices in x, y, z direction
                    write_packed(filepfb, ">i", int(ix))
                    write_packed(filepfb, ">i", int(iy))
                    write_packed(filepfb, ">i", int(iz))

                    # Write number of grid points in x, y and z direction for this subgrid
                    write_packed(filepfb, ">i", nnx)
                    write_packed(filepfb, ">i", nny)
                    write_packed(filepfb, ">i", nnz)

                    # Write the relative(to global) grid refinement in this subgrid
                    # 0=same resolution as global
                    write_packed(filepfb, ">i", 0)
                    write_packed(filepfb, ">i", 0)
                    write_packed(filepfb, ">i", 0)

                    # Assuming the data is stored in 3D array called varArray of global size nx*ny*nz
                    fmt = ">%dd" % (nnx * nny * nnz)
                    filepfb.write(
                        pack(fmt, *var[iz : iz + nnz, iy : iy + nny, ix : ix + nnx].flatten())
                    )

        filepfb.close()


# parflow reader


def read_packed(f, fmt, size):
    return unpack(fmt, f.read(size))


def read_pfb(filename, dtype=">f8", split=None, comm=ht.MPI_WORLD):
    """Reads .pfb file into a dndarray. Reading is only done by the
    root process and the data gets distributed afterwards.

    Parameters
    ----------
    filename : str
        File to read.
    split : int or None
        split-Axis of resulting dndarray (the default is None).
    comm : MPI.Communicator
        Communicator of resulting dndarray (the default is ht.MPI_WORLD).
    dtype : str or numpy dtype
        Dtype of the data in the pfb file (the default is ">f8").

    Returns
    -------
    dndarray
        Data in file.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    comm = ht.communication.sanitize_comm(comm)
    if comm.rank == 0:  # root-process does reading
        with open(filename, "rb") as f:
            # read meta informations of datafile
            meta_inf = np.fromfile(f, dtype=">f8", count=3)
            x1 = meta_inf[0]
            y1 = meta_inf[1]
            z1 = meta_inf[2]

            meta_inf = np.fromfile(f, dtype=">i4", count=3)
            nx = meta_inf[0]
            ny = meta_inf[1]
            nz = meta_inf[2]
            nn = nx * ny * nz

            meta_inf = np.fromfile(f, dtype=">f8", count=3)
            dx = meta_inf[0]
            dy = meta_inf[1]
            dz = meta_inf[2]

            meta_inf = np.fromfile(f, dtype=">i4", count=1)
            nsubgrid = meta_inf[0]

            data = np.ndarray(shape=(nz, ny, nx), dtype='>f8')
            #data = np.ndarray(shape=(nz, ny, nx), dtype=dtype)

            for s in range(nsubgrid):
                meta_inf = np.fromfile(f, dtype=">i4", count=9)
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

                data[iz : iz + nz, iy : iy + ny, ix : ix + nx] = np.fromfile(
                    f, dtype=">f8", count=nn
                ).reshape((nz, ny, nx))

            data = data[None, :]  # expand by empty dimension

            size = ht.torch.tensor(len(data.shape), dtype=ht.torch.int32)
            comm.Bcast(size)
            shape = ht.torch.tensor(data.shape, dtype=ht.torch.int32)
            comm.Bcast(shape)
            shape = shape.numpy()
    else:  # non-root processes
        size = ht.torch.empty(1, dtype=ht.torch.int32)
        comm.Bcast(size)
        shape = ht.torch.empty(size, dtype=ht.torch.int32)
        comm.Bcast(shape)
        shape = shape.numpy()
        shape[0] = 0  # empty dimension
        data = np.empty(shape, dtype=dtype)

    split = ht.sanitize_axis(shape[1:], split)
    if split is not None:  # account for added empty dimension
        split = split + 1
#    try:
#    return (
#            ht.array(data, is_split=0, comm=comm).resplit_(split).squeeze(0)
#        )  # reduce by empty dimension and split data
#    except ValueError:  # If byteorder is non native because HeAt needs native byteorder
    data = data.astype(data.dtype.newbyteorder("="))
    return (
            ht.array(data, is_split=0, comm=comm).resplit_(split).squeeze(0)
        )  # reduce by empty dimension and split data


def read_nc4(dict, dir=".", split=None):
    """    dict: files -> variables of that file;
             multiple variables can be given in a list
        TODO: - Allow different splitAxes
              - Add support for same Variable for every file
              - Introduce names, using filenames is not feasible for users
              - Combine names with Option to flatten output:
                        name -> dndarray
                (might use named tensors for it (torch 1.4))

    Parameters
    ----------
    dict : Dictionary
        files -> variables of that file;
                 multiple variables can be given in a list.
                 filenames and variable-names as strings
    dir : str
        Directory of the files (the default is '.').
    split : type
        split-Axis used by HeAT (the default is None).
        If the specified split-Axis is not usable, it is set to None.

    Returns
    -------
    nested dictionary
        file -> variable -> dndarray
        access data using: data[file][variable]

    Raises
    -------
    ExceptionName
        Why the exception is raised.

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
            if not file.startswith("/") or dir != ".":
                path = dir + "/" + file
            try:
                return_dict[file][var] = ht.load_netcdf(path, split=split, variable=var)
            except ValueError:
                return_dict[file][var] = ht.load_netcdf(path, split=None, variable=var)
    return return_dict


def create_nc4_Construct(
    filename, tstart=0, loni=0, lati=0, delta=1.0, ntimes=None, nlevels=None, nlons=None, nlats=None
):
    """Generates, if given, the ('time', 'level', 'lon', 'lat')-Dimensions of a NetCDF4-File. The time-dimension is set unlimited.

    Parameters
    ----------
    filename : str
        Filename.
    tstart : int
        Offset in Time
    loni : int
        Offset in Longitude
    lati : int
        Offset in Latitude
    delta : float
        Stepsize in both Lat and Lon Direction
    ntimes : int
        Size of time-dimension
    nlevels : int
        Size of level-dimension
    nlons : int
        Size of lon-dimension
    nlats : int
        Size of lat-dimension



    Returns
    -------
    None
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    mode = "w"
    if ntimes is not None:
        time_out = tstart + ht.arange(ntimes, dtype=ht.float32, split=0)
        ht.save_netcdf(time_out, filename, "time", mode, ["time"], is_unlimited=True)
        mode = "r+"
    if nlevels is not None:
        lev_out = 1.0 + ht.arange(nlevels, dtype=ht.float32, split=0)
        ht.save_netcdf(lev_out, filename, "level", mode, ["level"])
        mode = "r+"
    if nlons is not None:
        lons_out = loni + delta * ht.arange(nlons, dtype=ht.float32, split=0)
        ht.save_netcdf(lons_out, filename, "lon", mode, ["lon"])
        mode = "r+"
    if nlats is not None:
        lats_out = lati + delta * ht.arange(nlats, dtype=ht.float32, split=0)
        ht.save_netcdf(lats_out, filename, "lat", mode, ["lat"])
        mode = "r+"


def write_nc4(dict, filename, missval=None, slices=[slice(None)]):
    """
    dic : {varName : dndarray}
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
    dims3D = ("time", "level", "lon", "lat")
    dims2D = ("time", "lon", "lat")
    mode = "r+"
    if not os.path.isfile(filename):
        mode = "w"
    for varName, data in dict.items():
        if len(data.shape) == 4:  # 4 because time is included
            dim_names = dims3D
        else:
            dim_names = dims2D
        ht.save_netcdf(
            data,
            filename,
            varName,
            mode,
            dim_names,
            file_slices=slices,
            fill_value=missval,
            zlib=False,
            least_significant_digit=6,
        )  # zlib should be True
        mode = "r+"
