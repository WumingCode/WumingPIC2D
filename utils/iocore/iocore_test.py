#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess
import numpy as np
import h5py

def create_hdf5(h5file, rawfile):
    # open and read small data
    with open(rawfile, 'rb') as fp:
        nproc  = np.fromfile(fp, 'i4', 1)[0]
        dims   = np.fromfile(fp, 'i4', 2)
        nx     = np.fromfile(fp, 'i4', 1)[0]
        ny     = np.fromfile(fp, 'i4', 1)[0]
        var_i4 = np.fromfile(fp, 'i4', 1)[0]
        var_i8 = np.fromfile(fp, 'i8', 1)[0]
        var_r4 = np.fromfile(fp, 'f4', 1)[0]
        var_r8 = np.fromfile(fp, 'f8', 1)[0]
        ldesc  = np.fromfile(fp, 'i4', 1)[0]
        desc   = fp.read(ldesc).decode('utf-8')
        offset = fp.tell()
        print('{:30s} : {:4d}'.format('# MPI process', nproc))
        print('{:30s} : {:4d}'.format('# MPI process in x', dims[0]))
        print('{:30s} : {:4d}'.format('# MPI process in y', dims[1]))
        print('{:30s} : {:4d}'.format('# local grid in x', nx))
        print('{:30s} : {:4d}'.format('# local grid in y', nx))
        print('{:30s} : {:4d}'.format('# global grid in x', nx*dims[0]))
        print('{:30s} : {:4d}'.format('# global grid in y', nx*dims[1]))
        print('{:30s} : {:4d}'.format('# integer(4)', var_i4))
        print('{:30s} : {:4d}'.format('# integer(8)', var_i8))
        print('{:30s} : {:10.4f}'.format('# real(4)', var_r4))
        print('{:30s} : {:10.4f}'.format('# real(8)', var_r8))
        print('{:30s} : {:4d}'.format('# character(*) len', ldesc))
        print('{:30s} : {:s}'.format('# character(*)', desc))

    # create hdf5 and write dataset for large data
    with h5py.File(h5file, 'w') as fp:
        n = nx*ny * nproc
        s = (ny*dims[1], nx*dims[0],)

        # int32
        print('{:30s} : {:8d}'.format('# offset for integer(4)', offset))
        size = 4 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_i4', shape=s, dtype='i4', external=ext)
        offset = offset + size

        # int64
        print('{:30s} : {:8d}'.format('# offset for integer(8)', offset))
        size = 8 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_i8', shape=s, dtype='i8', external=ext)
        offset = offset + size

        # float32
        print('{:30s} : {:8d}'.format('# offset for real(4)', offset))
        size = 4 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_r4', shape=s, dtype='f4', external=ext)
        offset = offset + size

        # int64
        print('{:30s} : {:8d}'.format('# offset for real(8)', offset))
        size = 8 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_r8', shape=s, dtype='f8', external=ext)
        offset = offset + size


def h5dump(h5file):
    command = ["h5dump", h5file]
    result = subprocess.run(command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    print('\n\n*** {} ***'.format('result of h5dump start here'))
    print(result.stdout.decode('utf-8'))
    print('*** {} ***\n\n'.format('result of h5dump end here'))


if __name__ == '__main__':
    h5file  = 'iocore_test.h5'
    rawfile = 'iocore_test.raw'
    create_hdf5(h5file, rawfile)
    h5dump(h5file)
