#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import io
import subprocess
import difflib
import numpy as np
import h5py

MPIEXEC = 'mpiexec'
TEST1   = 'iocore_test1'
TEST2   = 'iocore_test2'
hdffile = 'iocore_test.h5'
rawfile = 'iocore_test.raw'


def unittest():
    # execute MPI program
    PIPE = subprocess.PIPE
    test1 = [MPIEXEC, '-n', '8', './' + TEST1]
    test2 = [MPIEXEC, '-n', '8', './' + TEST2]
    result1 = subprocess.run(test1, stdout=PIPE, stderr=PIPE)
    result2 = subprocess.run(test2, stdout=PIPE, stderr=PIPE)
    stdout1 = result1.stdout.decode('utf-8')
    stdout2 = result2.stdout.decode('utf-8')
    stderr1 = result1.stderr.decode('utf-8')
    stderr2 = result2.stderr.decode('utf-8')

    # check stderr
    if not (stderr1 == '' and stderr2 == ''):
        if not stderr1 == '':
            print('*** Error detected in executing {} ***'.format(TEST1))
            print(stderr1)
            raise RuntimeError
        if not stderr2 == '':
            print('*** Error detected in executing {} ***'.format(TEST2))
            print(stderr2)
            raise RuntimeError

    # read in python
    with io.StringIO() as buf:
        sys.stdout = buf
        create_hdf5(hdffile, rawfile)
        stdout0 = buf.getvalue()
    sys.stdout = sys.__stdout__

    # check output
    diff01 = [diff for diff in difflib.unified_diff(stdout0, stdout1)]
    diff12 = [diff for diff in difflib.unified_diff(stdout1, stdout2)]
    if len(diff01) != 0:
        print('Error: python output and {} output differ!'.format(TEST1))
        print('*** begin diff ***')
        sys.stdout.writelines(diff01)
        print('\n*** end diff ***')
        raise RuntimeError
    if len(diff12) != 0:
        print('Error: {} output and {} output differ!'.format(TEST1, TEST2))
        print('*** begin diff ***')
        sys.stdout.writelines(diff12)
        print('\n*** end diff ***')
        raise RuntimeError

    print('\n***')
    print('iocore test seems to be succesful !')
    print('HDF5 file {} has been generated'.format(hdffile))
    print('***\n')


def create_hdf5(hdffile, rawfile):
    # open and read small data
    with open(rawfile, 'rb') as fp:
        endian = np.fromfile(fp, 'i4', 1)[0]
        if endian == 1:
            bo = '<'
        elif endian == 16777216:
            bo = '>'
        else:
            errmsg = 'Error: unrecognized endian flag: {}'.format(endian)
            raise RuntimeError(errmsg)
        nproc  = np.fromfile(fp, bo+'i4', 1)[0]
        dims   = np.fromfile(fp, bo+'i4', 2)
        nx     = np.fromfile(fp, bo+'i4', 1)[0]
        ny     = np.fromfile(fp, bo+'i4', 1)[0]
        var_i4 = np.fromfile(fp, bo+'i4', 1)[0]
        var_i8 = np.fromfile(fp, bo+'i8', 1)[0]
        var_r4 = np.fromfile(fp, bo+'f4', 1)[0]
        var_r8 = np.fromfile(fp, bo+'f8', 1)[0]
        ldesc  = np.fromfile(fp, bo+'i4', 1)[0]
        desc   = fp.read(ldesc).decode('utf-8')
        offset = fp.tell()
        print('{:30s} : {:4d}'.format('# endian flag', endian))
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
    with h5py.File(hdffile, 'w') as fp:
        n = nx*ny * nproc
        s = (ny*dims[1], nx*dims[0],)

        # int32
        print('{:30s} : {:8d}'.format('# offset for integer(4)', offset))
        size = 4 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_i4', shape=s, dtype=bo+'i4', external=ext)
        offset = offset + size

        # int64
        print('{:30s} : {:8d}'.format('# offset for integer(8)', offset))
        size = 8 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_i8', shape=s, dtype=bo+'i8', external=ext)
        offset = offset + size

        # float32
        print('{:30s} : {:8d}'.format('# offset for real(4)', offset))
        size = 4 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_r4', shape=s, dtype=bo+'f4', external=ext)
        offset = offset + size

        # int64
        print('{:30s} : {:8d}'.format('# offset for real(8)', offset))
        size = 8 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_r8', shape=s, dtype=bo+'f8', external=ext)
        offset = offset + size


if __name__ == '__main__':
    unittest()
