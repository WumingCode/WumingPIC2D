# -*- Makefile -*-

# HDF5
#HDF5DIR = $(shell spack location -i hdf5 %gcc)

# compilers and arguments
AR      = ar
CC      = mpicc -fopenmp
FC      = mpif90 -fopenmp
FCFLAGS = -cpp -I$(BASEDIR)/include
LDFLAGS = -L$(BASEDIR)/lib

