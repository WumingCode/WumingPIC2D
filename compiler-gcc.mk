# -*- Makefile -*-

# base directory
BASEDIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

# include common configuration
include $(BASEDIR)/common.mk

# wuming
WUMING_LIB_COMMON = wuming2d_common

# HDF5
HDF5DIR = $(shell spack location -i hdf5 %gcc)

# compilers and arguments
AR     = ar
CC     = mpicc -fopenmp
FC     = mpif90 -fopenmp
F90    = mpif90 -fopenmp
FFLAGS = -I$(BASEDIR)/common -I$(HDF5DIR)/include
LFLAGS = -L$(BASEDIR)/common -L$(HDF5DIR)/lib -l$(WUMING_LIB_COMMON) -lhdf5 -lhdf5_fortran
