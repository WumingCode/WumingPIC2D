# -*- Makefile -*-

# base directory
BASEDIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

# include common configuration
include $(BASEDIR)/common.mk

# HDF5
HDF5DIR = $(shell spack location -i hdf5 %gcc)

# compilers and arguments
AR     = ar
CC     = mpicc -fopenmp
FC     = mpif90 -fopenmp
F90    = mpif90 -fopenmp
FFLAGS = -cpp -I$(BASEDIR)/common -I$(BASEDIR)/include -I$(HDF5DIR)/include
LFLAGS = -L$(BASEDIR)/lib -L$(HDF5DIR)/lib -l$(WUMING_LIB_COMMON) -lwuming_utils -lhdf5 -lhdf5_fortran
