# -*- Makefile -*-

# base directory
BASEDIR := $(dir $(lastword $(MAKEFILE_LIST)))

# include common configuration
include $(BASEDIR)common.mk

# HDF5
HDF5DIR = $(shell spack location -i hdf5 %intel)

# compilers and arguments
AR     = ar
CC     = mpiicc -qopenmp
FC     = mpiifort -qopenmp
F90    = mpiifort -qopenmp
FFLAGS = -I$(BASEDIR)core  -I$(HDF5DIR)/include
LFLAGS = -L$(BASEDIR)core -L$(HDF5DIR)/lib -lwumingcore2d -lhdf5 -lhdf5_fortran
