EXEDIR = ./
HDF5DIR = /vol0004/apps/oss/spack-v0.16/opt/spack/linux-rhel8-a64fx/fj-4.3.1/hdf5-1.10.7-7cxah5foptv7l7osleavke2kozq3y3g7/
FC = mpifrtpx
FFLAGS = -Kfast,openmp -I$(HDF5DIR)/include
LFLAGS = -L$(HDF5DIR)/lib/ -lhdf5 -lhdf5_fortran

OBJS = fio.o particle.o field.o boundary.o mpi_set.o const.o init.o main.o sort.o mom_calc.o h5io.o

.PHONY : all
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) -c $< $(FFLAGS)

all:
	@echo 'make test or make 1 or make 2 or make 3'

test: $(OBJS)
	$(FC) -o $(EXEDIR)em2dsk_$@.out $(FFLAGS) $(LFLAGS) $(OBJS)

1: $(OBJS)
	$(FC) -o $(EXEDIR)em2dsk$@.out $(FFLAGS) $(LFLAGS) $(OBJS)

2: $(OBJS)
	$(FC) -o $(EXEDIR)em2dsk$@.out $(FFLAGS) $(LFLAGS) $(OBJS)

3: $(OBJS)
	$(FC) -o $(EXEDIR)em2dsk$@.out $(FFLAGS) $(LFLAGS) $(OBJS)

# Dependencies
field.o : boundary.o
main.o  : init.o const.o mpi_set.o boundary.o fio.o particle.o field.o sort.o mom_calc.o h5io.o
init.o  : const.o mpi_set.o boundary.o fio.o sort.o mom_calc.o h5io.o
h5io.o  : mpi_set.o

clean :
	rm -f $(OBJS) $(TARGET) *.mod *.out
