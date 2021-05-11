EXEDIR = ./
HDF5DIR = /usr/local/Cellar/hdf5-mpi/1.12.0_1
FC = mpif90
FFLAGS = -O3 -fopenmp -I$(HDF5DIR)/include/ -fallow-argument-mismatch -fbacktrace -fbounds-check
#-fallow-argument-mismatch: to avoid bugs caused by gcc 10
LFLAGS = -L$(HDF5DIR)/lib/ -lhdf5 -lhdf5_fortran

OBJS = fio.o particle.o field.o boundary.o mpi_set.o const.o init.o main.o sort.o mom_calc.o h5util.o h5io.o

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
main.o : init.o const.o mpi_set.o boundary.o fio.o particle.o field.o sort.o mom_calc.o h5io.o h5util.o
init.o : const.o mpi_set.o boundary.o fio.o sort.o mom_calc.o h5io.o
h5io.o : h5util.o mpi_set.o

clean :
	rm -f $(OBJS) $(TARGET) *.mod *.out
