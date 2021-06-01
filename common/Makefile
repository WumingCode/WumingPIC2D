# -*- Makefile -*-
include ../compiler.mk

SRCS   = wuming2d.f90 boundary.f90 field.f90 particle.f90 mom_calc.f90 sort.f90 mpi_set.f90 fio.f90 h5io.f90
OBJS   = $(SRCS:%.f90=%.o)
TARGET = lib$(WUMING_LIB_COMMON).a

$(TARGET) : $(OBJS)
	$(AR) r $@ $(OBJS)

clean :
	rm -f $(TARGET) $(OBJS) *.mod *.out

# Dependencies
field.o: boundary.o
h5io.o: mpi_set.o
wuming2d.o: boundary.o field.o particle.o mom_calc.o sort.o mpi_set.o fio.o h5io.o
