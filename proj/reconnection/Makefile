# -*- Makefile -*-
include ../../common.mk

SRCS   = main.f90 app.f90 boundary_reconnection.f90
OBJS   = $(SRCS:%.f90=%.o)
TARGET = main.out

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS) -lwuming2d_common -lwuming_utils $(HDF5LIB)

clean :
	rm -f $(TARGET) $(OBJS) *.i *.mod *.out

# Dependencies
main.o: app.o
app.o: boundary_reconnection.o
