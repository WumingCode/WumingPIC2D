# -*- Makefile -*-
include ../../compiler.mk

SRCS   = iocore.f90
OBJS   = $(SRCS:%.f90=%.o)
MODS   = $(SRCS:%.f90=%.mod)
EXES   = iocore_test1 iocore_test2

default: $(OBJS)
	$(AR) r $(WUMING_LIB)/libwuming_utils.a $(OBJS)
	cp $(MODS) $(WUMING_INCLUDE)

test: $(EXES)
	python iocore_test.py

iocore_test1: iocore_test1.o iocore_test.o iocore.o
	$(FC) -o $@ $^

iocore_test2: iocore_test2.o iocore_test.o iocore.o
	$(FC) -o $@ $^

clean :
	rm -f $(EXES) $(OBJS) *.o *.mod *.out iocore_test.raw iocore_test.h5

# Dependencies
iocore_test1.o: iocore.o iocore_test.o
iocore_test2.o: iocore.o iocore_test.o
