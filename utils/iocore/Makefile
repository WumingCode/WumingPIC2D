# -*- Makefile -*-
include ../../compiler.mk

SRCS   = iocore.f90
OBJS   = $(SRCS:%.f90=%.o)
EXES   = iocore_test1 iocore_test2

default: test

iocore_test1: iocore_test1.o iocore_test.o iocore.o
	$(FC) -o $@ $^

iocore_test2: iocore_test2.o iocore_test.o iocore.o
	$(FC) -o $@ $^

test: $(EXES)
	python iocore_test.py

clean :
	rm -f $(EXES) $(OBJS) *.o *.mod *.out

# Dependencies
iocore_test1.o: iocore.o iocore_test.o
iocore_test2.o: iocore.o iocore_test.o
