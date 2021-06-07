# -*- Makefile -*-
include ../../compiler.mk

SRCS   = iocore.f90
OBJS   = $(SRCS:%.f90=%.o)
MODS   = $(SRCS:%.f90=%.mod)
EXES   = jsonio_test iocore_test1 iocore_test2

default: $(OBJS)
	$(AR) r $(WUMING_LIB)/libwuming_utils.a $(OBJS)
	cp $(MODS) $(WUMING_INCLUDE)

test: $(EXES)
	python iocore_test.py

iocore_test1: iocore_test1.o iocore_test.o iocore.o
	$(FC) -o $@ $^

iocore_test2: iocore_test2.o iocore_test.o iocore.o
	$(FC) -o $@ $^

jsonio_test: jsonio.o jsonio_test.o
	$(FC) -o $@ $^ $(LFLAGS)

clean :
	rm -f $(EXES) $(OBJS) *.o *.mod *.out iocore_test.raw iocore_test.h5

# Dependencies
jsonio_test.o: jsonio.o
iocore_test1.o: iocore.o iocore_test.o
iocore_test2.o: iocore.o iocore_test.o
