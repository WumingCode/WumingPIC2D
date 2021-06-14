# -*- Makefile -*-
include ../../compiler.mk

SRCS   = main.f90 app.f90
OBJS   = $(SRCS:%.f90=%.o)
TARGET = main.out

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(LFLAGS)

clean :
	rm -f $(TARGET) $(OBJS) *.mod *.out

# Dependencies
main.o: app.o
