.PHONY : all
.PHONY : clean

.SUFFIXES :
.SUFFIXES : .o .f90

.f90.o:
	$(FC) -c $(FFLAGS) $< -o $@
