.PHONY : all
.PHONY : clean

.SUFFIXES :
.SUFFIXES : .o .f90

%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

%.o: %.F90
	$(FC) -c $(FFLAGS) $< -o $@


# Wuming
WUMING_INCLUDE    = $(BASEDIR)/include
WUMING_LIB        = $(BASEDIR)/lib
WUMING_LIB_COMMON = wuming2d_common
