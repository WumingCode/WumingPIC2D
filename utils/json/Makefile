# -*- Makefile -*-
include ../../compiler.mk

SRCS   = \
	json_file_module.F90 \
	json_module.F90 \
	json_string_utilities.F90 \
	json_kinds.F90 \
	json_parameters.F90 \
	json_value_module.F90
OBJS   = $(SRCS:%.F90=%.o)
MODS   = $(SRCS:%.F90=%.mod)

default: $(OBJS)
	$(AR) r $(WUMING_LIB)/libwuming_utils.a $(OBJS)
	cp $(MODS) $(WUMING_INCLUDE)

clean :
	rm -f $(OBJS) *.o *.mod *.out

# dependencies: modified from output of 'gfortran -cpp -M'
json_file_module.o: json_file_module.F90 \
 json_macros.inc json_initialize_dummy_arguments.inc \
 json_initialize_arguments.inc json_kinds.o json_parameters.o \
 json_string_utilities.o json_value_module.o
json_kinds.o: json_kinds.F90
json_module.o: json_module.F90 json_kinds.o \
 json_parameters.o json_value_module.o json_file_module.o
json_parameters.o: json_parameters.F90 json_kinds.o
json_string_utilities.o: \
 json_string_utilities.F90 \
 json_kinds.o json_parameters.o
json_value_module.o: json_value_module.F90 \
 json_macros.inc json_initialize_dummy_arguments.inc \
 json_initialize_arguments.inc json_get_scalar_by_path.inc \
 json_get_vec_by_path.inc json_get_vec_by_path_alloc.inc \
 json_kinds.o json_parameters.o json_string_utilities.o
