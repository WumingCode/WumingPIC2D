# -*- Makefile -*-
include compiler.mk

SUBDIRS = common

default:
	$(MAKE) -C common

clean :
	rm -f $(OBJS) $(TARGET) *.mod *.out
	# clean subdirectories
	for dir in $(SUBDIRS); do \
		$(MAKE) clean -C $$dir; \
	done
