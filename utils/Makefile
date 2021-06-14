# -*- Makefile -*-
include ../common.mk

SUBDIRS = json iocore

default:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done

clean :
	rm -f $(OBJS) $(TARGET) *.mod *.out
	# clean subdirectories
	for dir in $(SUBDIRS); do \
		$(MAKE) clean -C $$dir; \
	done

json: iocore
