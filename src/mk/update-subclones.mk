#
# Makefile for subclonedb
#

PROGS = $(O)/update_subclones

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

UPDATE_SUBCLONESBIN=$(O)

USOBJ	=\
	$(UPDATE_SUBCLONESBIN)/update_subclones.o

$(O)/update_subclones: $(USOBJ)
	$(CLD) -o $@ $(USOBJ) $(MISC_LIB) $(LIBSC)

DEPEND_OBJ = $(USOBJ)


install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
