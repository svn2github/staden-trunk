#
# Makefile for lookup program
#
PROGS=$(O)/pregap_lookup
SRCROOT=..

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

OBJ = $(O)/lookup.o

$(O)/pregap_lookup: $(OBJ)
	$(CLD) -o $@ $(OBJ) $(LIBSC)

DEPEND_OBJ = $(OBJ)

install:
	cp $(PROGS) $(INSTALLBIN)
