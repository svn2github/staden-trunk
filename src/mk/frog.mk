#
# Makefile for frog
#

PROGS = $(O)/frog $(O)/toad

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC)

FROGBIN=$(O)

FROGOBJ =\
	$(FROGBIN)/frog.o\
	$(FROGBIN)/mach-io.o

$(O)/frog : $(FROGOBJ)
	$(CLD) -o $@ $(FROGOBJ) $(LIBSC)


DEPS=\
	$(IOLIBBIN)/libread.a \
	$(IOLIBBIN)/libio-utils.a \
	$(MISCBIN)/libmisc.a

LIBS=\
	$(IOLIB_LIB) \
	$(MISC_LIB)

TOADOBJ =\
	$(FROGBIN)/toad.o\
	$(FROGBIN)/mach-io.o

$(O)/toad : $(TOADOBJ)
	$(CLD) -o $@ $(TOADOBJ) $(LIBS) $(LIBSC)

DEPEND_OBJ = $(FROGOBJ) $(TOADOBJ)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
