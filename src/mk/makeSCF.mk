#
# Makefile for makeSCF
#

PROGS = $(O)/makeSCF

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC)

MAKESCFSRC=$(O)

OBJS	= $(MAKESCFSRC)/makeSCF.o

$(O)/makeSCF: $(OBJS)
	$(CLD) -o $@ $(OBJS) $(IOLIB_LIB) $(LIBSC) $(IOLIB_DEP)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
