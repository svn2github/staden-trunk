#
# Makefile for eba (Estimate Base Accuracy)
#

PROGS = $(O)/eba

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC)

EBASRC=$(O)

EBAOBJS	= $(EBASRC)/qual.o $(EBASRC)/conf.o

$(O)/eba: $(EBAOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(EBAOBJS) $(IOLIB_LIB) $(MISC_LIB) $(LIBSC) $(LIBM) $(MATH_LIB)

DEPEND_OBJ = $(EBAOBJS) $(SNROBJS)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
