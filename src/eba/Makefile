#
# Makefile for eba (Estimate Base Accuracy)
#

PROGS = eba

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC)

EBAOBJS	= qual.o conf.o

eba: $(EBAOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(EBAOBJS) $(IOLIB_LIB) $(MISC_LIB) $(LIBSC) $(LIBM) $(MATH_LIB)

DEPEND_OBJ = $(EBAOBJS) $(SNROBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it.

conf.o: $(SRCROOT)/eba/conf.h
qual.o: $(SRCROOT)/eba/conf.h
