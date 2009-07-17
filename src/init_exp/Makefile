#
# Makefile for init_exp
#

PROGS = init_exp

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(IOLIB_INC)

OBJS=init_exp.o

init_exp:	$(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(IOLIB_LIB) $(MISC_LIB) $(LIBSC)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
