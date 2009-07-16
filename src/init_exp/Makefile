#
# Makefile for init_exp
#

PROGS = $(O)/init_exp

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(IOLIB_INC)

INITEXPSRC=$(O)

OBJS=$(INITEXPSRC)/init_exp.o

$(O)/init_exp:	$(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(IOLIB_LIB) $(MISC_LIB) $(LIBSC)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
