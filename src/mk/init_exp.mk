#
# Makefile for init_exp
#

PROGS = $(O)/init_exp

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC)

INITEXPSRC=$(O)

OBJS=$(INITEXPSRC)/init_exp.o

# 7/1/99 johnt - abstractions to support Visual C++
$(O)/init_exp:	$(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(IOLIB_LIB) $(LIBSC) $(IOLIB_DEP)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
