#
# Makefile for eba (Estimate Base Accuracy)
#

PROGS = $(O)/eba

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC)

EBASRC=$(O)

OBJS	= $(EBASRC)/qual.o $(EBASRC)/conf.o

$(O)/eba: $(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJS) $(IOLIB_LIB) $(LIBSC) $(IOLIB_DEP)

DEPEND_OBJ = $(OBJS)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
