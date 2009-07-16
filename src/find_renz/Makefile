#
# Makefile for find_renz
#

PROGS = $(O)/find_renz

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(SEQUTILS_INC) $(MISC_INC)

FRSRC=$(O)

OBJS	= \
	$(FRSRC)/find_renz.o

LIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB)

$(O)/find_renz: $(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(LIBS) $(LIBSC) 


DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
