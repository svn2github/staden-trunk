#
# Makefile for alf
#

PROGS = $(O)/alfsplit

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC)

ALFBIN=$(O)

#
# Make the programs
#
ALFOBJ = $(ALFBIN)/alfsplit.o

$(O)/alfsplit: $(ALFOBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(ALFOBJ) $(LIBSC)

RSOBJ = $(ALFBIN)/2rs.o

$(O)/2rs: $(RSOBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(RSOBJ) $(LIBSC)

DEPEND_OBJ = $(ALFOBJ) $(RSOBJ)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
