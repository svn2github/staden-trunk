#
# Makefile for alf
#

PROGS = $(O)/alfsplit

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

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
