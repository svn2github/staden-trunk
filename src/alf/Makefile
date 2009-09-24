#
# Makefile for alf
#

PROGS = alfsplit

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC)

#
# Make the programs
#
ALFOBJ = alfsplit.o

alfsplit: $(ALFOBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(ALFOBJ) $(LIBSC)

RSOBJ = $(ALFBIN)/2rs.o

2rs: $(RSOBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(RSOBJ) $(LIBSC)

DEPEND_OBJ = $(ALFOBJ) $(RSOBJ)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it.

alfsplit.o: $(SRCROOT)/Misc/os.h
alfsplit.o: $(PWD)/staden_config.h
