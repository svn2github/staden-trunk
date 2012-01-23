#
# Makefile for mutscan
#

PROGS = mutscan

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(IOLIB_INC) $(MUT_INC)


HLIBS=\
	$(MUT_LIB) \
	$(IOLIB_LIB) \
	$(TK_LIB) \
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB)

OBJSD	= \
	main.o \
	pathutil.o \
	stringlist.o


mutscan: $(OBJSD)
	$(CXXLD) $(CXXLDFLAGS) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(HLIBS) $(LIBSC)


distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.[ch]pp $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

DEPEND_OBJ = $(OBJS) $(OBJS) $(OBJSA)

# DO NOT DELETE THIS LINE -- make depend depends on it
