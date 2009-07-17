#
# Makefile for mutscan
#

PROGS = mutscan

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC) $(MUT_INC)


HLIBS=\
	$(TEXTUTILS_LIB) \
	$(MUT_LIB) \
	$(IOLIB_LIB) \
	$(TK_LIB)

OBJSD	= \
	main.o \
	pathutil.o \
	stringlist.o


mutscan: $(OBJSD)
	$(CXXLD) $(CXXLDFLAGS) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(HLIBS) $(LIBSC)


distsrc: distsrc_dirs
	-cp -R *.[ch] *.[ch]pp Makefile dependencies $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

DEPEND_OBJ = $(OBJS) $(OBJS) $(OBJSA)

include dependencies




