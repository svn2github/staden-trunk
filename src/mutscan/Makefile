#
# Makefile for mutscan
#

PROGS = $(O)/mutscan

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

MUTSCANBIN=$(O)


INCLUDES_E += $(IOLIB_INC) $(MUT_INC)


HLIBS=\
	$(TEXTUTILS_LIB) \
	$(MUT_LIB) \
	$(IOLIB_LIB) \
	$(TK_LIB)

OBJSD	= \
	$(MUTSCANBIN)/main.o \
	$(MUTSCANBIN)/pathutil.o \
	$(MUTSCANBIN)/stringlist.o


$(O)/mutscan: $(OBJSD)
	$(CXXLD) $(CXXLDFLAGS) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(HLIBS) $(LIBSC)


distsrc: distsrc_dirs
	-cp -R *.[ch] *.[ch]pp Makefile dependencies $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

DEPEND_OBJ = $(OBJS) $(OBJS) $(OBJSA)

include dependencies




