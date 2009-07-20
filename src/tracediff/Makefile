#
# Makefile for tracediff
#

PROGS = tracediff

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk


INCLUDES_E += $(IOLIB_INC) $(MUT_INC)


HLIBS=\
	$(TEXTUTILS_LIB) \
	$(MUT_LIB) \
	$(IOLIB_LIB)

OBJSD	= main.o


tracediff: $(OBJSD)
	$(CXXLD) $(CXXLDFLAGS) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(HLIBS) $(LIBSC)


pathutil.cpp:
	ln -s ../mutlib/pathutil.cpp .

pathutil.h:
	ln -s ../mutlib/pathutil.h .

stringlist.cpp:
	ln -s ../mutlib/stringlist.cpp .

stringlist.hpp:
	ln -s ../mutlib/stringlist.hpp .


distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.[ch]pp $(S)/Makefile $(S)/dependencies $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

DEPEND_OBJ = $(OBJSD)

include dependencies




