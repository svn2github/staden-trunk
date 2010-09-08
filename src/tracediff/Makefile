#
# Makefile for tracediff
#

PROGS = tracediff

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk


INCLUDES_E += $(MISC_INC) $(IOLIB_INC) $(MUT_INC)


HLIBS=\
	$(TEXTUTILS_LIB) \
	$(MUT_LIB) \
	$(IOLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(MISC_LIB)	

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
	cp $(S)/*.[ch] $(S)/*.[ch]pp $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

DEPEND_OBJ = $(OBJSD)

# DO NOT DELETE THIS LINE -- make depend depends on it.

main.o: $(SRCROOT)/Misc/os.h
main.o: $(PWD)/staden_config.h
main.o: $(SRCROOT)/mutlib/mutlib.h
main.o: $(SRCROOT)/mutlib/pathutil.h
main.o: $(SRCROOT)/mutlib/stringlist.hpp
main.o: $(SRCROOT)/tracediff/staden.h
