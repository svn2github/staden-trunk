#
# Makefile for tracediff
#

PROGS = $(O)/tracediff

#SRCROOT=$(STADENROOT)/src
SRCROOT=..

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

TRACEDIFFBIN=$(O)


INCLUDES_E += $(IOLIB_INC) $(MUT_INC)


HLIBS=\
	$(TEXTUTILS_LIB) \
	$(MUT_LIB) \
	$(IOLIB_LIB)

OBJSD	= \
	$(TRACEDIFFBIN)/main.o \
	$(TRACEDIFFBIN)/pathutil.o \
	$(TRACEDIFFBIN)/stringlist.o


$(O)/tracediff: $(OBJSD)
	$(CXXLD) $(CXXLDFLAGS) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(HLIBS) $(LIBSC)


pathutil.cpp:
	ln -s ../mutlib/pathutil.cpp

pathutil.h:
	ln -s ../mutlib/pathutil.h

stringlist.cpp:
	ln -s ../mutlib/stringlist.cpp

stringlist.hpp:
	ln -s ../mutlib/stringlist.hpp


distsrc: distsrc_dirs
	-cp -R *.[ch] *.[ch]pp Makefile dependencies $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

DEPEND_OBJ = $(OBJS) $(OBJS) $(OBJSA)

include dependencies




