#
# Makefile for screen_seq
#

PROGS = $(O)/screen_seq

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(SEQLIB_INC) $(TEXTUTILS_INC) $(IOLIB_INC) $(SEQUTILS_INC) $(MISC_INC)

#EXTRA_LIBS += -ldmalloc

DEPS=\
	$(IOLIB_DEP) \
	$(SEQUTILS_DEP)


CLPSRC=$(O)

OBJS	= \
	$(CLPSRC)/screen_seq.o

HLIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(IOUTILS_LIB) \
	$(MISC_LIB) \
	$(IOLIB_LIB) \
	$(MATH_LIB)

#30/6/99 johnt -added abstracttions to support Visual C++ linker
$(O)/screen_seq: $(OBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(HLIBS) $(DEPS) $(LIBSC)


DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
