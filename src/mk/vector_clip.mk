#
# Makefile for vector_clip
#

PROGS = $(O)/vector_clip

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(SEQLIB_INC) $(IOLIB_INC) $(SEQUTILS_INC) $(MISC_INC) $(TEXTUTILS_INC)

#EXTRA_LIBS += -ldmalloc

DEPS=\
	$(IOLIB_DEP) \
	$(SEQUTILS_DEP)


CLPSRC=$(O)

OBJS	= \
	$(CLPSRC)/vector_clip.o

HLIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(IOUTILS_LIB) \
	$(MISC_LIB) \
	$(IOLIB_LIB) \
	$(MATH_LIB)

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(O)/vector_clip: $(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(HLIBS) $(LIBSC) 


DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
