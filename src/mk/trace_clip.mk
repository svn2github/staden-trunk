#
# Makefile for trace_clip and scale_trace_clip
#

PROGS = $(O)/trace_clip $(O)/scale_trace_clip

SRCROOT=..

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC)

CLPSRC=$(O)

OBJS	= $(CLPSRC)/trace_clip.o

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(O)/trace_clip: $(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(IOLIB_LIB) $(MISC_LIB) $(LIBSC)

OBJSS	= $(CLPSRC)/scale_trace_clip.o

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(O)/scale_trace_clip: $(OBJSS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSS) $(IOLIB_LIB) $(MISC_LIB) $(LIBSC)

DEPEND_OBJ = $(OBJS) $(OBJSS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
