#
# Makefile for trace_diff
#

PROGS = $(O)/trace_diff

#SRCROOT=$(STADENROOT)/src
SRCROOT=..

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

TRACEDIFFBIN=$(O)

# next 4 lines for static linking
#L = $(SRCROOT)/lib-static/$(O)
#SHLIB_SUFFIX=.a
#SHLIB_LD=ar
#SHLIB_LDFLAGS=rv

INCLUDES_E += $(SEQLIB_INC) $(IOLIB_INC) $(SEQUTILS_INC) $(MISC_INC) $(TEXTUTILS_INC)


HLIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(IOUTILS_LIB) \
	$(MISC_LIB) \
	$(IOLIB_LIB) \
	$(MATH_LIB)



OBJSD	= \
	$(TRACEDIFFBIN)/main.o \
	$(TRACEDIFFBIN)/rsalign.o \
	$(TRACEDIFFBIN)/diff_trace.o \
	$(TRACEDIFFBIN)/find_mutations.o

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(O)/trace_diff: $(OBJSD)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(HLIBS) $(LIBSC)

install:
	cp $(PROGS) $(INSTALLBIN)

DEPEND_OBJ = $(OBJS) $(OBJS) $(OBJSA)

include dependencies




