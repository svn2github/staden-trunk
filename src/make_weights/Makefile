#
# Makefile for make_weights
#

PROGS = $(O)/make_weights

SRCROOT=..
#SRCROOT=/home5/pubseq/share/src
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk
#COPTDEBUG	= -O2


INCLUDES_E += $(SEQUTILS_INC) $(MISC_INC) $(TEXTUTILS_INC)

#EXTRA_LIBS += -ldmalloc

MWSRC=$(O)

OBJS	= \
	$(MWSRC)/make_weights.o

HLIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB) \
	$(MATH_LIB)


$(O)/make_weights: $(OBJS) $(DEPS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJS) $(HLIBS) $(LIBSC) 
#	purify $(CLD) -o $@ $(OBJS) $(HLIBS) $(LIBSC) 



DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

#include dependencies
