# Makefile for the composition 'package' to add to gap4.

LIBS = composition
PROGS= $(LIBS)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INSTALLDIR  = ./install

INCLUDES_E += $(TCL_INC) $(TKUTILS_INC) $(GAP4_INC) $(G_INC)
CFLAGS	   += $(SHLIB_CFLAGS)

TESTBIN	    = $(O)
L	    = $(INSTALLDIR)/$(O)

# Objects
OBJS = \
	$(TESTBIN)/composition.o

CLIBS = \
	$(G_LIB) \
	$(TKUTILS_LIB) \
	$(TCL_LIB)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS)
	-mkdir $(INSTALLDIR)
	-mkdir $(INSTALLDIR)/$(O)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $@ $(OBJS) $(CLIBS)

DEPEND_OBJ = $(OBJS)

install: $(LIBS)
	cp tclIndex composition.tcl compositionrc composition.topic \
	composition.index composition.html $(INSTALLDIR)

include dependencies
