#
# Makefile for text_utils routines
#

LIBS = text_utils
PROGS= $(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(TCLUTILS_INC)
CFLAGS += $(SHLIB_CFLAGS)
DEFINES += -DCHECK_LICENCE

TEXTOUTPUTBIN=$(O)

#
# Objects
#
OBJS = \
	$(TEXTOUTPUTBIN)/text_output_stubs.o \
	$(TEXTOUTPUTBIN)/text_output_stubs2.o\
	$(TEXTOUTPUTBIN)/check_licence.o\
	$(TEXTOUTPUTBIN)/licence_utils.o\
	$(TEXTOUTPUTBIN)/valid_seq.o\
	$(TEXTOUTPUTBIN)/md52.o

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for WINNT
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
# DO NOT DELETE
