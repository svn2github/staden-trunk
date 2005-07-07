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
#DEFINES += -DCHECK_LICENCE

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
$(LIBS) : .links $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(MISC_LIB)

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

.PRECIOUS:	.links
.links:
	-ln -s ../licence/boxes.h .
	-ln -s ../licence/licence.h .
	-ln -s ../licence/check_licence.c .
	-ln -s ../licence/licence_utils.c .
	-ln -s ../licence/licence_utils.h .
	-ln -s ../licence/md52.c .
	-ln -s ../licence/md52.h .
	-ln -s ../licence/valid_seq.c .
	touch .links
	$(MAKE) $(MFLAGS)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
# DO NOT DELETE
