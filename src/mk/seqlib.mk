LIBS	= seqlib
PROGS	= $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include	$(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

SEQLIBBIN=$(O)

SRS_HOME=$(SRCROOT)/srs
#SRS_HOME=/home6/kfs/srs
SRS_BIN=$(SRS_HOME)/bin/$(O)
SRS_INC=-I$(SRS_HOME)/src -I$(SRS_BIN)
#SRS_LIB=$(LINK_PATHFLAG)$(SRS_BIN) $(LINK_LIBFLAG)srs$(LIB_EXT)
SRS_LIB=$(SRS_BIN)/libsrs.a

INCLUDES_E += $(SEQLIB_INC) $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) \
	      $(MISC_INC) $(SRS_INC) $(TK_INC)

CFLAGS += $(SHLIB_CFLAGS)

# SRS support
# 04/03/98 johnt - made SRS support optional
ifdef NO_SRS
SEQCLIBSRSOBJS =
else
SEQCLIBSRSOBJS =\
	$(SEQLIBBIN)/seq_browse_srs.o \
	$(SEQLIBBIN)/seqlib_file_srs.o \
	$(SRS_BIN)/local.o\
	$(SRS_LIB)
endif


# C sequence library IO
SEQCLIBOBJS = \
	$(SEQLIBBIN)/seq_browse.o \
	$(SEQLIBBIN)/seqlib_file.o \
	$(SEQLIBBIN)/seq_hits.o \
	$(SEQLIBBIN)/seq_result.o \
	$(SEQLIBBIN)/seqlib_globals.o \
	$(GCGBITS) \
	$(SEQCLIBSRSOBJS)

# Library + tcl interface
SEQLIBOBJS = \
	$(SEQCLIBOBJS) \
	$(SEQLIBBIN)/busy_dialog.o\
	$(SEQLIBBIN)/seq_cmds.o \
	$(SEQLIBSRSOBJS)

SEQLIB2OBJS = \
	$(SRS_BIN)/local.o\
	$(SEQLIBBIN)/seq_cmds2.o

GCGBITS=\
	$(SEQLIBBIN)/bit.o

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

# 7/1/99 johnt - added SEQUTILS_DEP ( for Visual C++)
DEPS = \
	$(TCL_DEP) \
	$(TKUTILS_DEP) \
	$(IOUTILS_DEP) \
	$(MISC_DEP) \
	$(SEQUTILS_DEP)

SEQLIB_LIBS = \
	$(SEQUTILS_LIB)

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(L)/$(SHLIB_PREFIX)seqlib$(SHLIB_SUFFIX): $(SEQLIBOBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(SEQLIBOBJS) $(DEPS) $(LINK_PATHFLAG)$(L) $(SEQLIB_LIBS)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for WINNT
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(SEQLIBOBJS)
	$(MKDEFL) $@ $(SEQLIBOBJS)


DEPEND_OBJ=$(SEQLIBOBJS)

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)
	-mkdir $(INSTALLLIB)/seqlib
	cp *.tcl tclIndex $(INSTALLLIB)/seqlib

distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl tclIndex Makefile slim dependencies \
		$(DIRNAME)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.

