#LIBS = sip
#PROGS = $(LIBS)

SRCROOT=..
include	 $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk
#PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SIP4BIN=$(O)

INCLUDES_E += $(SEQLIB_INC) $(TK_INC) $(IOLIB_INC) $(TKUTILS_INC) \
	      $(SEQUTILS_INC) $(SPIN_INC) $(MISC_INC)

SIPOBJS = \
	$(SIP4BIN)/compare_spans.o \
	$(SIP4BIN)/sip_hash.o \
	$(SIP4BIN)/sip_similar_spans.o\
	$(SIP4BIN)/sip_find_identity.o\
	$(SIP4BIN)/sip_quick_scan.o\
	$(SIP4BIN)/sip_align.o \
	$(SIP4BIN)/sip_globals.o \
	$(SIP4BIN)/sip_cmds.o \
	$(SIP4BIN)/sip_results.o \
	$(SIP4BIN)/rescan_matches.o \
	$(SIP4BIN)/probs.o \
	$(SIP4BIN)/sim.o \
	$(SIP4BIN)/sip_sim.o\
	$(SIP4BIN)/readpam.o

# 7/1/99 johnt - added additional dependencies ( for Visual C++)
DEPS=\
	$(IOLIB_DEP) \
	$(MISC_DEP) \
	$(SEQUTILS_DEP) \
	$(SPIN_DEP) \
	$(TKUTILS_DEP) \
	$(IOLIB_DEP) \
	$(SEQLIB_DEP) \
        $(TCL_DEP) \
	$(TK_DEP) \
	$(IOUTILS_DEP) \
	$(MISC_DEP) \
	$(SPIN_DEP)


SIP4_LIBS=\
	$(SEQLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(SPIN_LIB) \
	$(TKUTILS_LIB) \
	$(IOLIB_LIB) \
	$(TK_LIB) \
	$(IOUTILS_LIB) \
	$(MISC_LIB)


$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@


$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(SIPOBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(SIPOBJS) $(DEPS)


# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for WINNT
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(SIPOBJS)
	$(MKDEFL) $@ $(SIPOBJS)


DEPEND_OBJ = $(SIP4OBJS)

distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl tclIndex Makefile dependencies \
		$(DIRNAME)
	-rm $(DIRNAME)/resource.h

install:
	-mkdir $(INSTALLLIB)/sip
	cp *.tcl tclIndex $(INSTALLLIB)/sip
#	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.

