LIBS = copy_reads
PROGS = cplib
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

# Local comment: Comment out next line for remote compilation

INCLUDES_E := $(MISC_INC) $(GAP4_INC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(P3_INC) $(INCLUDES_E)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	copy_reads_globals.o\
	copy_reads.o

COPIER_LIB=\
	$(IOLIB_LIB) \
	$(TK_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(MISC_LIB) \
	$(GAP_LIB)

cplib : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	-$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(COPIER_LIB) $(SHLIB_DEP) $(F77_DEP)

S=$(SRCROOT)/copy_reads
install:
	-mkdir $(INSTALLTCL)/copy_reads
	-cp $(S)/copy_reads $(INSTALLBIN)
	cp $(S)/copy_reads.tcl $(S)/copy_reads_main.tcl $(S)/tclIndex $(INSTALLTCL)/copy_reads
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/copy_readsrc $(INSTALLETC)

distsrc: distsrc_dirs
	-cp -R *.[ch] copy_reads.tcl copy_reads_main.tcl tclIndex \
		copy_reads Makefile dependencies \
		$(DIRNAME)

DEPEND_OBJ = $(OBJS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
