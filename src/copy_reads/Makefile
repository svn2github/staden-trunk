LIBS = copy_reads
PROGS = cplib
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

#SRCROOT=/nfs/arran/home5/pubseq/share/src
SRCROOT=..

include	 $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

# Local comment: Comment out next line for remote compilation
COPYREADSBIN=$(O)

#L=$(O)

INCLUDES_E := $(GAP4_INC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(P3_INC) $(INCLUDES_E)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	$(COPYREADSBIN)/copy_reads_globals.o\
	$(COPYREADSBIN)/copy_reads.o

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

install:
	-mkdir $(INSTALLLIB)/copy_reads
	-cp copy_reads $(INSTALLSCRIPT)
	cp copy_reads.tcl copy_reads_main.tcl tclIndex $(INSTALLLIB)/copy_reads
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

distsrc: distsrc_dirs
	-cp -R *.[ch] copy_reads.tcl copy_reads_main.tcl tclIndex \
		copy_reads Makefile dependencies \
		$(DIRNAME)

DEPEND_OBJ = $(OBJS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
