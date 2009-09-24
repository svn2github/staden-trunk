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


install:
	-mkdir $(INSTALLTCL)/copy_reads
	cp $(S)/copy_reads $(INSTALLBIN)
	cp $(S)/copy_reads.tcl $(S)/copy_reads_main.tcl $(S)/tclIndex $(INSTALLTCL)/copy_reads
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/copy_readsrc $(INSTALLETC)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/tclIndex $(S)/copy_reads \
	   $(S)/Makefile $(S)/copy_readsrc $(DIRNAME)

DEPEND_OBJ = $(OBJS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

copy_reads.o: $(SRCROOT)/Misc/array.h
copy_reads.o: $(SRCROOT)/Misc/bitmap.h
copy_reads.o: $(SRCROOT)/Misc/dstring.h
copy_reads.o: $(SRCROOT)/Misc/misc.h
copy_reads.o: $(SRCROOT)/Misc/os.h
copy_reads.o: $(SRCROOT)/Misc/xalloc.h
copy_reads.o: $(SRCROOT)/Misc/xerror.h
copy_reads.o: $(SRCROOT)/g/freetree.h
copy_reads.o: $(SRCROOT)/g/g-defs.h
copy_reads.o: $(SRCROOT)/g/g-error.h
copy_reads.o: $(SRCROOT)/g/g-filedefs.h
copy_reads.o: $(SRCROOT)/g/g-os.h
copy_reads.o: $(SRCROOT)/g/g-struct.h
copy_reads.o: $(SRCROOT)/gap4/IO.h
copy_reads.o: $(SRCROOT)/gap4/active_tags.h
copy_reads.o: $(SRCROOT)/gap4/assemble_direct.h
copy_reads.o: $(SRCROOT)/gap4/complement.h
copy_reads.o: $(SRCROOT)/gap4/consen.h
copy_reads.o: $(SRCROOT)/gap4/edStructs.h
copy_reads.o: $(SRCROOT)/gap4/edUtils.h
copy_reads.o: $(SRCROOT)/gap4/fij.h
copy_reads.o: $(SRCROOT)/gap4/fort.h
copy_reads.o: $(SRCROOT)/gap4/fortran.h
copy_reads.o: $(SRCROOT)/gap4/gap-dbstruct.h
copy_reads.o: $(SRCROOT)/gap4/gap-if.h
copy_reads.o: $(SRCROOT)/gap4/gap_cli_arg.h
copy_reads.o: $(SRCROOT)/gap4/gap_globals.h
copy_reads.o: $(SRCROOT)/gap4/hash_lib.h
copy_reads.o: $(SRCROOT)/gap4/io-reg.h
copy_reads.o: $(SRCROOT)/gap4/io_handle.h
copy_reads.o: $(SRCROOT)/gap4/io_utils.h
copy_reads.o: $(SRCROOT)/gap4/list.h
copy_reads.o: $(SRCROOT)/gap4/list_proc.h
copy_reads.o: $(SRCROOT)/gap4/primlib.h
copy_reads.o: $(SRCROOT)/gap4/qual.h
copy_reads.o: $(SRCROOT)/gap4/seqInfo.h
copy_reads.o: $(SRCROOT)/gap4/tagDefs.h
copy_reads.o: $(SRCROOT)/gap4/tagUtils.h
copy_reads.o: $(SRCROOT)/gap4/tagdb.h
copy_reads.o: $(SRCROOT)/gap4/template.h
copy_reads.o: $(SRCROOT)/gap4/tkEdNames.h
copy_reads.o: $(SRCROOT)/gap4/tkEditor.h
copy_reads.o: $(SRCROOT)/gap4/tman_display.h
copy_reads.o: $(SRCROOT)/gap4/undo.h
copy_reads.o: $(SRCROOT)/primer3/src/dpal.h
copy_reads.o: $(SRCROOT)/primer3/src/primer3.h
copy_reads.o: $(SRCROOT)/seq_utils/align.h
copy_reads.o: $(SRCROOT)/seq_utils/align_lib.h
copy_reads.o: $(SRCROOT)/seq_utils/align_lib_old.h
copy_reads.o: $(SRCROOT)/seq_utils/dna_utils.h
copy_reads.o: $(SRCROOT)/tk_utils/cli_arg.h
copy_reads.o: $(SRCROOT)/tk_utils/postscript.h
copy_reads.o: $(SRCROOT)/tk_utils/sheet.h
copy_reads.o: $(SRCROOT)/tk_utils/tcl_utils.h
copy_reads.o: $(SRCROOT)/tk_utils/tkSheet.h
copy_reads.o: $(SRCROOT)/tk_utils/tkSheet_common.h
copy_reads.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
copy_reads.o: $(SRCROOT)/tk_utils/tkTrace.h
copy_reads.o: $(PWD)/staden_config.h
copy_reads.o: copy_reads.h
copy_reads.o: copy_reads_globals.h
