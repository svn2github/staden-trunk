LIBS = prefinish
PROGS = x$(LIBS) insert_size_dist read_size_dist read_conf_dist
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E := $(MISC_INC) $(GAP4_INC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(P3_INC) \
	      $(INCLUDES_E)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	finish.o \
	finish_hash.o \
	finish_long.o \
	finish_main.o \
	finish_utils.o \
	finish_walk.o \
	finish_reverse.o \
	finish_filter.o \
	finish_pcr.o \
	dust.o

FIN_DEP=\
	$(TK_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(MISC_LIB) \
	$(GAP_LIB) \
	$(P3_LIB)

x$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE) $(L)/.dir
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(FIN_DEP) $(SHLIB_DEP)

# Most of these are only needed because GAP_LIB has a dependency
ISLIBS=\
	$(G_LIB) \
	$(GAP_LIB) \
	$(TKUTILS_LIB) \
	$(P3_LIB) \
	$(MUT_LIB) \
	$(SEQUTILS_LIB) \
	$(MISC_LIB) \
	$(TCL_LIB) \
	$(IOLIB_LIB)

ISOBJS= insert_size_dist.o

RSOBJS=	read_size_dist.o

RCOBJS=	read_conf_dist.o

insert_size_dist: $(ISOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(ISOBJS) $(ISLIBS) $(LIBSC)

read_size_dist: $(RSOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(RSOBJS) $(ISLIBS) $(LIBSC)

read_conf_dist: $(RCOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(RCOBJS) $(ISLIBS) $(LIBSC)

S=$(SRCROOT)/prefinish
install:
	-mkdir $(INSTALLTCL)/prefinish
	cp $(S)/*.tcl $(S)/*.itcl $(S)/*.itk $(S)/help $(S)/args.template $(S)/tclIndex $(INSTALLTCL)/prefinish
	cp $(S)/config.prefinish $(INSTALLTCL)/prefinish
	cp $(PROGLIBS) $(INSTALLLIB)
	$(INSTALL) prefinish $(INSTALLBIN)
	cp $(S)/finish_sanger $(S)/finish_cDNA $(S)/finish_cDNA_ends_only \
	    $(INSTALLSCRIPTS)
	-mkdir $(INSTALLTCL)/prefinish/rules
	cp $(S)/rules/*.rule $(INSTALLTCL)/prefinish/rules
	cp $(S)/config.prefinish $(S)/prefinishrc $(INSTALLETC)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/*.itcl $(S)/*.itk $(S)/help \
	   $(S)/args.template $(S)/Makefile $(S)/prefinish $(S)/tclIndex \
	   $(S)/config.prefinish $(S)/finish_sanger $(S)/finish_cDNA \
	   $(S)/finish_cDNA_ends_only $(S)/prefinishrc \
	   $(DIRNAME)
	-mkdir $(DIRNAME)/rules
	cp $(S)/rules/* $(DIRNAME)/rules
	-mkdir $(DIRNAME)/docs
	cp $(S)/docs/* $(DIRNAME)/docs

DEPEND_OBJ = $(OBJS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

dust.o: $(SRCROOT)/seq_utils/dna_utils.h
finish.o: $(SRCROOT)/Misc/array.h
finish.o: $(SRCROOT)/Misc/bitmap.h
finish.o: $(SRCROOT)/Misc/dstring.h
finish.o: $(SRCROOT)/Misc/misc.h
finish.o: $(SRCROOT)/Misc/os.h
finish.o: $(SRCROOT)/Misc/reg_exp.h
finish.o: $(SRCROOT)/Misc/xalloc.h
finish.o: $(SRCROOT)/Misc/xerror.h
finish.o: $(PWD)/staden_config.h
finish.o: $(SRCROOT)/g/freetree.h
finish.o: $(SRCROOT)/g/g-defs.h
finish.o: $(SRCROOT)/g/g-error.h
finish.o: $(SRCROOT)/g/g-filedefs.h
finish.o: $(SRCROOT)/g/g-os.h
finish.o: $(SRCROOT)/g/g-struct.h
finish.o: $(SRCROOT)/gap4/IO.h
finish.o: $(SRCROOT)/gap4/consen.h
finish.o: $(SRCROOT)/gap4/edStructs.h
finish.o: $(SRCROOT)/gap4/edUtils.h
finish.o: $(SRCROOT)/gap4/fort.h
finish.o: $(SRCROOT)/gap4/fortran.h
finish.o: $(SRCROOT)/gap4/gap-dbstruct.h
finish.o: $(SRCROOT)/gap4/gap-if.h
finish.o: $(SRCROOT)/gap4/gap_globals.h
finish.o: $(SRCROOT)/gap4/hash_lib.h
finish.o: $(SRCROOT)/gap4/io-reg.h
finish.o: $(SRCROOT)/gap4/io_handle.h
finish.o: $(SRCROOT)/gap4/io_utils.h
finish.o: $(SRCROOT)/gap4/list.h
finish.o: $(SRCROOT)/gap4/list_proc.h
finish.o: $(SRCROOT)/gap4/primlib.h
finish.o: $(SRCROOT)/gap4/qual.h
finish.o: $(SRCROOT)/gap4/tagDefs.h
finish.o: $(SRCROOT)/gap4/tagUtils.h
finish.o: $(SRCROOT)/gap4/template.h
finish.o: $(SRCROOT)/gap4/tkEdNames.h
finish.o: $(SRCROOT)/gap4/tkEditor.h
finish.o: $(SRCROOT)/gap4/tman_display.h
finish.o: $(SRCROOT)/gap4/undo.h
finish.o: $(SRCROOT)/gap4/vseqs.h
finish.o: $(SRCROOT)/prefinish/dust.h
finish.o: $(SRCROOT)/prefinish/finish.h
finish.o: $(SRCROOT)/prefinish/finish_filter.h
finish.o: $(SRCROOT)/prefinish/finish_hash.h
finish.o: $(SRCROOT)/prefinish/finish_long.h
finish.o: $(SRCROOT)/prefinish/finish_main.h
finish.o: $(SRCROOT)/prefinish/finish_pcr.h
finish.o: $(SRCROOT)/prefinish/finish_utils.h
finish.o: $(SRCROOT)/prefinish/finish_walk.h
finish.o: $(SRCROOT)/prefinish/gap_cli_arg.h
finish.o: $(SRCROOT)/primer3/src/dpal.h
finish.o: $(SRCROOT)/primer3/src/primer3.h
finish.o: $(SRCROOT)/seq_utils/align_lib.h
finish.o: $(SRCROOT)/seq_utils/dna_utils.h
finish.o: $(SRCROOT)/seq_utils/sequence_formats.h
finish.o: $(SRCROOT)/tk_utils/cli_arg.h
finish.o: $(SRCROOT)/tk_utils/intrinsic_type.h
finish.o: $(SRCROOT)/tk_utils/postscript.h
finish.o: $(SRCROOT)/tk_utils/sheet.h
finish.o: $(SRCROOT)/tk_utils/tcl_utils.h
finish.o: $(SRCROOT)/tk_utils/tkSheet.h
finish.o: $(SRCROOT)/tk_utils/tkSheet_common.h
finish.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
finish.o: $(SRCROOT)/tk_utils/tkTrace.h
finish_filter.o: $(SRCROOT)/Misc/array.h
finish_filter.o: $(SRCROOT)/Misc/bitmap.h
finish_filter.o: $(SRCROOT)/Misc/os.h
finish_filter.o: $(SRCROOT)/Misc/xalloc.h
finish_filter.o: $(SRCROOT)/Misc/xerror.h
finish_filter.o: $(PWD)/staden_config.h
finish_filter.o: $(SRCROOT)/g/freetree.h
finish_filter.o: $(SRCROOT)/g/g-defs.h
finish_filter.o: $(SRCROOT)/g/g-error.h
finish_filter.o: $(SRCROOT)/g/g-filedefs.h
finish_filter.o: $(SRCROOT)/g/g-os.h
finish_filter.o: $(SRCROOT)/g/g-struct.h
finish_filter.o: $(SRCROOT)/gap4/IO.h
finish_filter.o: $(SRCROOT)/gap4/gap-dbstruct.h
finish_filter.o: $(SRCROOT)/gap4/gap-if.h
finish_filter.o: $(SRCROOT)/gap4/hash_lib.h
finish_filter.o: $(SRCROOT)/gap4/io_handle.h
finish_filter.o: $(SRCROOT)/gap4/io_utils.h
finish_filter.o: $(SRCROOT)/gap4/list.h
finish_filter.o: $(SRCROOT)/gap4/primlib.h
finish_filter.o: $(SRCROOT)/gap4/qual.h
finish_filter.o: $(SRCROOT)/gap4/template.h
finish_filter.o: $(SRCROOT)/gap4/vseqs.h
finish_filter.o: $(SRCROOT)/prefinish/dust.h
finish_filter.o: $(SRCROOT)/prefinish/finish.h
finish_filter.o: $(SRCROOT)/prefinish/finish_filter.h
finish_filter.o: $(SRCROOT)/prefinish/finish_hash.h
finish_filter.o: $(SRCROOT)/primer3/src/dpal.h
finish_filter.o: $(SRCROOT)/primer3/src/primer3.h
finish_filter.o: $(SRCROOT)/seq_utils/align_lib.h
finish_filter.o: $(SRCROOT)/seq_utils/dna_utils.h
finish_filter.o: $(SRCROOT)/seq_utils/filter_words.h
finish_hash.o: $(SRCROOT)/Misc/misc.h
finish_hash.o: $(SRCROOT)/Misc/os.h
finish_hash.o: $(SRCROOT)/Misc/xalloc.h
finish_hash.o: $(PWD)/staden_config.h
finish_hash.o: $(SRCROOT)/gap4/hash_lib.h
finish_hash.o: $(SRCROOT)/prefinish/finish_hash.h
finish_hash.o: $(SRCROOT)/seq_utils/align_lib.h
finish_hash.o: $(SRCROOT)/seq_utils/dna_utils.h
finish_long.o: $(SRCROOT)/Misc/array.h
finish_long.o: $(SRCROOT)/Misc/bitmap.h
finish_long.o: $(SRCROOT)/Misc/misc.h
finish_long.o: $(SRCROOT)/Misc/os.h
finish_long.o: $(SRCROOT)/Misc/xalloc.h
finish_long.o: $(SRCROOT)/Misc/xerror.h
finish_long.o: $(PWD)/staden_config.h
finish_long.o: $(SRCROOT)/g/freetree.h
finish_long.o: $(SRCROOT)/g/g-defs.h
finish_long.o: $(SRCROOT)/g/g-error.h
finish_long.o: $(SRCROOT)/g/g-filedefs.h
finish_long.o: $(SRCROOT)/g/g-os.h
finish_long.o: $(SRCROOT)/g/g-struct.h
finish_long.o: $(SRCROOT)/gap4/IO.h
finish_long.o: $(SRCROOT)/gap4/gap-dbstruct.h
finish_long.o: $(SRCROOT)/gap4/gap-if.h
finish_long.o: $(SRCROOT)/gap4/hash_lib.h
finish_long.o: $(SRCROOT)/gap4/io_handle.h
finish_long.o: $(SRCROOT)/gap4/io_utils.h
finish_long.o: $(SRCROOT)/gap4/list.h
finish_long.o: $(SRCROOT)/gap4/primlib.h
finish_long.o: $(SRCROOT)/gap4/qual.h
finish_long.o: $(SRCROOT)/gap4/template.h
finish_long.o: $(SRCROOT)/gap4/vseqs.h
finish_long.o: $(SRCROOT)/prefinish/finish.h
finish_long.o: $(SRCROOT)/prefinish/finish_hash.h
finish_long.o: $(SRCROOT)/prefinish/finish_long.h
finish_long.o: $(SRCROOT)/prefinish/finish_utils.h
finish_long.o: $(SRCROOT)/primer3/src/dpal.h
finish_long.o: $(SRCROOT)/primer3/src/primer3.h
finish_long.o: $(SRCROOT)/seq_utils/align_lib.h
finish_main.o: $(SRCROOT)/Misc/array.h
finish_main.o: $(SRCROOT)/Misc/bitmap.h
finish_main.o: $(SRCROOT)/Misc/dstring.h
finish_main.o: $(SRCROOT)/Misc/misc.h
finish_main.o: $(SRCROOT)/Misc/os.h
finish_main.o: $(SRCROOT)/Misc/xalloc.h
finish_main.o: $(SRCROOT)/Misc/xerror.h
finish_main.o: $(PWD)/staden_config.h
finish_main.o: $(SRCROOT)/g/freetree.h
finish_main.o: $(SRCROOT)/g/g-defs.h
finish_main.o: $(SRCROOT)/g/g-error.h
finish_main.o: $(SRCROOT)/g/g-filedefs.h
finish_main.o: $(SRCROOT)/g/g-os.h
finish_main.o: $(SRCROOT)/g/g-struct.h
finish_main.o: $(SRCROOT)/gap4/IO.h
finish_main.o: $(SRCROOT)/gap4/edStructs.h
finish_main.o: $(SRCROOT)/gap4/edUtils.h
finish_main.o: $(SRCROOT)/gap4/find_fragments.h
finish_main.o: $(SRCROOT)/gap4/fort.h
finish_main.o: $(SRCROOT)/gap4/fortran.h
finish_main.o: $(SRCROOT)/gap4/gap-dbstruct.h
finish_main.o: $(SRCROOT)/gap4/gap-if.h
finish_main.o: $(SRCROOT)/gap4/gap_globals.h
finish_main.o: $(SRCROOT)/gap4/hash_lib.h
finish_main.o: $(SRCROOT)/gap4/io-reg.h
finish_main.o: $(SRCROOT)/gap4/io_handle.h
finish_main.o: $(SRCROOT)/gap4/io_utils.h
finish_main.o: $(SRCROOT)/gap4/list.h
finish_main.o: $(SRCROOT)/gap4/primlib.h
finish_main.o: $(SRCROOT)/gap4/qual.h
finish_main.o: $(SRCROOT)/gap4/tagDefs.h
finish_main.o: $(SRCROOT)/gap4/tagUtils.h
finish_main.o: $(SRCROOT)/gap4/template.h
finish_main.o: $(SRCROOT)/gap4/tkEdNames.h
finish_main.o: $(SRCROOT)/gap4/tkEditor.h
finish_main.o: $(SRCROOT)/gap4/tman_display.h
finish_main.o: $(SRCROOT)/gap4/undo.h
finish_main.o: $(SRCROOT)/gap4/vseqs.h
finish_main.o: $(SRCROOT)/prefinish/finish.h
finish_main.o: $(SRCROOT)/prefinish/finish_filter.h
finish_main.o: $(SRCROOT)/prefinish/finish_hash.h
finish_main.o: $(SRCROOT)/prefinish/finish_long.h
finish_main.o: $(SRCROOT)/prefinish/finish_reverse.h
finish_main.o: $(SRCROOT)/prefinish/finish_utils.h
finish_main.o: $(SRCROOT)/prefinish/finish_walk.h
finish_main.o: $(SRCROOT)/primer3/src/dpal.h
finish_main.o: $(SRCROOT)/primer3/src/primer3.h
finish_main.o: $(SRCROOT)/seq_utils/align_lib.h
finish_main.o: $(SRCROOT)/tk_utils/cli_arg.h
finish_main.o: $(SRCROOT)/tk_utils/intrinsic_type.h
finish_main.o: $(SRCROOT)/tk_utils/postscript.h
finish_main.o: $(SRCROOT)/tk_utils/sheet.h
finish_main.o: $(SRCROOT)/tk_utils/tcl_utils.h
finish_main.o: $(SRCROOT)/tk_utils/tkSheet.h
finish_main.o: $(SRCROOT)/tk_utils/tkSheet_common.h
finish_main.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
finish_main.o: $(SRCROOT)/tk_utils/tkTrace.h
finish_pcr.o: $(SRCROOT)/Misc/array.h
finish_pcr.o: $(SRCROOT)/Misc/bitmap.h
finish_pcr.o: $(SRCROOT)/Misc/dstring.h
finish_pcr.o: $(SRCROOT)/Misc/misc.h
finish_pcr.o: $(SRCROOT)/Misc/os.h
finish_pcr.o: $(SRCROOT)/Misc/xalloc.h
finish_pcr.o: $(SRCROOT)/Misc/xerror.h
finish_pcr.o: $(PWD)/staden_config.h
finish_pcr.o: $(SRCROOT)/g/freetree.h
finish_pcr.o: $(SRCROOT)/g/g-defs.h
finish_pcr.o: $(SRCROOT)/g/g-error.h
finish_pcr.o: $(SRCROOT)/g/g-filedefs.h
finish_pcr.o: $(SRCROOT)/g/g-os.h
finish_pcr.o: $(SRCROOT)/g/g-struct.h
finish_pcr.o: $(SRCROOT)/gap4/IO.h
finish_pcr.o: $(SRCROOT)/gap4/gap-dbstruct.h
finish_pcr.o: $(SRCROOT)/gap4/gap-if.h
finish_pcr.o: $(SRCROOT)/gap4/gap_globals.h
finish_pcr.o: $(SRCROOT)/gap4/hash_lib.h
finish_pcr.o: $(SRCROOT)/gap4/io_handle.h
finish_pcr.o: $(SRCROOT)/gap4/io_utils.h
finish_pcr.o: $(SRCROOT)/gap4/list.h
finish_pcr.o: $(SRCROOT)/gap4/primlib.h
finish_pcr.o: $(SRCROOT)/gap4/qual.h
finish_pcr.o: $(SRCROOT)/gap4/template.h
finish_pcr.o: $(SRCROOT)/gap4/vseqs.h
finish_pcr.o: $(SRCROOT)/prefinish/finish.h
finish_pcr.o: $(SRCROOT)/prefinish/finish_filter.h
finish_pcr.o: $(SRCROOT)/prefinish/finish_hash.h
finish_pcr.o: $(SRCROOT)/prefinish/finish_pcr.h
finish_pcr.o: $(SRCROOT)/prefinish/finish_utils.h
finish_pcr.o: $(SRCROOT)/prefinish/gap_cli_arg.h
finish_pcr.o: $(SRCROOT)/primer3/src/dpal.h
finish_pcr.o: $(SRCROOT)/primer3/src/primer3.h
finish_pcr.o: $(SRCROOT)/seq_utils/align_lib.h
finish_pcr.o: $(SRCROOT)/seq_utils/dna_utils.h
finish_pcr.o: $(SRCROOT)/tk_utils/cli_arg.h
finish_pcr.o: $(SRCROOT)/tk_utils/tcl_utils.h
finish_pcr.o: $(SRCROOT)/tk_utils/text_output.h
finish_reverse.o: $(SRCROOT)/Misc/array.h
finish_reverse.o: $(SRCROOT)/Misc/bitmap.h
finish_reverse.o: $(SRCROOT)/Misc/misc.h
finish_reverse.o: $(SRCROOT)/Misc/os.h
finish_reverse.o: $(SRCROOT)/Misc/xalloc.h
finish_reverse.o: $(SRCROOT)/Misc/xerror.h
finish_reverse.o: $(PWD)/staden_config.h
finish_reverse.o: $(SRCROOT)/g/freetree.h
finish_reverse.o: $(SRCROOT)/g/g-defs.h
finish_reverse.o: $(SRCROOT)/g/g-error.h
finish_reverse.o: $(SRCROOT)/g/g-filedefs.h
finish_reverse.o: $(SRCROOT)/g/g-os.h
finish_reverse.o: $(SRCROOT)/g/g-struct.h
finish_reverse.o: $(SRCROOT)/gap4/IO.h
finish_reverse.o: $(SRCROOT)/gap4/gap-dbstruct.h
finish_reverse.o: $(SRCROOT)/gap4/gap-if.h
finish_reverse.o: $(SRCROOT)/gap4/hash_lib.h
finish_reverse.o: $(SRCROOT)/gap4/io_handle.h
finish_reverse.o: $(SRCROOT)/gap4/io_utils.h
finish_reverse.o: $(SRCROOT)/gap4/list.h
finish_reverse.o: $(SRCROOT)/gap4/primlib.h
finish_reverse.o: $(SRCROOT)/gap4/qual.h
finish_reverse.o: $(SRCROOT)/gap4/template.h
finish_reverse.o: $(SRCROOT)/gap4/vseqs.h
finish_reverse.o: $(SRCROOT)/prefinish/finish.h
finish_reverse.o: $(SRCROOT)/prefinish/finish_hash.h
finish_reverse.o: $(SRCROOT)/prefinish/finish_reverse.h
finish_reverse.o: $(SRCROOT)/prefinish/finish_utils.h
finish_reverse.o: $(SRCROOT)/primer3/src/dpal.h
finish_reverse.o: $(SRCROOT)/primer3/src/primer3.h
finish_reverse.o: $(SRCROOT)/seq_utils/align_lib.h
finish_utils.o: $(SRCROOT)/Misc/array.h
finish_utils.o: $(SRCROOT)/Misc/bitmap.h
finish_utils.o: $(SRCROOT)/Misc/dstring.h
finish_utils.o: $(SRCROOT)/Misc/misc.h
finish_utils.o: $(SRCROOT)/Misc/os.h
finish_utils.o: $(SRCROOT)/Misc/xalloc.h
finish_utils.o: $(SRCROOT)/Misc/xerror.h
finish_utils.o: $(PWD)/staden_config.h
finish_utils.o: $(SRCROOT)/g/freetree.h
finish_utils.o: $(SRCROOT)/g/g-defs.h
finish_utils.o: $(SRCROOT)/g/g-error.h
finish_utils.o: $(SRCROOT)/g/g-filedefs.h
finish_utils.o: $(SRCROOT)/g/g-os.h
finish_utils.o: $(SRCROOT)/g/g-struct.h
finish_utils.o: $(SRCROOT)/gap4/IO.h
finish_utils.o: $(SRCROOT)/gap4/edStructs.h
finish_utils.o: $(SRCROOT)/gap4/edUtils.h
finish_utils.o: $(SRCROOT)/gap4/fort.h
finish_utils.o: $(SRCROOT)/gap4/fortran.h
finish_utils.o: $(SRCROOT)/gap4/gap-dbstruct.h
finish_utils.o: $(SRCROOT)/gap4/gap-if.h
finish_utils.o: $(SRCROOT)/gap4/hash_lib.h
finish_utils.o: $(SRCROOT)/gap4/io-reg.h
finish_utils.o: $(SRCROOT)/gap4/io_handle.h
finish_utils.o: $(SRCROOT)/gap4/io_utils.h
finish_utils.o: $(SRCROOT)/gap4/list.h
finish_utils.o: $(SRCROOT)/gap4/primlib.h
finish_utils.o: $(SRCROOT)/gap4/qual.h
finish_utils.o: $(SRCROOT)/gap4/tagDefs.h
finish_utils.o: $(SRCROOT)/gap4/tagUtils.h
finish_utils.o: $(SRCROOT)/gap4/template.h
finish_utils.o: $(SRCROOT)/gap4/tkEdNames.h
finish_utils.o: $(SRCROOT)/gap4/tkEditor.h
finish_utils.o: $(SRCROOT)/gap4/tman_display.h
finish_utils.o: $(SRCROOT)/gap4/undo.h
finish_utils.o: $(SRCROOT)/gap4/vseqs.h
finish_utils.o: $(SRCROOT)/prefinish/finish.h
finish_utils.o: $(SRCROOT)/prefinish/finish_hash.h
finish_utils.o: $(SRCROOT)/prefinish/finish_utils.h
finish_utils.o: $(SRCROOT)/primer3/src/dpal.h
finish_utils.o: $(SRCROOT)/primer3/src/primer3.h
finish_utils.o: $(SRCROOT)/seq_utils/align_lib.h
finish_utils.o: $(SRCROOT)/tk_utils/cli_arg.h
finish_utils.o: $(SRCROOT)/tk_utils/intrinsic_type.h
finish_utils.o: $(SRCROOT)/tk_utils/postscript.h
finish_utils.o: $(SRCROOT)/tk_utils/sheet.h
finish_utils.o: $(SRCROOT)/tk_utils/tkSheet.h
finish_utils.o: $(SRCROOT)/tk_utils/tkSheet_common.h
finish_utils.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
finish_utils.o: $(SRCROOT)/tk_utils/tkTrace.h
finish_walk.o: $(SRCROOT)/Misc/array.h
finish_walk.o: $(SRCROOT)/Misc/bitmap.h
finish_walk.o: $(SRCROOT)/Misc/misc.h
finish_walk.o: $(SRCROOT)/Misc/os.h
finish_walk.o: $(SRCROOT)/Misc/xalloc.h
finish_walk.o: $(SRCROOT)/Misc/xerror.h
finish_walk.o: $(PWD)/staden_config.h
finish_walk.o: $(SRCROOT)/g/freetree.h
finish_walk.o: $(SRCROOT)/g/g-defs.h
finish_walk.o: $(SRCROOT)/g/g-error.h
finish_walk.o: $(SRCROOT)/g/g-filedefs.h
finish_walk.o: $(SRCROOT)/g/g-os.h
finish_walk.o: $(SRCROOT)/g/g-struct.h
finish_walk.o: $(SRCROOT)/gap4/IO.h
finish_walk.o: $(SRCROOT)/gap4/gap-dbstruct.h
finish_walk.o: $(SRCROOT)/gap4/gap-if.h
finish_walk.o: $(SRCROOT)/gap4/gap_globals.h
finish_walk.o: $(SRCROOT)/gap4/hash_lib.h
finish_walk.o: $(SRCROOT)/gap4/io_handle.h
finish_walk.o: $(SRCROOT)/gap4/io_utils.h
finish_walk.o: $(SRCROOT)/gap4/list.h
finish_walk.o: $(SRCROOT)/gap4/primlib.h
finish_walk.o: $(SRCROOT)/gap4/qual.h
finish_walk.o: $(SRCROOT)/gap4/template.h
finish_walk.o: $(SRCROOT)/gap4/vseqs.h
finish_walk.o: $(SRCROOT)/prefinish/finish.h
finish_walk.o: $(SRCROOT)/prefinish/finish_filter.h
finish_walk.o: $(SRCROOT)/prefinish/finish_hash.h
finish_walk.o: $(SRCROOT)/prefinish/finish_utils.h
finish_walk.o: $(SRCROOT)/prefinish/finish_walk.h
finish_walk.o: $(SRCROOT)/primer3/src/dpal.h
finish_walk.o: $(SRCROOT)/primer3/src/primer3.h
finish_walk.o: $(SRCROOT)/seq_utils/align_lib.h
finish_walk.o: $(SRCROOT)/seq_utils/dna_utils.h
finish_walk.o: $(SRCROOT)/seq_utils/search_utils.h
finish_walk.o: $(SRCROOT)/tk_utils/tcl_utils.h
