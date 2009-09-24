LIBS = haplo
PROGS = x$(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E := $(MISC_INC) $(GAP4_INC) $(G_INC) $(TK_INC) $(ITK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(INCLUDES_E)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	haplo_tcl.o \
	haplo_snps.o \
	haplo_split.o \
	haplo_cons.o

HAPLO_DEP=\
	$(TK_LIB) \
	$(TKUTILS_LIB) \
	$(MISC_LIB) \
	$(GAP_LIB)

x$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(HAPLO_DEP) $(SHLIB_DEP)

install:
	-mkdir $(INSTALLTCL)/haplo
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/haplo
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/haplorc $(INSTALLETC)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/Makefile $(S)/tclIndex $(S)/haplorc \
	   $(DIRNAME)

DEPEND_OBJ = $(OBJS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

haplo_cons.o: $(SRCROOT)/Misc/array.h
haplo_cons.o: $(SRCROOT)/Misc/bitmap.h
haplo_cons.o: $(SRCROOT)/Misc/dstring.h
haplo_cons.o: $(SRCROOT)/Misc/os.h
haplo_cons.o: $(SRCROOT)/Misc/xalloc.h
haplo_cons.o: $(SRCROOT)/Misc/xerror.h
haplo_cons.o: $(SRCROOT)/g/freetree.h
haplo_cons.o: $(SRCROOT)/g/g-defs.h
haplo_cons.o: $(SRCROOT)/g/g-error.h
haplo_cons.o: $(SRCROOT)/g/g-filedefs.h
haplo_cons.o: $(SRCROOT)/g/g-os.h
haplo_cons.o: $(SRCROOT)/g/g-struct.h
haplo_cons.o: $(SRCROOT)/gap4/IO.h
haplo_cons.o: $(SRCROOT)/gap4/gap-dbstruct.h
haplo_cons.o: $(SRCROOT)/gap4/gap-if.h
haplo_cons.o: $(SRCROOT)/gap4/gap_globals.h
haplo_cons.o: $(SRCROOT)/gap4/io_handle.h
haplo_cons.o: $(SRCROOT)/gap4/io_utils.h
haplo_cons.o: $(SRCROOT)/gap4/qual.h
haplo_cons.o: $(SRCROOT)/gap4/vseqs.h
haplo_cons.o: $(SRCROOT)/tk_utils/tcl_utils.h
haplo_cons.o: $(PWD)/staden_config.h
haplo_cons.o: haplo.h
haplo_snps.o: $(SRCROOT)/Misc/array.h
haplo_snps.o: $(SRCROOT)/Misc/bitmap.h
haplo_snps.o: $(SRCROOT)/Misc/dstring.h
haplo_snps.o: $(SRCROOT)/Misc/misc.h
haplo_snps.o: $(SRCROOT)/Misc/os.h
haplo_snps.o: $(SRCROOT)/Misc/xalloc.h
haplo_snps.o: $(SRCROOT)/Misc/xerror.h
haplo_snps.o: $(SRCROOT)/g/freetree.h
haplo_snps.o: $(SRCROOT)/g/g-defs.h
haplo_snps.o: $(SRCROOT)/g/g-error.h
haplo_snps.o: $(SRCROOT)/g/g-filedefs.h
haplo_snps.o: $(SRCROOT)/g/g-os.h
haplo_snps.o: $(SRCROOT)/g/g-struct.h
haplo_snps.o: $(SRCROOT)/gap4/IO.h
haplo_snps.o: $(SRCROOT)/gap4/gap-dbstruct.h
haplo_snps.o: $(SRCROOT)/gap4/gap-if.h
haplo_snps.o: $(SRCROOT)/gap4/gap_globals.h
haplo_snps.o: $(SRCROOT)/gap4/io_handle.h
haplo_snps.o: $(SRCROOT)/gap4/io_utils.h
haplo_snps.o: $(SRCROOT)/gap4/list.h
haplo_snps.o: $(SRCROOT)/gap4/notes.h
haplo_snps.o: $(SRCROOT)/gap4/qual.h
haplo_snps.o: $(SRCROOT)/gap4/template.h
haplo_snps.o: $(SRCROOT)/tk_utils/tcl_utils.h
haplo_snps.o: $(PWD)/staden_config.h
haplo_snps.o: haplo.h
haplo_split.o: $(SRCROOT)/Misc/array.h
haplo_split.o: $(SRCROOT)/Misc/bitmap.h
haplo_split.o: $(SRCROOT)/Misc/dstring.h
haplo_split.o: $(SRCROOT)/Misc/os.h
haplo_split.o: $(SRCROOT)/Misc/xalloc.h
haplo_split.o: $(SRCROOT)/Misc/xerror.h
haplo_split.o: $(SRCROOT)/g/freetree.h
haplo_split.o: $(SRCROOT)/g/g-defs.h
haplo_split.o: $(SRCROOT)/g/g-error.h
haplo_split.o: $(SRCROOT)/g/g-filedefs.h
haplo_split.o: $(SRCROOT)/g/g-os.h
haplo_split.o: $(SRCROOT)/g/g-struct.h
haplo_split.o: $(SRCROOT)/gap4/IO.h
haplo_split.o: $(SRCROOT)/gap4/gap-dbstruct.h
haplo_split.o: $(SRCROOT)/gap4/gap-if.h
haplo_split.o: $(SRCROOT)/gap4/io_handle.h
haplo_split.o: $(SRCROOT)/gap4/io_utils.h
haplo_split.o: $(PWD)/staden_config.h
haplo_split.o: haplo.h
haplo_tcl.o: $(SRCROOT)/Misc/array.h
haplo_tcl.o: $(SRCROOT)/Misc/bitmap.h
haplo_tcl.o: $(SRCROOT)/Misc/dstring.h
haplo_tcl.o: $(SRCROOT)/Misc/os.h
haplo_tcl.o: $(SRCROOT)/Misc/xalloc.h
haplo_tcl.o: $(SRCROOT)/Misc/xerror.h
haplo_tcl.o: $(SRCROOT)/g/freetree.h
haplo_tcl.o: $(SRCROOT)/g/g-defs.h
haplo_tcl.o: $(SRCROOT)/g/g-error.h
haplo_tcl.o: $(SRCROOT)/g/g-filedefs.h
haplo_tcl.o: $(SRCROOT)/g/g-os.h
haplo_tcl.o: $(SRCROOT)/g/g-struct.h
haplo_tcl.o: $(SRCROOT)/gap4/IO.h
haplo_tcl.o: $(SRCROOT)/gap4/gap-dbstruct.h
haplo_tcl.o: $(SRCROOT)/gap4/gap-if.h
haplo_tcl.o: $(SRCROOT)/gap4/gap_cli_arg.h
haplo_tcl.o: $(SRCROOT)/gap4/io_handle.h
haplo_tcl.o: $(SRCROOT)/gap4/io_utils.h
haplo_tcl.o: $(SRCROOT)/gap4/list_proc.h
haplo_tcl.o: $(SRCROOT)/tk_utils/cli_arg.h
haplo_tcl.o: $(PWD)/staden_config.h
haplo_tcl.o: haplo.h
