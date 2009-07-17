LIBS = gap
ifeq ($(MACHINE),windows)
PROGS = $(LIBS) copy_db 
else
PROGS = $(LIBS) copy_db 
endif
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

# Local comment: Comment out next line for remote compilation

INCLUDES_E += $(MISC_INC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(MUT_INC) \
	      $(F77_INC) $(P3_INC)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)


#
# Building the programs
# This should be just a linking phase because all of the object
# files and library files are generated using implicit rules.
# We use the fortran compiler to do linking.
#

GAPDB_LOW=\
	gap-local.o\
	gap-remote.o\
	gap-if.o\
	gap-init.o\
	gap-dbstruct.o\
	gap-create.o\
	gap-error.o\
	stack_dump.o\
	gap-io.o

GAPDB_MID=\
        IO.o \
        io-reg.o \
        actf.o

GAPDB_UTILS=\
	io_handle.o \
	io_utils.o

# GAPDB_EXT_OBJS is basically the low and mid level files plus the
# text-io-reg.o object. This is simply stub routines to enable safe linking
# and isn't needed by gap4 itself.

GAPDB_EXT_OBJS=\
	$(GAPDB_LOW) \
	$(GAPDB_MID) \
	$(GAPDB_UTILS) \
	text-io-reg.o

GAPDB_EXT_INC=$(G_INC) -I$(GAP4SRC)

GAPDB_EXT_DEPS=\
	$(G_DEP) \
	$(TEXTUTILS_DEP) \
	$(MISC_DEP) \
	$(TCL_DEP)

GAPDB_EXT_LIBS=\
	$(G_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB) \
	$(TCL_LIB)

GAPDB=\
	$(GAPDB_LOW)\
	$(GAPDB_MID)\
	IO2.o\
	seqInfo.o\
	parse_ft.o\
	IO3.o \
	io_utils.o \
	io_handle.o \
	gap-tcl.o \
	tk-io-reg.o

CONSEN=\
	consen.o\
	gap_hash.o\
	qual.o\
	qualIO.o

TCLBITS=\
	gap_cli_arg.o \
	gap_globals.o \
	newgap_cmds.o \
	init.o

EDITOR = \
	tkEditor.o\
	tkEdNames.o\
	tkEdUtils.o\
	edInterface.o\
	edUtils2.o \
	tagU1.o\
	undo.o\
	edExtend.o\
	edCommands.o\
	edMutations.o\
	tagEditor.o\
	searchUtils.o\
	tman_interface.o\
	tman_display.o\
	tman_cons.o\
	tman_diff.o\
	contigEditor.o\
	join.o\
	oligo.o

GAPFUNCS = \
	bubbl3.o\
	tagdb.o\
	notedb.o\
	active_tags.o\
	dbcheck.o\
	clones.o\
	extract.o\
	preass.o \
	list.o \
	reactions.o \
	probe.o \
	template.o \
	template_display.o \
	ruler_display.o \
	gap_canvas_box.o \
	hash.o \
	gap_array.o \
	show_relationships.o \
	fij.o \
	hash_lib.o \
	do_fij.o \
	auto_assemble.o \
	dis_readings.o \
	find_repeats.o \
	break_contig.o \
	quality_plot.o \
	readpair.o \
	contig_selector.o \
	complement.o \
	cs-object.o\
	list_proc.o\
	dstrand.o\
	oligo_sel.o\
	primlib.o\
	alter_rel.o\
	restriction_enzymes.o \
	stop_codon.o \
	assemble_direct.o\
	check_assembly.o\
	tagU2.o\
	mess.o\
	find_oligo.o\
	copy_db.o \
	contig_order.o\
	clip.o\
	notes.o\
	consistency_display.o\
	consistency_canvas_box.o\
	confidence_graph.o\
	reading_coverage.o\
	readpair_coverage.o\
	strand_coverage.o\
	find_fragments.o\
	vseqs.o\
	$(GAP4_LEGACY)\
	f2c.o \
	shuffle_pads.o \
	auto_break.o

GAPOBJS=\
	$(COMMONOBJS)\
	$(GAPDB)\
	$(TCLBITS)\
	$(CONSEN)\
	$(GAPFUNCS)\
	$(EDITOR)

GAP4_LIBS=\
	$(IOLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(BIOLIMS_LIB) \
	$(BIOGAP_LIB) \
	$(G_LIB) \
	$(MUT_LIB) \
	$(MISC_LIB) \
	$(P3_LIB)

GAP4SH_LIBS=\
	$(IOLIB_LIB) \
	$(G_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(MISC_LIB) \
	$(MUT_LIB) \
	$(P3_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(GAPOBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LINK_PATHFLAG)$(L) $(GAPOBJS) $(GAP4_LIBS) $(CXX_DEP) $(SHLIB_DEP) $(F77_DEP)

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(GAPOBJS)
	$(MKDEFL) $@ $(GAPOBJS)

COPYDBOBJS=\
        $(GAPDB_LOW) \
        $(GAPDB_MID) \
        copy_db.o \
        copy_db_main.o \
        text-io-reg.o\
	io_utils.o \
	io_handle.o

CPLIB=\
	$(G_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB) \
	$(TCL_LIB)

copy_db:   $(COPYDBOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(COPYDBOBJS) $(CPLIB) $(LIBSC)

THRASHOBJS=\
        $(GAPDB_LOW) \
        $(GAPDB_MID) \
        gap-thrash3.o \
        text-io-reg.o

THRASHBOBJS=\
        $(GAPDB_LOW) \
        $(GAPDB_MID) \
        gap-thrash2bug.o \
        text-io-reg.o

TLIB=\
	$(G_LIB) \
	$(TEXTUTILS_LIB) \
	./libmisc.a

thrash:   $(THRASHOBJS)
	$(CLD) -o $@ $(THRASHOBJS) $(TLIB) $(LIBSC) -ldmalloc
thrashb:   $(THRASHBOBJS)
	$(CLD) -o $@ $(THRASHBOBJS) $(TLIB) $(LIBSC) -ldmalloc

DEPEND_OBJ = $(GAPOBJS) $(COPYDBOBJS)

S=$(SRCROOT)/gap4
install:
	cp copy_db $(INSTALLBIN)
	-cp $(S)/gap4 $(INSTALLBIN)
	-mkdir $(INSTALLTCL)/gap
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/gap
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/gaprc $(S)/gaprc_menu_full $(S)/gaprc_menu_mito $(INSTALLETC)
	cp $(S)/GTAGDB $(S)/NOTEDB $(INSTALLETC)
	cp $(S)/cap2rc $(S)/cap3rc $(S)/fak2rc $(S)/phraprc $(INSTALLETC)
	cp $(S)/*.bitmap $(S)/*.gif $(INSTALLETC)

distsrc: distsrc_dirs
	-cp -R *.[chf] *.tcl tclIndex Makefile gap4 gap4viewer gap4.bat \
		dependencies $(DIRNAME)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
