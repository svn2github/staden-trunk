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

install:
	cp copy_db $(INSTALLBIN)
	cp $(S)/gap4 $(INSTALLBIN)
	-mkdir $(INSTALLTCL)/gap
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/gap
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/gaprc $(S)/gaprc_menu_full $(S)/gaprc_menu_mito $(INSTALLETC)
	cp $(S)/GTAGDB $(S)/NOTEDB $(INSTALLETC)
	cp $(S)/cap2rc $(S)/cap3rc $(S)/fak2rc $(S)/phraprc $(INSTALLETC)
	cp $(S)/*.bitmap $(S)/*.gif $(INSTALLETC)

distsrc: distsrc_dirs
	cp $(S)/*.[chf] $(S)/*.tcl $(S)/tclIndex $(S)/Makefile $(S)/gap4 \
	   $(S)/gap4.bat $(S)/GTAGDB $(S)/NOTEDB $(S)/*.bitmap $(S)/*.gif \
	   $(S)/gaprc $(S)/gaprc_menu_full $(S)/gaprc_menu_mito \
	   $(S)/cap2rc $(S)/cap3rc $(S)/fak2rc $(S)/phraprc $(DIRNAME)

# DO NOT DELETE THIS LINE -- make depend depends on it.

IO.o: $(SRCROOT)/Misc/FtoC.h
IO.o: $(SRCROOT)/Misc/array.h
IO.o: $(SRCROOT)/Misc/bitmap.h
IO.o: $(SRCROOT)/Misc/dstring.h
IO.o: $(SRCROOT)/Misc/misc.h
IO.o: $(SRCROOT)/Misc/os.h
IO.o: $(SRCROOT)/Misc/xalloc.h
IO.o: $(SRCROOT)/Misc/xerror.h
IO.o: $(SRCROOT)/g/freetree.h
IO.o: $(SRCROOT)/g/g-defs.h
IO.o: $(SRCROOT)/g/g-error.h
IO.o: $(SRCROOT)/g/g-filedefs.h
IO.o: $(SRCROOT)/g/g-os.h
IO.o: $(SRCROOT)/g/g-struct.h
IO.o: $(SRCROOT)/primer3/src/dpal.h
IO.o: $(SRCROOT)/primer3/src/primer3.h
IO.o: $(PWD)/staden_config.h
IO.o: $(SRCROOT)/tk_utils/cli_arg.h
IO.o: $(SRCROOT)/tk_utils/postscript.h
IO.o: $(SRCROOT)/tk_utils/sheet.h
IO.o: $(SRCROOT)/tk_utils/text_output.h
IO.o: $(SRCROOT)/tk_utils/tkSheet.h
IO.o: $(SRCROOT)/tk_utils/tkSheet_common.h
IO.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
IO.o: $(SRCROOT)/tk_utils/tkTrace.h
IO.o: actf.h edUtils.h qual.h
IO.o: edStructs.h fort.h
IO.o: fortran.h
IO.o: gap-create.h gap-io.h
IO.o: gap-error.h IO.h
IO.o: gap-if.h
IO.o: gap-init.h gap-dbstruct.h
IO.o: io_handle.h io_utils.h
IO.o: notes.h
IO.o: primlib.h
IO.o: tagUtils.h tagDefs.h
IO.o: template.h list.h undo.h
IO.o: tkEdNames.h io-reg.h
IO.o: tkEditor.h
IO.o: tman_display.h
IO2.o: $(SRCROOT)/Misc/FtoC.h
IO2.o: $(SRCROOT)/Misc/array.h
IO2.o: $(SRCROOT)/Misc/bitmap.h
IO2.o: $(SRCROOT)/Misc/dstring.h
IO2.o: $(SRCROOT)/Misc/misc.h
IO2.o: $(SRCROOT)/Misc/os.h
IO2.o: $(SRCROOT)/Misc/xalloc.h
IO2.o: $(SRCROOT)/Misc/xerror.h
IO2.o: $(SRCROOT)/g/freetree.h
IO2.o: $(SRCROOT)/g/g-defs.h
IO2.o: $(SRCROOT)/g/g-error.h
IO2.o: $(SRCROOT)/g/g-filedefs.h
IO2.o: $(SRCROOT)/g/g-os.h
IO2.o: $(SRCROOT)/g/g-struct.h
IO2.o: $(SRCROOT)/primer3/src/dpal.h
IO2.o: $(SRCROOT)/primer3/src/primer3.h
IO2.o: $(SRCROOT)/seq_utils/dna_utils.h
IO2.o: $(PWD)/staden_config.h
IO2.o: $(SRCROOT)/tk_utils/canvas_box.h
IO2.o: $(SRCROOT)/tk_utils/cli_arg.h
IO2.o: $(SRCROOT)/tk_utils/postscript.h
IO2.o: $(SRCROOT)/tk_utils/sheet.h
IO2.o: $(SRCROOT)/tk_utils/text_output.h
IO2.o: $(SRCROOT)/tk_utils/tkSheet.h
IO2.o: $(SRCROOT)/tk_utils/tkSheet_common.h
IO2.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
IO2.o: $(SRCROOT)/tk_utils/tkTrace.h
IO2.o: IO.h
IO2.o: IO2.h tagUtils.h
IO2.o: clones.h seqInfo.h notes.h
IO2.o: contig_selector.h
IO2.o: cs-object.h
IO2.o: edStructs.h
IO2.o: edUtils.h qual.h
IO2.o: fort.h tkEditor.h
IO2.o: gap-create.h
IO2.o: gap-dbstruct.h
IO2.o: gap-if.h
IO2.o: io_handle.h io_utils.h
IO2.o: primlib.h
IO2.o: tagDefs.h fortran.h
IO2.o: template.h list.h undo.h
IO2.o: tkEdNames.h io-reg.h
IO2.o: tman_display.h
IO3.o: $(SRCROOT)/Misc/FtoC.h
IO3.o: $(SRCROOT)/Misc/array.h
IO3.o: $(SRCROOT)/Misc/bitmap.h
IO3.o: $(SRCROOT)/Misc/dstring.h
IO3.o: $(SRCROOT)/Misc/misc.h
IO3.o: $(SRCROOT)/Misc/os.h
IO3.o: $(SRCROOT)/Misc/xalloc.h
IO3.o: $(SRCROOT)/Misc/xerror.h
IO3.o: $(SRCROOT)/g/freetree.h
IO3.o: $(SRCROOT)/g/g-defs.h
IO3.o: $(SRCROOT)/g/g-error.h
IO3.o: $(SRCROOT)/g/g-filedefs.h
IO3.o: $(SRCROOT)/g/g-os.h
IO3.o: $(SRCROOT)/g/g-struct.h
IO3.o: $(SRCROOT)/primer3/src/dpal.h
IO3.o: $(SRCROOT)/primer3/src/primer3.h
IO3.o: $(PWD)/staden_config.h
IO3.o: $(SRCROOT)/tk_utils/cli_arg.h
IO3.o: $(SRCROOT)/tk_utils/postscript.h
IO3.o: $(SRCROOT)/tk_utils/sheet.h
IO3.o: $(SRCROOT)/tk_utils/tkSheet.h
IO3.o: $(SRCROOT)/tk_utils/tkSheet_common.h
IO3.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
IO3.o: $(SRCROOT)/tk_utils/tkTrace.h
IO3.o: IO.h
IO3.o: edStructs.h fort.h
IO3.o: edUtils.h qual.h
IO3.o: fortran.h
IO3.o: gap-dbstruct.h
IO3.o: gap-error.h
IO3.o: gap-if.h
IO3.o: io_handle.h io_utils.h
IO3.o: primlib.h
IO3.o: tagUtils.h tagDefs.h
IO3.o: template.h list.h undo.h
IO3.o: tkEdNames.h io-reg.h
IO3.o: tkEditor.h
IO3.o: tman_display.h
actf.o: $(SRCROOT)/Misc/array.h
actf.o: $(SRCROOT)/Misc/bitmap.h
actf.o: $(SRCROOT)/Misc/misc.h
actf.o: $(SRCROOT)/Misc/os.h
actf.o: $(SRCROOT)/Misc/xalloc.h
actf.o: $(SRCROOT)/Misc/xerror.h
actf.o: $(SRCROOT)/g/freetree.h
actf.o: $(SRCROOT)/g/g-defs.h
actf.o: $(SRCROOT)/g/g-error.h
actf.o: $(SRCROOT)/g/g-filedefs.h
actf.o: $(SRCROOT)/g/g-os.h
actf.o: $(SRCROOT)/g/g-struct.h
actf.o: $(PWD)/staden_config.h
actf.o: $(SRCROOT)/tk_utils/text_output.h
actf.o: IO.h
actf.o: actf.h
actf.o: fort.h
actf.o: gap-dbstruct.h
actf.o: gap-if.h
actf.o: io_handle.h
actf.o: io_utils.h
active_tags.o: $(SRCROOT)/Misc/array.h
active_tags.o: $(SRCROOT)/Misc/bitmap.h
active_tags.o: $(SRCROOT)/Misc/dstring.h
active_tags.o: $(SRCROOT)/Misc/misc.h
active_tags.o: $(SRCROOT)/Misc/os.h
active_tags.o: $(SRCROOT)/Misc/xalloc.h
active_tags.o: $(SRCROOT)/Misc/xerror.h
active_tags.o: $(SRCROOT)/g/freetree.h
active_tags.o: $(SRCROOT)/g/g-defs.h
active_tags.o: $(SRCROOT)/g/g-error.h
active_tags.o: $(SRCROOT)/g/g-filedefs.h
active_tags.o: $(SRCROOT)/g/g-os.h
active_tags.o: $(SRCROOT)/g/g-struct.h
active_tags.o: $(SRCROOT)/primer3/src/dpal.h
active_tags.o: $(SRCROOT)/primer3/src/primer3.h
active_tags.o: $(PWD)/staden_config.h
active_tags.o: $(SRCROOT)/tk_utils/cli_arg.h
active_tags.o: $(SRCROOT)/tk_utils/postscript.h
active_tags.o: $(SRCROOT)/tk_utils/sheet.h
active_tags.o: $(SRCROOT)/tk_utils/tkSheet.h
active_tags.o: $(SRCROOT)/tk_utils/tkSheet_common.h
active_tags.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
active_tags.o: $(SRCROOT)/tk_utils/tkTrace.h
active_tags.o: IO.h
active_tags.o: contigEditor.h
active_tags.o: fort.h edStructs.h
active_tags.o: fortran.h
active_tags.o: gap-dbstruct.h
active_tags.o: gap-if.h
active_tags.o: io-reg.h primlib.h
active_tags.o: io_handle.h
active_tags.o: io_utils.h
active_tags.o: list_proc.h
active_tags.o: qual.h tagUtils.h
active_tags.o: tagDefs.h
active_tags.o: tagdb.h
active_tags.o: template.h list.h
active_tags.o: tkEdNames.h
active_tags.o: tkEditor.h
active_tags.o: tman_display.h
active_tags.o: undo.h edUtils.h
alter_rel.o: $(SRCROOT)/Misc/array.h
alter_rel.o: $(SRCROOT)/Misc/bitmap.h
alter_rel.o: $(SRCROOT)/Misc/dstring.h
alter_rel.o: $(SRCROOT)/Misc/misc.h
alter_rel.o: $(SRCROOT)/Misc/os.h
alter_rel.o: $(SRCROOT)/Misc/xalloc.h
alter_rel.o: $(SRCROOT)/Misc/xerror.h
alter_rel.o: $(SRCROOT)/g/freetree.h
alter_rel.o: $(SRCROOT)/g/g-defs.h
alter_rel.o: $(SRCROOT)/g/g-error.h
alter_rel.o: $(SRCROOT)/g/g-filedefs.h
alter_rel.o: $(SRCROOT)/g/g-os.h
alter_rel.o: $(SRCROOT)/g/g-struct.h
alter_rel.o: $(SRCROOT)/primer3/src/dpal.h
alter_rel.o: $(SRCROOT)/primer3/src/primer3.h
alter_rel.o: $(PWD)/staden_config.h
alter_rel.o: $(SRCROOT)/tk_utils/cli_arg.h
alter_rel.o: $(SRCROOT)/tk_utils/postscript.h
alter_rel.o: $(SRCROOT)/tk_utils/sheet.h
alter_rel.o: $(SRCROOT)/tk_utils/tkSheet.h
alter_rel.o: $(SRCROOT)/tk_utils/tkSheet_common.h
alter_rel.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
alter_rel.o: $(SRCROOT)/tk_utils/tkTrace.h
alter_rel.o: IO.h
alter_rel.o: dis_readings.h
alter_rel.o: edStructs.h
alter_rel.o: gap-dbstruct.h
alter_rel.o: gap-error.h fort.h
alter_rel.o: gap-if.h
alter_rel.o: io-reg.h alter_rel.h
alter_rel.o: io_handle.h
alter_rel.o: io_utils.h
alter_rel.o: primlib.h
alter_rel.o: qual.h
alter_rel.o: tagDefs.h fortran.h
alter_rel.o: tagUtils.h
alter_rel.o: template.h list.h
alter_rel.o: tkEdNames.h
alter_rel.o: tkEditor.h
alter_rel.o: tman_display.h
alter_rel.o: undo.h edUtils.h
assemble_direct.o: $(SRCROOT)/Misc/array.h
assemble_direct.o: $(SRCROOT)/Misc/bitmap.h
assemble_direct.o: $(SRCROOT)/Misc/dstring.h
assemble_direct.o: $(SRCROOT)/Misc/misc.h
assemble_direct.o: $(SRCROOT)/Misc/os.h
assemble_direct.o: $(SRCROOT)/Misc/xalloc.h
assemble_direct.o: $(SRCROOT)/Misc/xerror.h
assemble_direct.o: $(SRCROOT)/g/freetree.h
assemble_direct.o: $(SRCROOT)/g/g-defs.h
assemble_direct.o: $(SRCROOT)/g/g-error.h
assemble_direct.o: $(SRCROOT)/g/g-filedefs.h
assemble_direct.o: $(SRCROOT)/g/g-os.h
assemble_direct.o: $(SRCROOT)/g/g-struct.h
assemble_direct.o: $(SRCROOT)/primer3/src/dpal.h
assemble_direct.o: $(SRCROOT)/primer3/src/primer3.h
assemble_direct.o: $(SRCROOT)/seq_utils/align.h
assemble_direct.o: $(SRCROOT)/seq_utils/align_lib_old.h
assemble_direct.o: $(SRCROOT)/seq_utils/dna_utils.h
assemble_direct.o: $(PWD)/staden_config.h
assemble_direct.o: $(SRCROOT)/tk_utils/cli_arg.h
assemble_direct.o: $(SRCROOT)/tk_utils/postscript.h
assemble_direct.o: $(SRCROOT)/tk_utils/sheet.h
assemble_direct.o: $(SRCROOT)/tk_utils/tcl_utils.h
assemble_direct.o: $(SRCROOT)/tk_utils/text_output.h
assemble_direct.o: $(SRCROOT)/tk_utils/tkSheet.h
assemble_direct.o: $(SRCROOT)/tk_utils/tkSheet_common.h
assemble_direct.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
assemble_direct.o: $(SRCROOT)/tk_utils/tkTrace.h
assemble_direct.o: IO.h
assemble_direct.o: IO2.h
assemble_direct.o: assemble_direct.h
assemble_direct.o: clones.h
assemble_direct.o: consen.h
assemble_direct.o: dis_readings.h
assemble_direct.o: edStructs.h
assemble_direct.o: edUtils.h
assemble_direct.o: fort.h
assemble_direct.o: fortran.h
assemble_direct.o: gap-dbstruct.h
assemble_direct.o: gap-if.h
assemble_direct.o: gap_globals.h
assemble_direct.o: io-reg.h
assemble_direct.o: io_handle.h
assemble_direct.o: io_utils.h
assemble_direct.o: list.h undo.h
assemble_direct.o: list_proc.h
assemble_direct.o: notes.h
assemble_direct.o: primlib.h
assemble_direct.o: qual.h
assemble_direct.o: seqInfo.h
assemble_direct.o: tagDefs.h
assemble_direct.o: tagUtils.h
assemble_direct.o: template.h
assemble_direct.o: tkEdNames.h
assemble_direct.o: tkEditor.h
assemble_direct.o: tman_display.h
auto_assemble.o: $(SRCROOT)/Misc/array.h
auto_assemble.o: $(SRCROOT)/Misc/bitmap.h
auto_assemble.o: $(SRCROOT)/Misc/misc.h
auto_assemble.o: $(SRCROOT)/Misc/os.h
auto_assemble.o: $(SRCROOT)/Misc/xalloc.h
auto_assemble.o: $(SRCROOT)/Misc/xerror.h
auto_assemble.o: $(SRCROOT)/g/freetree.h
auto_assemble.o: $(SRCROOT)/g/g-defs.h
auto_assemble.o: $(SRCROOT)/g/g-error.h
auto_assemble.o: $(SRCROOT)/g/g-filedefs.h
auto_assemble.o: $(SRCROOT)/g/g-os.h
auto_assemble.o: $(SRCROOT)/g/g-struct.h
auto_assemble.o: $(PWD)/staden_config.h
auto_assemble.o: $(SRCROOT)/tk_utils/tcl_utils.h
auto_assemble.o: IO.h
auto_assemble.o: complement.h
auto_assemble.o: fort.h
auto_assemble.o: gap-dbstruct.h
auto_assemble.o: gap-if.h
auto_assemble.o: gap_globals.h
auto_assemble.o: io_handle.h
auto_assemble.o: io_utils.h
auto_assemble.o: list_proc.h
auto_break.o: $(SRCROOT)/Misc/array.h
auto_break.o: $(SRCROOT)/Misc/bitmap.h
auto_break.o: $(SRCROOT)/Misc/dstring.h
auto_break.o: $(SRCROOT)/Misc/misc.h
auto_break.o: $(SRCROOT)/Misc/os.h
auto_break.o: $(SRCROOT)/Misc/xalloc.h
auto_break.o: $(SRCROOT)/Misc/xerror.h
auto_break.o: $(SRCROOT)/g/freetree.h
auto_break.o: $(SRCROOT)/g/g-defs.h
auto_break.o: $(SRCROOT)/g/g-error.h
auto_break.o: $(SRCROOT)/g/g-filedefs.h
auto_break.o: $(SRCROOT)/g/g-os.h
auto_break.o: $(SRCROOT)/g/g-struct.h
auto_break.o: $(SRCROOT)/seq_utils/filter_words.h
auto_break.o: $(PWD)/staden_config.h
auto_break.o: $(SRCROOT)/tk_utils/tcl_utils.h
auto_break.o: $(SRCROOT)/tk_utils/text_output.h
auto_break.o: IO.h
auto_break.o: auto_break.h
auto_break.o: gap-dbstruct.h
auto_break.o: gap-if.h
auto_break.o: gap_globals.h
auto_break.o: io_handle.h
auto_break.o: io_utils.h
auto_break.o: qual.h qualIO.h
auto_break.o: template.h list.h
break_contig.o: $(SRCROOT)/Misc/array.h
break_contig.o: $(SRCROOT)/Misc/bitmap.h
break_contig.o: $(SRCROOT)/Misc/misc.h
break_contig.o: $(SRCROOT)/Misc/os.h
break_contig.o: $(SRCROOT)/Misc/xalloc.h
break_contig.o: $(SRCROOT)/Misc/xerror.h
break_contig.o: $(SRCROOT)/g/freetree.h
break_contig.o: $(SRCROOT)/g/g-defs.h
break_contig.o: $(SRCROOT)/g/g-error.h
break_contig.o: $(SRCROOT)/g/g-filedefs.h
break_contig.o: $(SRCROOT)/g/g-os.h
break_contig.o: $(SRCROOT)/g/g-struct.h
break_contig.o: $(PWD)/staden_config.h
break_contig.o: IO.h
break_contig.o: dis_readings.h
break_contig.o: gap-dbstruct.h
break_contig.o: gap-if.h
break_contig.o: io-reg.h fort.h
break_contig.o: io_handle.h
break_contig.o: io_utils.h
bubbl3.o: $(SRCROOT)/Misc/os.h
bubbl3.o: $(PWD)/staden_config.h
check_assembly.o: $(SRCROOT)/Misc/array.h
check_assembly.o: $(SRCROOT)/Misc/bitmap.h
check_assembly.o: $(SRCROOT)/Misc/dstring.h
check_assembly.o: $(SRCROOT)/Misc/misc.h
check_assembly.o: $(SRCROOT)/Misc/os.h
check_assembly.o: $(SRCROOT)/Misc/xalloc.h
check_assembly.o: $(SRCROOT)/Misc/xerror.h
check_assembly.o: $(SRCROOT)/g/freetree.h
check_assembly.o: $(SRCROOT)/g/g-defs.h
check_assembly.o: $(SRCROOT)/g/g-error.h
check_assembly.o: $(SRCROOT)/g/g-filedefs.h
check_assembly.o: $(SRCROOT)/g/g-os.h
check_assembly.o: $(SRCROOT)/g/g-struct.h
check_assembly.o: $(SRCROOT)/primer3/src/dpal.h
check_assembly.o: $(SRCROOT)/primer3/src/primer3.h
check_assembly.o: $(SRCROOT)/seq_utils/align.h
check_assembly.o: $(SRCROOT)/seq_utils/align_lib_old.h
check_assembly.o: $(SRCROOT)/seq_utils/dna_utils.h
check_assembly.o: $(PWD)/staden_config.h
check_assembly.o: $(SRCROOT)/tk_utils/canvas_box.h
check_assembly.o: $(SRCROOT)/tk_utils/cli_arg.h
check_assembly.o: $(SRCROOT)/tk_utils/postscript.h
check_assembly.o: $(SRCROOT)/tk_utils/sheet.h
check_assembly.o: $(SRCROOT)/tk_utils/tclXkeylist.h
check_assembly.o: $(SRCROOT)/tk_utils/tcl_utils.h
check_assembly.o: $(SRCROOT)/tk_utils/text_output.h
check_assembly.o: $(SRCROOT)/tk_utils/tkSheet.h
check_assembly.o: $(SRCROOT)/tk_utils/tkSheet_common.h
check_assembly.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
check_assembly.o: $(SRCROOT)/tk_utils/tkTrace.h
check_assembly.o: IO.h
check_assembly.o: consen.h
check_assembly.o: contigEditor.h
check_assembly.o: contig_selector.h
check_assembly.o: cs-object.h
check_assembly.o: edStructs.h
check_assembly.o: edUtils.h
check_assembly.o: fort.h
check_assembly.o: fortran.h
check_assembly.o: gap-dbstruct.h
check_assembly.o: gap-if.h
check_assembly.o: gap_globals.h
check_assembly.o: io-reg.h
check_assembly.o: io_handle.h
check_assembly.o: io_utils.h
check_assembly.o: list.h undo.h
check_assembly.o: primlib.h
check_assembly.o: qual.h
check_assembly.o: tagDefs.h
check_assembly.o: tagUtils.h
check_assembly.o: template.h
check_assembly.o: tkEdNames.h
check_assembly.o: tkEditor.h
check_assembly.o: tman_display.h
clip.o: $(SRCROOT)/Misc/array.h
clip.o: $(SRCROOT)/Misc/bitmap.h
clip.o: $(SRCROOT)/Misc/dstring.h
clip.o: $(SRCROOT)/Misc/misc.h
clip.o: $(SRCROOT)/Misc/os.h
clip.o: $(SRCROOT)/Misc/xalloc.h
clip.o: $(SRCROOT)/Misc/xerror.h
clip.o: $(SRCROOT)/g/freetree.h
clip.o: $(SRCROOT)/g/g-defs.h
clip.o: $(SRCROOT)/g/g-error.h
clip.o: $(SRCROOT)/g/g-filedefs.h
clip.o: $(SRCROOT)/g/g-os.h
clip.o: $(SRCROOT)/g/g-struct.h
clip.o: $(SRCROOT)/primer3/src/dpal.h
clip.o: $(SRCROOT)/primer3/src/primer3.h
clip.o: $(SRCROOT)/seq_utils/dna_utils.h
clip.o: $(PWD)/staden_config.h
clip.o: $(SRCROOT)/tk_utils/cli_arg.h
clip.o: $(SRCROOT)/tk_utils/postscript.h
clip.o: $(SRCROOT)/tk_utils/sheet.h
clip.o: $(SRCROOT)/tk_utils/tcl_utils.h
clip.o: $(SRCROOT)/tk_utils/tkSheet.h
clip.o: $(SRCROOT)/tk_utils/tkSheet_common.h
clip.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
clip.o: $(SRCROOT)/tk_utils/tkTrace.h
clip.o: clip.h IO.h
clip.o: dis_readings.h
clip.o: edUtils.h tman_display.h
clip.o: fort.h tkEditor.h
clip.o: gap-dbstruct.h
clip.o: gap-if.h
clip.o: gap_globals.h
clip.o: io_handle.h io_utils.h
clip.o: primlib.h
clip.o: qual.h
clip.o: tagDefs.h fortran.h
clip.o: tagUtils.h edStructs.h
clip.o: template.h list.h undo.h
clip.o: tkEdNames.h io-reg.h
clones.o: $(SRCROOT)/Misc/array.h
clones.o: $(SRCROOT)/Misc/bitmap.h
clones.o: $(SRCROOT)/Misc/dstring.h
clones.o: $(SRCROOT)/Misc/os.h
clones.o: $(SRCROOT)/Misc/xalloc.h
clones.o: $(SRCROOT)/Misc/xerror.h
clones.o: $(SRCROOT)/g/freetree.h
clones.o: $(SRCROOT)/g/g-defs.h
clones.o: $(SRCROOT)/g/g-error.h
clones.o: $(SRCROOT)/g/g-filedefs.h
clones.o: $(SRCROOT)/g/g-os.h
clones.o: $(SRCROOT)/g/g-struct.h
clones.o: $(SRCROOT)/primer3/src/dpal.h
clones.o: $(SRCROOT)/primer3/src/primer3.h
clones.o: $(PWD)/staden_config.h
clones.o: $(SRCROOT)/tk_utils/cli_arg.h
clones.o: $(SRCROOT)/tk_utils/postscript.h
clones.o: $(SRCROOT)/tk_utils/sheet.h
clones.o: $(SRCROOT)/tk_utils/tkSheet.h
clones.o: $(SRCROOT)/tk_utils/tkSheet_common.h
clones.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
clones.o: $(SRCROOT)/tk_utils/tkTrace.h
clones.o: clones.h notes.h
clones.o: fort.h tkEditor.h
clones.o: gap-dbstruct.h
clones.o: gap-defaults.h IO.h
clones.o: gap-if.h
clones.o: io_handle.h io_utils.h
clones.o: primlib.h
clones.o: seqInfo.h
clones.o: tagDefs.h fortran.h
clones.o: tagUtils.h edStructs.h
clones.o: template.h list.h
clones.o: tkEdNames.h io-reg.h
clones.o: tman_display.h
clones.o: undo.h edUtils.h qual.h
complement.o: $(SRCROOT)/Misc/array.h
complement.o: $(SRCROOT)/Misc/bitmap.h
complement.o: $(SRCROOT)/Misc/dstring.h
complement.o: $(SRCROOT)/Misc/misc.h
complement.o: $(SRCROOT)/Misc/os.h
complement.o: $(SRCROOT)/Misc/xalloc.h
complement.o: $(SRCROOT)/Misc/xerror.h
complement.o: $(SRCROOT)/g/freetree.h
complement.o: $(SRCROOT)/g/g-defs.h
complement.o: $(SRCROOT)/g/g-error.h
complement.o: $(SRCROOT)/g/g-filedefs.h
complement.o: $(SRCROOT)/g/g-os.h
complement.o: $(SRCROOT)/g/g-struct.h
complement.o: $(SRCROOT)/primer3/src/dpal.h
complement.o: $(SRCROOT)/primer3/src/primer3.h
complement.o: $(SRCROOT)/seq_utils/dna_utils.h
complement.o: $(PWD)/staden_config.h
complement.o: $(SRCROOT)/tk_utils/cli_arg.h
complement.o: $(SRCROOT)/tk_utils/postscript.h
complement.o: $(SRCROOT)/tk_utils/sheet.h
complement.o: $(SRCROOT)/tk_utils/tcl_utils.h
complement.o: $(SRCROOT)/tk_utils/tkSheet.h
complement.o: $(SRCROOT)/tk_utils/tkSheet_common.h
complement.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
complement.o: $(SRCROOT)/tk_utils/tkTrace.h
complement.o: IO.h
complement.o: complement.h
complement.o: edStructs.h
complement.o: fort.h tagUtils.h
complement.o: gap-dbstruct.h
complement.o: gap-if.h
complement.o: gap_cli_arg.h
complement.o: io-reg.h
complement.o: io_handle.h
complement.o: io_utils.h
complement.o: list_proc.h
complement.o: primlib.h
complement.o: qual.h
complement.o: tagDefs.h fortran.h
complement.o: template.h list.h
complement.o: tkEdNames.h
complement.o: tkEditor.h
complement.o: tman_display.h
complement.o: undo.h edUtils.h
confidence_graph.o: $(SRCROOT)/Misc/array.h
confidence_graph.o: $(SRCROOT)/Misc/bitmap.h
confidence_graph.o: $(SRCROOT)/Misc/misc.h
confidence_graph.o: $(SRCROOT)/Misc/os.h
confidence_graph.o: $(SRCROOT)/Misc/xalloc.h
confidence_graph.o: $(SRCROOT)/Misc/xerror.h
confidence_graph.o: $(SRCROOT)/g/freetree.h
confidence_graph.o: $(SRCROOT)/g/g-defs.h
confidence_graph.o: $(SRCROOT)/g/g-error.h
confidence_graph.o: $(SRCROOT)/g/g-filedefs.h
confidence_graph.o: $(SRCROOT)/g/g-os.h
confidence_graph.o: $(SRCROOT)/g/g-struct.h
confidence_graph.o: $(PWD)/staden_config.h
confidence_graph.o: $(SRCROOT)/tk_utils/canvas_box.h
confidence_graph.o: $(SRCROOT)/tk_utils/tclXkeylist.h
confidence_graph.o: $(SRCROOT)/tk_utils/tcl_utils.h
confidence_graph.o: $(SRCROOT)/tk_utils/text_output.h
confidence_graph.o: IO.h
confidence_graph.o: confidence_graph.h
confidence_graph.o: consen.h
confidence_graph.o: consistency_display.h
confidence_graph.o: gap-dbstruct.h
confidence_graph.o: gap-if.h
confidence_graph.o: gap_canvas_box.h
confidence_graph.o: gap_globals.h
confidence_graph.o: io-reg.h
confidence_graph.o: io_handle.h
confidence_graph.o: io_utils.h
confidence_graph.o: list.h
confidence_graph.o: newgap_cmds.h
confidence_graph.o: qual.h
confidence_graph.o: ruler_display.h
consen.o: $(SRCROOT)/Misc/FtoC.h
consen.o: $(SRCROOT)/Misc/array.h
consen.o: $(SRCROOT)/Misc/bitmap.h
consen.o: $(SRCROOT)/Misc/dstring.h
consen.o: $(SRCROOT)/Misc/misc.h
consen.o: $(SRCROOT)/Misc/os.h
consen.o: $(SRCROOT)/Misc/xalloc.h
consen.o: $(SRCROOT)/Misc/xerror.h
consen.o: $(SRCROOT)/g/freetree.h
consen.o: $(SRCROOT)/g/g-defs.h
consen.o: $(SRCROOT)/g/g-error.h
consen.o: $(SRCROOT)/g/g-filedefs.h
consen.o: $(SRCROOT)/g/g-os.h
consen.o: $(SRCROOT)/g/g-struct.h
consen.o: $(SRCROOT)/primer3/src/dpal.h
consen.o: $(SRCROOT)/primer3/src/primer3.h
consen.o: $(SRCROOT)/seq_utils/align_lib.h
consen.o: $(SRCROOT)/seq_utils/dna_utils.h
consen.o: $(PWD)/staden_config.h
consen.o: $(SRCROOT)/tk_utils/cli_arg.h
consen.o: $(SRCROOT)/tk_utils/postscript.h
consen.o: $(SRCROOT)/tk_utils/sheet.h
consen.o: $(SRCROOT)/tk_utils/tcl_utils.h
consen.o: $(SRCROOT)/tk_utils/tkSheet.h
consen.o: $(SRCROOT)/tk_utils/tkSheet_common.h
consen.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
consen.o: $(SRCROOT)/tk_utils/tkTrace.h
consen.o: consen.h gap_globals.h
consen.o: edStructs.h tkEditor.h
consen.o: edUtils.h qual.h
consen.o: gap-dbstruct.h
consen.o: gap-if.h
consen.o: io_handle.h fort.h
consen.o: io_utils.h IO.h
consen.o: list_proc.h extract.h
consen.o: primlib.h
consen.o: tagDefs.h fortran.h
consen.o: template.h list.h
consen.o: tkEdNames.h io-reg.h
consen.o: tman_display.h
consen.o: undo.h tagUtils.h
consistency_canvas_box.o: $(SRCROOT)/Misc/array.h
consistency_canvas_box.o: $(SRCROOT)/Misc/bitmap.h
consistency_canvas_box.o: $(SRCROOT)/Misc/os.h
consistency_canvas_box.o: $(SRCROOT)/Misc/xerror.h
consistency_canvas_box.o: $(SRCROOT)/g/freetree.h
consistency_canvas_box.o: $(SRCROOT)/g/g-defs.h
consistency_canvas_box.o: $(SRCROOT)/g/g-error.h
consistency_canvas_box.o: $(SRCROOT)/g/g-filedefs.h
consistency_canvas_box.o: $(SRCROOT)/g/g-os.h
consistency_canvas_box.o: $(SRCROOT)/g/g-struct.h
consistency_canvas_box.o: $(PWD)/staden_config.h
consistency_canvas_box.o: $(SRCROOT)/tk_utils/canvas_box.h
consistency_canvas_box.o: IO.h
consistency_canvas_box.o: consistency_display.h
consistency_canvas_box.o: gap-dbstruct.h
consistency_canvas_box.o: gap-if.h
consistency_canvas_box.o: gap_canvas_box.h
consistency_canvas_box.o: io-reg.h
consistency_canvas_box.o: io_handle.h
consistency_canvas_box.o: io_utils.h
consistency_display.o: $(SRCROOT)/Misc/array.h
consistency_display.o: $(SRCROOT)/Misc/bitmap.h
consistency_display.o: $(SRCROOT)/Misc/misc.h
consistency_display.o: $(SRCROOT)/Misc/os.h
consistency_display.o: $(SRCROOT)/Misc/xalloc.h
consistency_display.o: $(SRCROOT)/Misc/xerror.h
consistency_display.o: $(SRCROOT)/g/freetree.h
consistency_display.o: $(SRCROOT)/g/g-defs.h
consistency_display.o: $(SRCROOT)/g/g-error.h
consistency_display.o: $(SRCROOT)/g/g-filedefs.h
consistency_display.o: $(SRCROOT)/g/g-os.h
consistency_display.o: $(SRCROOT)/g/g-struct.h
consistency_display.o: $(PWD)/staden_config.h
consistency_display.o: $(SRCROOT)/tk_utils/canvas_box.h
consistency_display.o: $(SRCROOT)/tk_utils/ruler_tick.h
consistency_display.o: $(SRCROOT)/tk_utils/tcl_utils.h
consistency_display.o: IO.h
consistency_display.o: consistency_canvas_box.h
consistency_display.o: consistency_display.h
consistency_display.o: gap-dbstruct.h
consistency_display.o: gap-if.h
consistency_display.o: gap_canvas_box.h
consistency_display.o: gap_globals.h
consistency_display.o: hash.h
consistency_display.o: io-reg.h
consistency_display.o: io_handle.h
consistency_display.o: io_utils.h
consistency_display.o: list.h
consistency_display.o: ruler_display.h
consistency_display.o: template.h
consistency_display.o: template_display.h
contigEditor.o: $(SRCROOT)/Misc/array.h
contigEditor.o: $(SRCROOT)/Misc/bitmap.h
contigEditor.o: $(SRCROOT)/Misc/dstring.h
contigEditor.o: $(SRCROOT)/Misc/locks.h
contigEditor.o: $(SRCROOT)/Misc/misc.h
contigEditor.o: $(SRCROOT)/Misc/os.h
contigEditor.o: $(SRCROOT)/Misc/xalloc.h
contigEditor.o: $(SRCROOT)/Misc/xerror.h
contigEditor.o: $(SRCROOT)/g/freetree.h
contigEditor.o: $(SRCROOT)/g/g-defs.h
contigEditor.o: $(SRCROOT)/g/g-error.h
contigEditor.o: $(SRCROOT)/g/g-filedefs.h
contigEditor.o: $(SRCROOT)/g/g-os.h
contigEditor.o: $(SRCROOT)/g/g-struct.h
contigEditor.o: $(SRCROOT)/primer3/src/dpal.h
contigEditor.o: $(SRCROOT)/primer3/src/primer3.h
contigEditor.o: $(PWD)/staden_config.h
contigEditor.o: $(SRCROOT)/tk_utils/cli_arg.h
contigEditor.o: $(SRCROOT)/tk_utils/postscript.h
contigEditor.o: $(SRCROOT)/tk_utils/sheet.h
contigEditor.o: $(SRCROOT)/tk_utils/tclXkeylist.h
contigEditor.o: $(SRCROOT)/tk_utils/tcl_utils.h
contigEditor.o: $(SRCROOT)/tk_utils/tkSheet.h
contigEditor.o: $(SRCROOT)/tk_utils/tkSheet_common.h
contigEditor.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
contigEditor.o: $(SRCROOT)/tk_utils/tkTrace.h
contigEditor.o: contigEditor.h
contigEditor.o: edStructs.h
contigEditor.o: fort.h tkEditor.h
contigEditor.o: fortran.h
contigEditor.o: gap-dbstruct.h
contigEditor.o: gap-if.h
contigEditor.o: gap_cli_arg.h
contigEditor.o: gap_globals.h
contigEditor.o: io-reg.h IO.h
contigEditor.o: io_handle.h
contigEditor.o: io_utils.h
contigEditor.o: list_proc.h
contigEditor.o: primlib.h
contigEditor.o: qual.h tagUtils.h
contigEditor.o: select.h
contigEditor.o: tagDefs.h
contigEditor.o: template.h list.h
contigEditor.o: tkEdNames.h
contigEditor.o: tman_display.h
contigEditor.o: tman_interface.h
contigEditor.o: undo.h edUtils.h
contig_order.o: $(SRCROOT)/Misc/array.h
contig_order.o: $(SRCROOT)/Misc/bitmap.h
contig_order.o: $(SRCROOT)/Misc/misc.h
contig_order.o: $(SRCROOT)/Misc/os.h
contig_order.o: $(SRCROOT)/Misc/xalloc.h
contig_order.o: $(SRCROOT)/Misc/xerror.h
contig_order.o: $(SRCROOT)/g/freetree.h
contig_order.o: $(SRCROOT)/g/g-defs.h
contig_order.o: $(SRCROOT)/g/g-error.h
contig_order.o: $(SRCROOT)/g/g-filedefs.h
contig_order.o: $(SRCROOT)/g/g-os.h
contig_order.o: $(SRCROOT)/g/g-struct.h
contig_order.o: $(PWD)/staden_config.h
contig_order.o: IO.h
contig_order.o: QueueInterface.h
contig_order.o: SeqQueueTypes.h
contig_order.o: gap-dbstruct.h
contig_order.o: gap-if.h
contig_order.o: io-reg.h
contig_order.o: io_handle.h
contig_order.o: io_utils.h
contig_order.o: template.h list.h
contig_selector.o: $(SRCROOT)/Misc/array.h
contig_selector.o: $(SRCROOT)/Misc/bitmap.h
contig_selector.o: $(SRCROOT)/Misc/dstring.h
contig_selector.o: $(SRCROOT)/Misc/misc.h
contig_selector.o: $(SRCROOT)/Misc/os.h
contig_selector.o: $(SRCROOT)/Misc/xalloc.h
contig_selector.o: $(SRCROOT)/Misc/xerror.h
contig_selector.o: $(SRCROOT)/g/freetree.h
contig_selector.o: $(SRCROOT)/g/g-defs.h
contig_selector.o: $(SRCROOT)/g/g-error.h
contig_selector.o: $(SRCROOT)/g/g-filedefs.h
contig_selector.o: $(SRCROOT)/g/g-os.h
contig_selector.o: $(SRCROOT)/g/g-struct.h
contig_selector.o: $(SRCROOT)/primer3/src/dpal.h
contig_selector.o: $(SRCROOT)/primer3/src/primer3.h
contig_selector.o: $(PWD)/staden_config.h
contig_selector.o: $(SRCROOT)/tk_utils/canvas_box.h
contig_selector.o: $(SRCROOT)/tk_utils/cli_arg.h
contig_selector.o: $(SRCROOT)/tk_utils/postscript.h
contig_selector.o: $(SRCROOT)/tk_utils/sheet.h
contig_selector.o: $(SRCROOT)/tk_utils/tclXkeylist.h
contig_selector.o: $(SRCROOT)/tk_utils/tcl_utils.h
contig_selector.o: $(SRCROOT)/tk_utils/text_output.h
contig_selector.o: $(SRCROOT)/tk_utils/tkSheet.h
contig_selector.o: $(SRCROOT)/tk_utils/tkSheet_common.h
contig_selector.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
contig_selector.o: $(SRCROOT)/tk_utils/tkTrace.h
contig_selector.o: IO.h
contig_selector.o: active_tags.h
contig_selector.o: complement.h
contig_selector.o: contigEditor.h
contig_selector.o: contig_selector.h
contig_selector.o: cs-object.h
contig_selector.o: edStructs.h
contig_selector.o: edUtils.h
contig_selector.o: fort.h
contig_selector.o: fortran.h
contig_selector.o: gap-dbstruct.h
contig_selector.o: gap-if.h
contig_selector.o: gap_globals.h
contig_selector.o: hash.h
contig_selector.o: io-reg.h
contig_selector.o: io_handle.h
contig_selector.o: io_utils.h
contig_selector.o: list.h
contig_selector.o: newgap_cmds.h
contig_selector.o: primlib.h
contig_selector.o: qual.h
contig_selector.o: ruler_display.h
contig_selector.o: tagDefs.h
contig_selector.o: tagUtils.h
contig_selector.o: tagdb.h
contig_selector.o: template.h
contig_selector.o: template_display.h
contig_selector.o: tkEdNames.h
contig_selector.o: tkEditor.h
contig_selector.o: tman_display.h
contig_selector.o: undo.h
copy_db.o: $(SRCROOT)/Misc/array.h
copy_db.o: $(SRCROOT)/Misc/bitmap.h
copy_db.o: $(SRCROOT)/Misc/misc.h
copy_db.o: $(SRCROOT)/Misc/os.h
copy_db.o: $(SRCROOT)/Misc/xalloc.h
copy_db.o: $(SRCROOT)/Misc/xerror.h
copy_db.o: $(SRCROOT)/g/freetree.h
copy_db.o: $(SRCROOT)/g/g-defs.h
copy_db.o: $(SRCROOT)/g/g-error.h
copy_db.o: $(SRCROOT)/g/g-filedefs.h
copy_db.o: $(SRCROOT)/g/g-os.h
copy_db.o: $(SRCROOT)/g/g-struct.h
copy_db.o: $(PWD)/staden_config.h
copy_db.o: $(SRCROOT)/tk_utils/text_output.h
copy_db.o: IO.h
copy_db.o: copy_db.h
copy_db.o: gap-dbstruct.h
copy_db.o: gap-if.h
copy_db.o: io_handle.h io_utils.h
copy_db_main.o: $(SRCROOT)/Misc/array.h
copy_db_main.o: $(SRCROOT)/Misc/bitmap.h
copy_db_main.o: $(SRCROOT)/Misc/misc.h
copy_db_main.o: $(SRCROOT)/Misc/os.h
copy_db_main.o: $(SRCROOT)/Misc/xalloc.h
copy_db_main.o: $(SRCROOT)/Misc/xerror.h
copy_db_main.o: $(SRCROOT)/g/freetree.h
copy_db_main.o: $(SRCROOT)/g/g-defs.h
copy_db_main.o: $(SRCROOT)/g/g-error.h
copy_db_main.o: $(SRCROOT)/g/g-filedefs.h
copy_db_main.o: $(SRCROOT)/g/g-os.h
copy_db_main.o: $(SRCROOT)/g/g-struct.h
copy_db_main.o: $(PWD)/staden_config.h
copy_db_main.o: copy_db.h IO.h
copy_db_main.o: gap-create.h
copy_db_main.o: gap-dbstruct.h
copy_db_main.o: gap-if.h
copy_db_main.o: io_handle.h
copy_db_main.o: io_utils.h
cs-object.o: $(SRCROOT)/Misc/array.h
cs-object.o: $(SRCROOT)/Misc/bitmap.h
cs-object.o: $(SRCROOT)/Misc/misc.h
cs-object.o: $(SRCROOT)/Misc/os.h
cs-object.o: $(SRCROOT)/Misc/xalloc.h
cs-object.o: $(SRCROOT)/Misc/xerror.h
cs-object.o: $(SRCROOT)/g/freetree.h
cs-object.o: $(SRCROOT)/g/g-defs.h
cs-object.o: $(SRCROOT)/g/g-error.h
cs-object.o: $(SRCROOT)/g/g-filedefs.h
cs-object.o: $(SRCROOT)/g/g-os.h
cs-object.o: $(SRCROOT)/g/g-struct.h
cs-object.o: $(PWD)/staden_config.h
cs-object.o: $(SRCROOT)/tk_utils/canvas_box.h
cs-object.o: $(SRCROOT)/tk_utils/tcl_utils.h
cs-object.o: $(SRCROOT)/tk_utils/text_output.h
cs-object.o: IO.h
cs-object.o: contig_selector.h
cs-object.o: cs-object.h
cs-object.o: gap-dbstruct.h
cs-object.o: gap-if.h
cs-object.o: io-reg.h
cs-object.o: io_handle.h
cs-object.o: io_utils.h
cs-object.o: list.h gap_globals.h
cs-object.o: newgap_cmds.h
cs-object.o: tk-io-reg.h
dbcheck.o: $(SRCROOT)/Misc/array.h
dbcheck.o: $(SRCROOT)/Misc/bitmap.h
dbcheck.o: $(SRCROOT)/Misc/misc.h
dbcheck.o: $(SRCROOT)/Misc/os.h
dbcheck.o: $(SRCROOT)/Misc/xalloc.h
dbcheck.o: $(SRCROOT)/Misc/xerror.h
dbcheck.o: $(SRCROOT)/g/freetree.h
dbcheck.o: $(SRCROOT)/g/g-defs.h
dbcheck.o: $(SRCROOT)/g/g-error.h
dbcheck.o: $(SRCROOT)/g/g-filedefs.h
dbcheck.o: $(SRCROOT)/g/g-os.h
dbcheck.o: $(SRCROOT)/g/g-struct.h
dbcheck.o: $(PWD)/staden_config.h
dbcheck.o: $(SRCROOT)/tk_utils/tcl_utils.h
dbcheck.o: $(SRCROOT)/tk_utils/text_output.h
dbcheck.o: IO.h
dbcheck.o: fortran.h fort.h
dbcheck.o: gap-dbstruct.h
dbcheck.o: gap-error.h dbcheck.h
dbcheck.o: gap-if.h
dbcheck.o: gap_globals.h
dbcheck.o: io_handle.h
dbcheck.o: io_utils.h
dis_readings.o: $(SRCROOT)/Misc/array.h
dis_readings.o: $(SRCROOT)/Misc/bitmap.h
dis_readings.o: $(SRCROOT)/Misc/dstring.h
dis_readings.o: $(SRCROOT)/Misc/misc.h
dis_readings.o: $(SRCROOT)/Misc/os.h
dis_readings.o: $(SRCROOT)/Misc/xalloc.h
dis_readings.o: $(SRCROOT)/Misc/xerror.h
dis_readings.o: $(SRCROOT)/g/freetree.h
dis_readings.o: $(SRCROOT)/g/g-defs.h
dis_readings.o: $(SRCROOT)/g/g-error.h
dis_readings.o: $(SRCROOT)/g/g-filedefs.h
dis_readings.o: $(SRCROOT)/g/g-os.h
dis_readings.o: $(SRCROOT)/g/g-struct.h
dis_readings.o: $(SRCROOT)/primer3/src/dpal.h
dis_readings.o: $(SRCROOT)/primer3/src/primer3.h
dis_readings.o: $(PWD)/staden_config.h
dis_readings.o: $(SRCROOT)/tk_utils/cli_arg.h
dis_readings.o: $(SRCROOT)/tk_utils/postscript.h
dis_readings.o: $(SRCROOT)/tk_utils/sheet.h
dis_readings.o: $(SRCROOT)/tk_utils/tkSheet.h
dis_readings.o: $(SRCROOT)/tk_utils/tkSheet_common.h
dis_readings.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
dis_readings.o: $(SRCROOT)/tk_utils/tkTrace.h
dis_readings.o: IO.h
dis_readings.o: IO2.h fort.h
dis_readings.o: dbcheck.h
dis_readings.o: dis_readings.h
dis_readings.o: edStructs.h
dis_readings.o: fortran.h notes.h
dis_readings.o: gap-dbstruct.h
dis_readings.o: gap-if.h
dis_readings.o: io-reg.h
dis_readings.o: io_handle.h
dis_readings.o: io_utils.h
dis_readings.o: list_proc.h
dis_readings.o: primlib.h
dis_readings.o: qual.h
dis_readings.o: tagDefs.h
dis_readings.o: tagUtils.h
dis_readings.o: template.h list.h
dis_readings.o: tkEdNames.h
dis_readings.o: tkEditor.h
dis_readings.o: tman_display.h
dis_readings.o: undo.h edUtils.h
do_fij.o: $(SRCROOT)/Misc/array.h
do_fij.o: $(SRCROOT)/Misc/bitmap.h
do_fij.o: $(SRCROOT)/Misc/misc.h
do_fij.o: $(SRCROOT)/Misc/os.h
do_fij.o: $(SRCROOT)/Misc/xalloc.h
do_fij.o: $(SRCROOT)/Misc/xerror.h
do_fij.o: $(SRCROOT)/g/freetree.h
do_fij.o: $(SRCROOT)/g/g-defs.h
do_fij.o: $(SRCROOT)/g/g-error.h
do_fij.o: $(SRCROOT)/g/g-filedefs.h
do_fij.o: $(SRCROOT)/g/g-os.h
do_fij.o: $(SRCROOT)/g/g-struct.h
do_fij.o: $(SRCROOT)/seq_utils/align_lib.h
do_fij.o: $(SRCROOT)/seq_utils/dna_utils.h
do_fij.o: $(PWD)/staden_config.h
do_fij.o: fij.h io_utils.h IO.h
do_fij.o: gap-dbstruct.h
do_fij.o: gap-if.h
do_fij.o: hash_lib.h
do_fij.o: io_handle.h consen.h
dstrand.o: $(SRCROOT)/Misc/FtoC.h
dstrand.o: $(SRCROOT)/Misc/array.h
dstrand.o: $(SRCROOT)/Misc/bitmap.h
dstrand.o: $(SRCROOT)/Misc/dstring.h
dstrand.o: $(SRCROOT)/Misc/misc.h
dstrand.o: $(SRCROOT)/Misc/os.h
dstrand.o: $(SRCROOT)/Misc/xalloc.h
dstrand.o: $(SRCROOT)/Misc/xerror.h
dstrand.o: $(SRCROOT)/g/freetree.h
dstrand.o: $(SRCROOT)/g/g-defs.h
dstrand.o: $(SRCROOT)/g/g-error.h
dstrand.o: $(SRCROOT)/g/g-filedefs.h
dstrand.o: $(SRCROOT)/g/g-os.h
dstrand.o: $(SRCROOT)/g/g-struct.h
dstrand.o: $(SRCROOT)/primer3/src/dpal.h
dstrand.o: $(SRCROOT)/primer3/src/primer3.h
dstrand.o: $(SRCROOT)/seq_utils/align.h
dstrand.o: $(SRCROOT)/seq_utils/align_lib_old.h
dstrand.o: $(PWD)/staden_config.h
dstrand.o: $(SRCROOT)/tk_utils/cli_arg.h
dstrand.o: $(SRCROOT)/tk_utils/postscript.h
dstrand.o: $(SRCROOT)/tk_utils/sheet.h
dstrand.o: $(SRCROOT)/tk_utils/tcl_utils.h
dstrand.o: $(SRCROOT)/tk_utils/text_output.h
dstrand.o: $(SRCROOT)/tk_utils/tkSheet.h
dstrand.o: $(SRCROOT)/tk_utils/tkSheet_common.h
dstrand.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
dstrand.o: $(SRCROOT)/tk_utils/tkTrace.h
dstrand.o: IO.h
dstrand.o: IO2.h complement.h
dstrand.o: dstrand.h
dstrand.o: edUtils.h edStructs.h
dstrand.o: fort.h tkEditor.h
dstrand.o: gap-dbstruct.h
dstrand.o: gap-if.h
dstrand.o: gap_globals.h
dstrand.o: io_handle.h qual.h
dstrand.o: io_utils.h
dstrand.o: primlib.h
dstrand.o: tagDefs.h fortran.h
dstrand.o: template.h list.h
dstrand.o: tkEdNames.h io-reg.h
dstrand.o: tman_display.h
dstrand.o: undo.h tagUtils.h
edCommands.o: $(SRCROOT)/Misc/array.h
edCommands.o: $(SRCROOT)/Misc/bitmap.h
edCommands.o: $(SRCROOT)/Misc/dstring.h
edCommands.o: $(SRCROOT)/Misc/misc.h
edCommands.o: $(SRCROOT)/Misc/os.h
edCommands.o: $(SRCROOT)/Misc/xalloc.h
edCommands.o: $(SRCROOT)/Misc/xerror.h
edCommands.o: $(SRCROOT)/g/freetree.h
edCommands.o: $(SRCROOT)/g/g-defs.h
edCommands.o: $(SRCROOT)/g/g-error.h
edCommands.o: $(SRCROOT)/g/g-filedefs.h
edCommands.o: $(SRCROOT)/g/g-os.h
edCommands.o: $(SRCROOT)/g/g-struct.h
edCommands.o: $(SRCROOT)/primer3/src/dpal.h
edCommands.o: $(SRCROOT)/primer3/src/primer3.h
edCommands.o: $(SRCROOT)/seq_utils/align.h
edCommands.o: $(SRCROOT)/seq_utils/align_lib_old.h
edCommands.o: $(PWD)/staden_config.h
edCommands.o: $(SRCROOT)/tk_utils/cli_arg.h
edCommands.o: $(SRCROOT)/tk_utils/postscript.h
edCommands.o: $(SRCROOT)/tk_utils/sheet.h
edCommands.o: $(SRCROOT)/tk_utils/tcl_utils.h
edCommands.o: $(SRCROOT)/tk_utils/tkSheet.h
edCommands.o: $(SRCROOT)/tk_utils/tkSheet_common.h
edCommands.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
edCommands.o: $(SRCROOT)/tk_utils/tkTrace.h
edCommands.o: IO.h
edCommands.o: edStructs.h
edCommands.o: edUtils.h qual.h
edCommands.o: fort.h io_utils.h
edCommands.o: gap-dbstruct.h
edCommands.o: gap-if.h
edCommands.o: gap_globals.h
edCommands.o: io-reg.h primlib.h
edCommands.o: io_handle.h
edCommands.o: select.h
edCommands.o: tagDefs.h fortran.h
edCommands.o: template.h list.h
edCommands.o: tkEdNames.h
edCommands.o: tkEditor.h
edCommands.o: tman_display.h
edCommands.o: tman_interface.h
edCommands.o: undo.h tagUtils.h
edExtend.o: $(SRCROOT)/Misc/array.h
edExtend.o: $(SRCROOT)/Misc/bitmap.h
edExtend.o: $(SRCROOT)/Misc/dstring.h
edExtend.o: $(SRCROOT)/Misc/os.h
edExtend.o: $(SRCROOT)/Misc/xalloc.h
edExtend.o: $(SRCROOT)/Misc/xerror.h
edExtend.o: $(SRCROOT)/g/freetree.h
edExtend.o: $(SRCROOT)/g/g-defs.h
edExtend.o: $(SRCROOT)/g/g-error.h
edExtend.o: $(SRCROOT)/g/g-filedefs.h
edExtend.o: $(SRCROOT)/g/g-os.h
edExtend.o: $(SRCROOT)/g/g-struct.h
edExtend.o: $(SRCROOT)/primer3/src/dpal.h
edExtend.o: $(SRCROOT)/primer3/src/primer3.h
edExtend.o: $(PWD)/staden_config.h
edExtend.o: $(SRCROOT)/tk_utils/cli_arg.h
edExtend.o: $(SRCROOT)/tk_utils/postscript.h
edExtend.o: $(SRCROOT)/tk_utils/sheet.h
edExtend.o: $(SRCROOT)/tk_utils/tkSheet.h
edExtend.o: $(SRCROOT)/tk_utils/tkSheet_common.h
edExtend.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
edExtend.o: $(SRCROOT)/tk_utils/tkTrace.h
edExtend.o: IO.h
edExtend.o: contigEditor.h
edExtend.o: edStructs.h
edExtend.o: edUtils.h qual.h
edExtend.o: fort.h io_utils.h
edExtend.o: gap-dbstruct.h
edExtend.o: gap-if.h
edExtend.o: io_handle.h
edExtend.o: primlib.h
edExtend.o: select.h extend.h
edExtend.o: tagDefs.h fortran.h
edExtend.o: template.h list.h
edExtend.o: tkEdNames.h io-reg.h
edExtend.o: tkEditor.h
edExtend.o: tman_display.h
edExtend.o: undo.h tagUtils.h
edInterface.o: $(SRCROOT)/Misc/array.h
edInterface.o: $(SRCROOT)/Misc/bitmap.h
edInterface.o: $(SRCROOT)/Misc/dstring.h
edInterface.o: $(SRCROOT)/Misc/misc.h
edInterface.o: $(SRCROOT)/Misc/os.h
edInterface.o: $(SRCROOT)/Misc/xalloc.h
edInterface.o: $(SRCROOT)/Misc/xerror.h
edInterface.o: $(SRCROOT)/g/freetree.h
edInterface.o: $(SRCROOT)/g/g-defs.h
edInterface.o: $(SRCROOT)/g/g-error.h
edInterface.o: $(SRCROOT)/g/g-filedefs.h
edInterface.o: $(SRCROOT)/g/g-os.h
edInterface.o: $(SRCROOT)/g/g-struct.h
edInterface.o: $(SRCROOT)/primer3/src/dpal.h
edInterface.o: $(SRCROOT)/primer3/src/primer3.h
edInterface.o: $(PWD)/staden_config.h
edInterface.o: $(SRCROOT)/tk_utils/canvas_box.h
edInterface.o: $(SRCROOT)/tk_utils/cli_arg.h
edInterface.o: $(SRCROOT)/tk_utils/postscript.h
edInterface.o: $(SRCROOT)/tk_utils/sheet.h
edInterface.o: $(SRCROOT)/tk_utils/tcl_utils.h
edInterface.o: $(SRCROOT)/tk_utils/tkSheet.h
edInterface.o: $(SRCROOT)/tk_utils/tkSheet_common.h
edInterface.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
edInterface.o: $(SRCROOT)/tk_utils/tkTrace.h
edInterface.o: IO.h
edInterface.o: active_tags.h
edInterface.o: contigEditor.h
edInterface.o: edStructs.h
edInterface.o: edUtils.h qual.h
edInterface.o: extend.h
edInterface.o: fort.h io_utils.h
edInterface.o: fortran.h
edInterface.o: gap-dbstruct.h
edInterface.o: gap-if.h
edInterface.o: gap_globals.h
edInterface.o: hash.h
edInterface.o: io-reg.h primlib.h
edInterface.o: io_handle.h
edInterface.o: qualIO.h notes.h
edInterface.o: ruler_display.h
edInterface.o: tagDefs.h
edInterface.o: tagdb.h
edInterface.o: template.h list.h
edInterface.o: template_display.h
edInterface.o: tkEdNames.h
edInterface.o: tkEditor.h
edInterface.o: tman_display.h
edInterface.o: tman_interface.h
edInterface.o: undo.h tagUtils.h
edMutations.o: $(SRCROOT)/Misc/array.h
edMutations.o: $(SRCROOT)/Misc/bitmap.h
edMutations.o: $(SRCROOT)/Misc/dstring.h
edMutations.o: $(SRCROOT)/Misc/misc.h
edMutations.o: $(SRCROOT)/Misc/os.h
edMutations.o: $(SRCROOT)/Misc/xalloc.h
edMutations.o: $(SRCROOT)/Misc/xerror.h
edMutations.o: $(SRCROOT)/g/freetree.h
edMutations.o: $(SRCROOT)/g/g-defs.h
edMutations.o: $(SRCROOT)/g/g-error.h
edMutations.o: $(SRCROOT)/g/g-filedefs.h
edMutations.o: $(SRCROOT)/g/g-os.h
edMutations.o: $(SRCROOT)/g/g-struct.h
edMutations.o: $(SRCROOT)/primer3/src/dpal.h
edMutations.o: $(SRCROOT)/primer3/src/primer3.h
edMutations.o: $(SRCROOT)/seq_utils/dna_utils.h
edMutations.o: $(SRCROOT)/seq_utils/genetic_code.h
edMutations.o: $(PWD)/staden_config.h
edMutations.o: $(SRCROOT)/tk_utils/cli_arg.h
edMutations.o: $(SRCROOT)/tk_utils/postscript.h
edMutations.o: $(SRCROOT)/tk_utils/sheet.h
edMutations.o: $(SRCROOT)/tk_utils/tkSheet.h
edMutations.o: $(SRCROOT)/tk_utils/tkSheet_common.h
edMutations.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
edMutations.o: $(SRCROOT)/tk_utils/tkTrace.h
edMutations.o: IO.h
edMutations.o: edStructs.h
edMutations.o: edUtils.h qual.h
edMutations.o: fort.h io_utils.h
edMutations.o: fortran.h
edMutations.o: gap-dbstruct.h
edMutations.o: gap-if.h
edMutations.o: io-reg.h primlib.h
edMutations.o: io_handle.h
edMutations.o: parse_ft.h
edMutations.o: tagDefs.h
edMutations.o: template.h list.h
edMutations.o: tkEdNames.h
edMutations.o: tkEditor.h
edMutations.o: tman_display.h
edMutations.o: tman_interface.h
edMutations.o: undo.h tagUtils.h
edUtils2.o: $(SRCROOT)/Misc/array.h
edUtils2.o: $(SRCROOT)/Misc/bitmap.h
edUtils2.o: $(SRCROOT)/Misc/dstring.h
edUtils2.o: $(SRCROOT)/Misc/locks.h
edUtils2.o: $(SRCROOT)/Misc/misc.h
edUtils2.o: $(SRCROOT)/Misc/os.h
edUtils2.o: $(SRCROOT)/Misc/xalloc.h
edUtils2.o: $(SRCROOT)/Misc/xerror.h
edUtils2.o: $(SRCROOT)/g/freetree.h
edUtils2.o: $(SRCROOT)/g/g-defs.h
edUtils2.o: $(SRCROOT)/g/g-error.h
edUtils2.o: $(SRCROOT)/g/g-filedefs.h
edUtils2.o: $(SRCROOT)/g/g-os.h
edUtils2.o: $(SRCROOT)/g/g-struct.h
edUtils2.o: $(SRCROOT)/primer3/src/dpal.h
edUtils2.o: $(SRCROOT)/primer3/src/primer3.h
edUtils2.o: $(PWD)/staden_config.h
edUtils2.o: $(SRCROOT)/tk_utils/cli_arg.h
edUtils2.o: $(SRCROOT)/tk_utils/postscript.h
edUtils2.o: $(SRCROOT)/tk_utils/sheet.h
edUtils2.o: $(SRCROOT)/tk_utils/tcl_utils.h
edUtils2.o: $(SRCROOT)/tk_utils/tkSheet.h
edUtils2.o: $(SRCROOT)/tk_utils/tkSheet_common.h
edUtils2.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
edUtils2.o: $(SRCROOT)/tk_utils/tkTrace.h
edUtils2.o: IO.h
edUtils2.o: IO2.h active_tags.h
edUtils2.o: contigEditor.h
edUtils2.o: edUtils.h edStructs.h
edUtils2.o: fort.h io_utils.h
edUtils2.o: gap-dbstruct.h
edUtils2.o: gap-if.h
edUtils2.o: io_handle.h
edUtils2.o: notes.h
edUtils2.o: primlib.h
edUtils2.o: qual.h
edUtils2.o: select.h extend.h
edUtils2.o: tagDefs.h fortran.h
edUtils2.o: tagdb.h gap_globals.h
edUtils2.o: template.h list.h
edUtils2.o: tkEdNames.h io-reg.h
edUtils2.o: tkEditor.h
edUtils2.o: tman_display.h
edUtils2.o: tman_interface.h
edUtils2.o: undo.h tagUtils.h
extract.o: $(SRCROOT)/Misc/FtoC.h
extract.o: $(SRCROOT)/Misc/array.h
extract.o: $(SRCROOT)/Misc/bitmap.h
extract.o: $(SRCROOT)/Misc/dstring.h
extract.o: $(SRCROOT)/Misc/misc.h
extract.o: $(SRCROOT)/Misc/os.h
extract.o: $(SRCROOT)/Misc/xalloc.h
extract.o: $(SRCROOT)/Misc/xerror.h
extract.o: $(SRCROOT)/g/freetree.h
extract.o: $(SRCROOT)/g/g-defs.h
extract.o: $(SRCROOT)/g/g-error.h
extract.o: $(SRCROOT)/g/g-filedefs.h
extract.o: $(SRCROOT)/g/g-os.h
extract.o: $(SRCROOT)/g/g-struct.h
extract.o: $(SRCROOT)/primer3/src/dpal.h
extract.o: $(SRCROOT)/primer3/src/primer3.h
extract.o: $(PWD)/staden_config.h
extract.o: $(SRCROOT)/tk_utils/cli_arg.h
extract.o: $(SRCROOT)/tk_utils/postscript.h
extract.o: $(SRCROOT)/tk_utils/sheet.h
extract.o: $(SRCROOT)/tk_utils/text_output.h
extract.o: $(SRCROOT)/tk_utils/tkSheet.h
extract.o: $(SRCROOT)/tk_utils/tkSheet_common.h
extract.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
extract.o: $(SRCROOT)/tk_utils/tkTrace.h
extract.o: IO.h
extract.o: edStructs.h
extract.o: edUtils.h qual.h
extract.o: extract.h
extract.o: fort.h
extract.o: gap-dbstruct.h
extract.o: gap-if.h
extract.o: io_handle.h io_utils.h
extract.o: notes.h
extract.o: primlib.h
extract.o: tagDefs.h fortran.h
extract.o: template.h list.h
extract.o: tkEdNames.h io-reg.h
extract.o: tkEditor.h
extract.o: tman_display.h
extract.o: undo.h tagUtils.h
f2c.o: $(SRCROOT)/Misc/misc.h
f2c.o: $(SRCROOT)/Misc/os.h
f2c.o: $(SRCROOT)/Misc/xalloc.h
f2c.o: $(PWD)/staden_config.h
f2c.o: f2c.h
fij.o: $(SRCROOT)/Misc/array.h
fij.o: $(SRCROOT)/Misc/bitmap.h
fij.o: $(SRCROOT)/Misc/dstring.h
fij.o: $(SRCROOT)/Misc/misc.h
fij.o: $(SRCROOT)/Misc/os.h
fij.o: $(SRCROOT)/Misc/xalloc.h
fij.o: $(SRCROOT)/Misc/xerror.h
fij.o: $(SRCROOT)/g/freetree.h
fij.o: $(SRCROOT)/g/g-defs.h
fij.o: $(SRCROOT)/g/g-error.h
fij.o: $(SRCROOT)/g/g-filedefs.h
fij.o: $(SRCROOT)/g/g-os.h
fij.o: $(SRCROOT)/g/g-struct.h
fij.o: $(SRCROOT)/primer3/src/dpal.h
fij.o: $(SRCROOT)/primer3/src/primer3.h
fij.o: $(SRCROOT)/seq_utils/align_lib.h
fij.o: $(PWD)/staden_config.h
fij.o: $(SRCROOT)/tk_utils/canvas_box.h
fij.o: $(SRCROOT)/tk_utils/cli_arg.h
fij.o: $(SRCROOT)/tk_utils/postscript.h
fij.o: $(SRCROOT)/tk_utils/sheet.h
fij.o: $(SRCROOT)/tk_utils/tclXkeylist.h
fij.o: $(SRCROOT)/tk_utils/tcl_utils.h
fij.o: $(SRCROOT)/tk_utils/text_output.h
fij.o: $(SRCROOT)/tk_utils/tkSheet.h
fij.o: $(SRCROOT)/tk_utils/tkSheet_common.h
fij.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
fij.o: $(SRCROOT)/tk_utils/tkTrace.h
fij.o: IO.h
fij.o: complement.h
fij.o: contigEditor.h edStructs.h
fij.o: cs-object.h
fij.o: edUtils.h qual.h
fij.o: fij.h consen.h hash_lib.h
fij.o: fortran.h
fij.o: gap-dbstruct.h
fij.o: gap-if.h
fij.o: gap_globals.h
fij.o: io-reg.h contig_selector.h
fij.o: io_handle.h fort.h
fij.o: io_utils.h
fij.o: newgap_cmds.h list.h
fij.o: tagUtils.h tagDefs.h
fij.o: template.h undo.h
fij.o: tkEdNames.h primlib.h
fij.o: tkEditor.h
fij.o: tman_display.h
find_fragments.o: $(SRCROOT)/Misc/array.h
find_fragments.o: $(SRCROOT)/Misc/bitmap.h
find_fragments.o: $(SRCROOT)/Misc/misc.h
find_fragments.o: $(SRCROOT)/Misc/os.h
find_fragments.o: $(SRCROOT)/Misc/xalloc.h
find_fragments.o: $(SRCROOT)/Misc/xerror.h
find_fragments.o: $(SRCROOT)/g/freetree.h
find_fragments.o: $(SRCROOT)/g/g-defs.h
find_fragments.o: $(SRCROOT)/g/g-error.h
find_fragments.o: $(SRCROOT)/g/g-filedefs.h
find_fragments.o: $(SRCROOT)/g/g-os.h
find_fragments.o: $(SRCROOT)/g/g-struct.h
find_fragments.o: $(PWD)/staden_config.h
find_fragments.o: IO.h
find_fragments.o: find_fragments.h
find_fragments.o: gap-dbstruct.h
find_fragments.o: gap-if.h
find_fragments.o: io_handle.h
find_fragments.o: io_utils.h
find_fragments.o: qual.h
find_oligo.o: $(SRCROOT)/Misc/array.h
find_oligo.o: $(SRCROOT)/Misc/bitmap.h
find_oligo.o: $(SRCROOT)/Misc/dstring.h
find_oligo.o: $(SRCROOT)/Misc/misc.h
find_oligo.o: $(SRCROOT)/Misc/os.h
find_oligo.o: $(SRCROOT)/Misc/xalloc.h
find_oligo.o: $(SRCROOT)/Misc/xerror.h
find_oligo.o: $(SRCROOT)/g/freetree.h
find_oligo.o: $(SRCROOT)/g/g-defs.h
find_oligo.o: $(SRCROOT)/g/g-error.h
find_oligo.o: $(SRCROOT)/g/g-filedefs.h
find_oligo.o: $(SRCROOT)/g/g-os.h
find_oligo.o: $(SRCROOT)/g/g-struct.h
find_oligo.o: $(SRCROOT)/primer3/src/dpal.h
find_oligo.o: $(SRCROOT)/primer3/src/primer3.h
find_oligo.o: $(SRCROOT)/seq_utils/dna_utils.h
find_oligo.o: $(SRCROOT)/seq_utils/search_utils.h
find_oligo.o: $(SRCROOT)/seq_utils/sequence_formats.h
find_oligo.o: $(PWD)/staden_config.h
find_oligo.o: $(SRCROOT)/tk_utils/canvas_box.h
find_oligo.o: $(SRCROOT)/tk_utils/cli_arg.h
find_oligo.o: $(SRCROOT)/tk_utils/postscript.h
find_oligo.o: $(SRCROOT)/tk_utils/sheet.h
find_oligo.o: $(SRCROOT)/tk_utils/tclXkeylist.h
find_oligo.o: $(SRCROOT)/tk_utils/tcl_utils.h
find_oligo.o: $(SRCROOT)/tk_utils/text_output.h
find_oligo.o: $(SRCROOT)/tk_utils/tkSheet.h
find_oligo.o: $(SRCROOT)/tk_utils/tkSheet_common.h
find_oligo.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
find_oligo.o: $(SRCROOT)/tk_utils/tkTrace.h
find_oligo.o: IO.h
find_oligo.o: complement.h
find_oligo.o: contigEditor.h
find_oligo.o: contig_selector.h
find_oligo.o: cs-object.h
find_oligo.o: edStructs.h
find_oligo.o: fort.h tkEditor.h
find_oligo.o: gap-dbstruct.h
find_oligo.o: gap-if.h
find_oligo.o: gap_globals.h
find_oligo.o: hash.h
find_oligo.o: io-reg.h primlib.h
find_oligo.o: io_handle.h tagdb.h
find_oligo.o: io_utils.h
find_oligo.o: qual.h
find_oligo.o: ruler_display.h
find_oligo.o: tagDefs.h fortran.h
find_oligo.o: tagUtils.h
find_oligo.o: template.h list.h
find_oligo.o: template_display.h
find_oligo.o: tkEdNames.h
find_oligo.o: tman_display.h
find_oligo.o: undo.h edUtils.h
find_repeats.o: $(SRCROOT)/Misc/array.h
find_repeats.o: $(SRCROOT)/Misc/bitmap.h
find_repeats.o: $(SRCROOT)/Misc/dstring.h
find_repeats.o: $(SRCROOT)/Misc/misc.h
find_repeats.o: $(SRCROOT)/Misc/os.h
find_repeats.o: $(SRCROOT)/Misc/xalloc.h
find_repeats.o: $(SRCROOT)/Misc/xerror.h
find_repeats.o: $(SRCROOT)/g/freetree.h
find_repeats.o: $(SRCROOT)/g/g-defs.h
find_repeats.o: $(SRCROOT)/g/g-error.h
find_repeats.o: $(SRCROOT)/g/g-filedefs.h
find_repeats.o: $(SRCROOT)/g/g-os.h
find_repeats.o: $(SRCROOT)/g/g-struct.h
find_repeats.o: $(SRCROOT)/primer3/src/dpal.h
find_repeats.o: $(SRCROOT)/primer3/src/primer3.h
find_repeats.o: $(PWD)/staden_config.h
find_repeats.o: $(SRCROOT)/tk_utils/canvas_box.h
find_repeats.o: $(SRCROOT)/tk_utils/cli_arg.h
find_repeats.o: $(SRCROOT)/tk_utils/postscript.h
find_repeats.o: $(SRCROOT)/tk_utils/sheet.h
find_repeats.o: $(SRCROOT)/tk_utils/tclXkeylist.h
find_repeats.o: $(SRCROOT)/tk_utils/tcl_utils.h
find_repeats.o: $(SRCROOT)/tk_utils/text_output.h
find_repeats.o: $(SRCROOT)/tk_utils/tkSheet.h
find_repeats.o: $(SRCROOT)/tk_utils/tkSheet_common.h
find_repeats.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
find_repeats.o: $(SRCROOT)/tk_utils/tkTrace.h
find_repeats.o: IO.h
find_repeats.o: complement.h
find_repeats.o: consen.h
find_repeats.o: contigEditor.h
find_repeats.o: contig_selector.h
find_repeats.o: cs-object.h
find_repeats.o: edStructs.h
find_repeats.o: find_repeats.h
find_repeats.o: fort.h
find_repeats.o: fortran.h
find_repeats.o: gap-dbstruct.h
find_repeats.o: gap-if.h
find_repeats.o: gap_globals.h
find_repeats.o: gap_hash.h
find_repeats.o: io-reg.h
find_repeats.o: io_handle.h
find_repeats.o: io_utils.h
find_repeats.o: primlib.h
find_repeats.o: qual.h tagUtils.h
find_repeats.o: tagDefs.h
find_repeats.o: template.h list.h
find_repeats.o: tkEdNames.h
find_repeats.o: tkEditor.h
find_repeats.o: tman_display.h
find_repeats.o: undo.h edUtils.h
gap-create.o: $(SRCROOT)/Misc/FtoC.h
gap-create.o: $(SRCROOT)/Misc/array.h
gap-create.o: $(SRCROOT)/Misc/bitmap.h
gap-create.o: $(SRCROOT)/Misc/os.h
gap-create.o: $(SRCROOT)/Misc/xerror.h
gap-create.o: $(SRCROOT)/g/freetree.h
gap-create.o: $(SRCROOT)/g/g-defs.h
gap-create.o: $(SRCROOT)/g/g-error.h
gap-create.o: $(SRCROOT)/g/g-filedefs.h
gap-create.o: $(SRCROOT)/g/g-io.h
gap-create.o: $(SRCROOT)/g/g-misc.h
gap-create.o: $(SRCROOT)/g/g-os.h
gap-create.o: $(SRCROOT)/g/g-struct.h
gap-create.o: $(PWD)/staden_config.h
gap-create.o: gap-create.h
gap-create.o: gap-dbstruct.h
gap-create.o: gap-defaults.h
gap-create.o: gap-error.h
gap-create.o: gap-if.h
gap-create.o: gap-init.h gap-io.h
gap-dbstruct.o: $(SRCROOT)/Misc/os.h
gap-dbstruct.o: $(SRCROOT)/Misc/xerror.h
gap-dbstruct.o: $(SRCROOT)/g/g-error.h
gap-dbstruct.o: $(SRCROOT)/g/g-os.h
gap-dbstruct.o: $(PWD)/staden_config.h
gap-dbstruct.o: gap-dbstruct.h
gap-error.o: $(SRCROOT)/Misc/array.h
gap-error.o: $(SRCROOT)/Misc/bitmap.h
gap-error.o: $(SRCROOT)/Misc/misc.h
gap-error.o: $(SRCROOT)/Misc/os.h
gap-error.o: $(SRCROOT)/Misc/xalloc.h
gap-error.o: $(SRCROOT)/Misc/xerror.h
gap-error.o: $(SRCROOT)/g/g-error.h
gap-error.o: $(PWD)/staden_config.h
gap-error.o: gap-error.h
gap-error.o: stack_dump.h
gap-if.o: $(SRCROOT)/Misc/array.h
gap-if.o: $(SRCROOT)/Misc/os.h
gap-if.o: $(SRCROOT)/Misc/xerror.h
gap-if.o: $(SRCROOT)/g/freetree.h
gap-if.o: $(SRCROOT)/g/g-error.h
gap-if.o: $(SRCROOT)/g/g-filedefs.h
gap-if.o: $(SRCROOT)/g/g-os.h
gap-if.o: $(SRCROOT)/g/g-struct.h
gap-if.o: $(PWD)/staden_config.h
gap-if.o: gap-if.h
gap-if.o: gap-local.h
gap-if.o: gap-remote.h gap-init.h
gap-init.o: $(SRCROOT)/Misc/array.h
gap-init.o: $(SRCROOT)/Misc/os.h
gap-init.o: $(SRCROOT)/Misc/xerror.h
gap-init.o: $(SRCROOT)/g/freetree.h
gap-init.o: $(SRCROOT)/g/g-error.h
gap-init.o: $(SRCROOT)/g/g-filedefs.h
gap-init.o: $(SRCROOT)/g/g-os.h
gap-init.o: $(SRCROOT)/g/g-struct.h
gap-init.o: $(PWD)/staden_config.h
gap-init.o: gap-init.h gap-if.h
gap-init.o: gap-io.h
gap-io.o: $(SRCROOT)/Misc/array.h
gap-io.o: $(SRCROOT)/Misc/os.h
gap-io.o: $(SRCROOT)/Misc/xalloc.h
gap-io.o: $(SRCROOT)/Misc/xerror.h
gap-io.o: $(SRCROOT)/g/freetree.h
gap-io.o: $(SRCROOT)/g/g-error.h
gap-io.o: $(SRCROOT)/g/g-filedefs.h
gap-io.o: $(SRCROOT)/g/g-os.h
gap-io.o: $(SRCROOT)/g/g-struct.h
gap-io.o: $(PWD)/staden_config.h
gap-io.o: gap-error.h gap-if.h
gap-local.o: $(SRCROOT)/Misc/array.h
gap-local.o: $(SRCROOT)/Misc/os.h
gap-local.o: $(SRCROOT)/Misc/xalloc.h
gap-local.o: $(SRCROOT)/Misc/xerror.h
gap-local.o: $(SRCROOT)/g/freetree.h
gap-local.o: $(SRCROOT)/g/g-connect.h
gap-local.o: $(SRCROOT)/g/g-db.h
gap-local.o: $(SRCROOT)/g/g-error.h
gap-local.o: $(SRCROOT)/g/g-filedefs.h
gap-local.o: $(SRCROOT)/g/g-os.h
gap-local.o: $(SRCROOT)/g/g-request.h
gap-local.o: $(SRCROOT)/g/g-struct.h
gap-local.o: $(PWD)/staden_config.h
gap-local.o: gap-dbstruct.h
gap-local.o: gap-local.h gap-if.h
gap-remote.o: $(SRCROOT)/Misc/array.h
gap-remote.o: $(SRCROOT)/Misc/os.h
gap-remote.o: $(SRCROOT)/Misc/xerror.h
gap-remote.o: $(SRCROOT)/g/freetree.h
gap-remote.o: $(SRCROOT)/g/g-connect.h
gap-remote.o: $(SRCROOT)/g/g-db.h
gap-remote.o: $(SRCROOT)/g/g-error.h
gap-remote.o: $(SRCROOT)/g/g-filedefs.h
gap-remote.o: $(SRCROOT)/g/g-os.h
gap-remote.o: $(SRCROOT)/g/g-request.h
gap-remote.o: $(SRCROOT)/g/g-struct.h
gap-remote.o: $(PWD)/staden_config.h
gap-remote.o: gap-dbstruct.h
gap-remote.o: gap-if.h
gap-remote.o: gap-remote.h
gap-tcl.o: $(SRCROOT)/Misc/array.h
gap-tcl.o: $(SRCROOT)/Misc/bitmap.h
gap-tcl.o: $(SRCROOT)/Misc/dstring.h
gap-tcl.o: $(SRCROOT)/Misc/os.h
gap-tcl.o: $(SRCROOT)/Misc/xalloc.h
gap-tcl.o: $(SRCROOT)/Misc/xerror.h
gap-tcl.o: $(SRCROOT)/g/freetree.h
gap-tcl.o: $(SRCROOT)/g/g-defs.h
gap-tcl.o: $(SRCROOT)/g/g-error.h
gap-tcl.o: $(SRCROOT)/g/g-filedefs.h
gap-tcl.o: $(SRCROOT)/g/g-os.h
gap-tcl.o: $(SRCROOT)/g/g-struct.h
gap-tcl.o: $(SRCROOT)/primer3/src/dpal.h
gap-tcl.o: $(SRCROOT)/primer3/src/primer3.h
gap-tcl.o: $(PWD)/staden_config.h
gap-tcl.o: $(SRCROOT)/tk_utils/cli_arg.h
gap-tcl.o: $(SRCROOT)/tk_utils/postscript.h
gap-tcl.o: $(SRCROOT)/tk_utils/sheet.h
gap-tcl.o: $(SRCROOT)/tk_utils/tclXkeylist.h
gap-tcl.o: $(SRCROOT)/tk_utils/tcl_utils.h
gap-tcl.o: $(SRCROOT)/tk_utils/tkSheet.h
gap-tcl.o: $(SRCROOT)/tk_utils/tkSheet_common.h
gap-tcl.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
gap-tcl.o: $(SRCROOT)/tk_utils/tkTrace.h
gap-tcl.o: IO.h
gap-tcl.o: clones.h seqInfo.h
gap-tcl.o: edStructs.h
gap-tcl.o: fort.h tkEditor.h
gap-tcl.o: gap-dbstruct.h
gap-tcl.o: gap-error.h tagUtils.h
gap-tcl.o: gap-if.h
gap-tcl.o: gap-tcl.h
gap-tcl.o: gap_globals.h
gap-tcl.o: io_handle.h io_utils.h
gap-tcl.o: primlib.h
gap-tcl.o: qual.h tman_display.h
gap-tcl.o: tagDefs.h fortran.h
gap-tcl.o: template.h list.h
gap-tcl.o: tkEdNames.h io-reg.h
gap-tcl.o: undo.h edUtils.h
gap_array.o: $(SRCROOT)/tk_utils/tcl_utils.h
gap_array.o: active_tags.h
gap_array.o: gap_array.h
gap_array.o: tagdb.h
gap_canvas_box.o: $(SRCROOT)/Misc/array.h
gap_canvas_box.o: $(SRCROOT)/Misc/bitmap.h
gap_canvas_box.o: $(SRCROOT)/Misc/os.h
gap_canvas_box.o: $(SRCROOT)/Misc/xerror.h
gap_canvas_box.o: $(SRCROOT)/g/freetree.h
gap_canvas_box.o: $(SRCROOT)/g/g-defs.h
gap_canvas_box.o: $(SRCROOT)/g/g-error.h
gap_canvas_box.o: $(SRCROOT)/g/g-filedefs.h
gap_canvas_box.o: $(SRCROOT)/g/g-os.h
gap_canvas_box.o: $(SRCROOT)/g/g-struct.h
gap_canvas_box.o: $(PWD)/staden_config.h
gap_canvas_box.o: $(SRCROOT)/tk_utils/canvas_box.h
gap_canvas_box.o: $(SRCROOT)/tk_utils/tcl_utils.h
gap_canvas_box.o: IO.h
gap_canvas_box.o: gap-dbstruct.h
gap_canvas_box.o: gap-if.h
gap_canvas_box.o: gap_globals.h
gap_canvas_box.o: io-reg.h
gap_canvas_box.o: io_handle.h
gap_canvas_box.o: io_utils.h
gap_cli_arg.o: $(SRCROOT)/Misc/array.h
gap_cli_arg.o: $(SRCROOT)/Misc/bitmap.h
gap_cli_arg.o: $(SRCROOT)/Misc/misc.h
gap_cli_arg.o: $(SRCROOT)/Misc/os.h
gap_cli_arg.o: $(SRCROOT)/Misc/xalloc.h
gap_cli_arg.o: $(SRCROOT)/Misc/xerror.h
gap_cli_arg.o: $(SRCROOT)/g/freetree.h
gap_cli_arg.o: $(SRCROOT)/g/g-defs.h
gap_cli_arg.o: $(SRCROOT)/g/g-error.h
gap_cli_arg.o: $(SRCROOT)/g/g-filedefs.h
gap_cli_arg.o: $(SRCROOT)/g/g-os.h
gap_cli_arg.o: $(SRCROOT)/g/g-struct.h
gap_cli_arg.o: $(PWD)/staden_config.h
gap_cli_arg.o: $(SRCROOT)/tk_utils/cli_arg.h
gap_cli_arg.o: IO.h
gap_cli_arg.o: gap-dbstruct.h
gap_cli_arg.o: gap-if.h
gap_cli_arg.o: gap_cli_arg.h
gap_cli_arg.o: io_handle.h
gap_cli_arg.o: io_utils.h
gap_globals.o: $(SRCROOT)/Misc/array.h
gap_globals.o: $(SRCROOT)/Misc/bitmap.h
gap_globals.o: $(SRCROOT)/Misc/misc.h
gap_globals.o: $(SRCROOT)/Misc/os.h
gap_globals.o: $(SRCROOT)/Misc/xalloc.h
gap_globals.o: $(SRCROOT)/Misc/xerror.h
gap_globals.o: $(SRCROOT)/g/freetree.h
gap_globals.o: $(SRCROOT)/g/g-defs.h
gap_globals.o: $(SRCROOT)/g/g-error.h
gap_globals.o: $(SRCROOT)/g/g-filedefs.h
gap_globals.o: $(SRCROOT)/g/g-os.h
gap_globals.o: $(SRCROOT)/g/g-struct.h
gap_globals.o: $(SRCROOT)/seq_utils/align.h
gap_globals.o: $(SRCROOT)/seq_utils/align_lib.h
gap_globals.o: $(SRCROOT)/seq_utils/align_lib_old.h
gap_globals.o: $(SRCROOT)/seq_utils/dna_utils.h
gap_globals.o: $(SRCROOT)/seq_utils/genetic_code.h
gap_globals.o: $(SRCROOT)/seq_utils/read_matrix.h
gap_globals.o: $(PWD)/staden_config.h
gap_globals.o: $(SRCROOT)/tk_utils/tclXkeylist.h
gap_globals.o: $(SRCROOT)/tk_utils/tcl_utils.h
gap_globals.o: consen.h IO.h
gap_globals.o: fort.h
gap_globals.o: gap-dbstruct.h
gap_globals.o: gap-if.h
gap_globals.o: gap_globals.h
gap_globals.o: hash_lib.h
gap_globals.o: io_handle.h
gap_globals.o: io_utils.h
gap_globals.o: notedb.h
gap_hash.o: $(SRCROOT)/Misc/array.h
gap_hash.o: $(SRCROOT)/Misc/bitmap.h
gap_hash.o: $(SRCROOT)/Misc/misc.h
gap_hash.o: $(SRCROOT)/Misc/os.h
gap_hash.o: $(SRCROOT)/Misc/xalloc.h
gap_hash.o: $(SRCROOT)/Misc/xerror.h
gap_hash.o: $(SRCROOT)/g/freetree.h
gap_hash.o: $(SRCROOT)/g/g-defs.h
gap_hash.o: $(SRCROOT)/g/g-error.h
gap_hash.o: $(SRCROOT)/g/g-filedefs.h
gap_hash.o: $(SRCROOT)/g/g-os.h
gap_hash.o: $(SRCROOT)/g/g-struct.h
gap_hash.o: $(SRCROOT)/seq_utils/align_lib.h
gap_hash.o: $(SRCROOT)/seq_utils/dna_utils.h
gap_hash.o: $(PWD)/staden_config.h
gap_hash.o: fort.h
gap_hash.o: gap-dbstruct.h
gap_hash.o: gap-if.h
gap_hash.o: gap_hash.h hash_lib.h
gap_hash.o: io_handle.h
gap_hash.o: io_utils.h IO.h
hash.o: $(SRCROOT)/Misc/xalloc.h
hash.o: hash.h
hash_lib.o: $(SRCROOT)/Misc/misc.h
hash_lib.o: $(SRCROOT)/Misc/os.h
hash_lib.o: $(SRCROOT)/Misc/xalloc.h
hash_lib.o: $(SRCROOT)/seq_utils/align_lib.h
hash_lib.o: $(SRCROOT)/seq_utils/dna_utils.h
hash_lib.o: $(PWD)/staden_config.h
hash_lib.o: hash_lib.h
init.o: $(SRCROOT)/Misc/array.h
init.o: $(SRCROOT)/Misc/bitmap.h
init.o: $(SRCROOT)/Misc/dstring.h
init.o: $(SRCROOT)/Misc/os.h
init.o: $(SRCROOT)/Misc/xalloc.h
init.o: $(SRCROOT)/Misc/xerror.h
init.o: $(SRCROOT)/g/freetree.h
init.o: $(SRCROOT)/g/g-defs.h
init.o: $(SRCROOT)/g/g-error.h
init.o: $(SRCROOT)/g/g-filedefs.h
init.o: $(SRCROOT)/g/g-os.h
init.o: $(SRCROOT)/g/g-struct.h
init.o: $(SRCROOT)/primer3/src/dpal.h
init.o: $(SRCROOT)/primer3/src/primer3.h
init.o: $(PWD)/staden_config.h
init.o: $(SRCROOT)/tk_utils/cli_arg.h
init.o: $(SRCROOT)/tk_utils/postscript.h
init.o: $(SRCROOT)/tk_utils/sheet.h
init.o: $(SRCROOT)/tk_utils/tcl_utils.h
init.o: $(SRCROOT)/tk_utils/tkSheet.h
init.o: $(SRCROOT)/tk_utils/tkSheet_common.h
init.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
init.o: $(SRCROOT)/tk_utils/tkTrace.h
init.o: IO.h
init.o: contigEditor.h fort.h
init.o: edUtils.h qual.h
init.o: fortran.h
init.o: gap-dbstruct.h
init.o: gap-error.h
init.o: gap-if.h
init.o: gap-tcl.h active_tags.h
init.o: gap_globals.h
init.o: io-reg.h primlib.h
init.o: io_handle.h edStructs.h
init.o: io_utils.h
init.o: newgap_cmds.h
init.o: tagUtils.h tagDefs.h
init.o: tagdb.h
init.o: template.h list.h undo.h
init.o: tkEdNames.h
init.o: tkEditor.h
init.o: tman_display.h
io-reg.o: $(SRCROOT)/Misc/array.h
io-reg.o: $(SRCROOT)/Misc/bitmap.h
io-reg.o: $(SRCROOT)/Misc/misc.h
io-reg.o: $(SRCROOT)/Misc/os.h
io-reg.o: $(SRCROOT)/Misc/xalloc.h
io-reg.o: $(SRCROOT)/Misc/xerror.h
io-reg.o: $(SRCROOT)/g/freetree.h
io-reg.o: $(SRCROOT)/g/g-defs.h
io-reg.o: $(SRCROOT)/g/g-error.h
io-reg.o: $(SRCROOT)/g/g-filedefs.h
io-reg.o: $(SRCROOT)/g/g-os.h
io-reg.o: $(SRCROOT)/g/g-struct.h
io-reg.o: $(PWD)/staden_config.h
io-reg.o: $(SRCROOT)/tk_utils/text_output.h
io-reg.o: IO.h
io-reg.o: gap-dbstruct.h
io-reg.o: gap-if.h
io-reg.o: io-reg.h
io-reg.o: io_handle.h io_utils.h
io-reg.o: tk-io-reg.h
io_handle.o: $(SRCROOT)/Misc/array.h
io_handle.o: $(SRCROOT)/Misc/bitmap.h
io_handle.o: $(SRCROOT)/Misc/dstring.h
io_handle.o: $(SRCROOT)/Misc/misc.h
io_handle.o: $(SRCROOT)/Misc/os.h
io_handle.o: $(SRCROOT)/Misc/xalloc.h
io_handle.o: $(SRCROOT)/Misc/xerror.h
io_handle.o: $(SRCROOT)/g/freetree.h
io_handle.o: $(SRCROOT)/g/g-defs.h
io_handle.o: $(SRCROOT)/g/g-error.h
io_handle.o: $(SRCROOT)/g/g-filedefs.h
io_handle.o: $(SRCROOT)/g/g-os.h
io_handle.o: $(SRCROOT)/g/g-struct.h
io_handle.o: $(SRCROOT)/primer3/src/dpal.h
io_handle.o: $(SRCROOT)/primer3/src/primer3.h
io_handle.o: $(PWD)/staden_config.h
io_handle.o: $(SRCROOT)/tk_utils/cli_arg.h
io_handle.o: $(SRCROOT)/tk_utils/postscript.h
io_handle.o: $(SRCROOT)/tk_utils/sheet.h
io_handle.o: $(SRCROOT)/tk_utils/text_output.h
io_handle.o: $(SRCROOT)/tk_utils/tkSheet.h
io_handle.o: $(SRCROOT)/tk_utils/tkSheet_common.h
io_handle.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
io_handle.o: $(SRCROOT)/tk_utils/tkTrace.h
io_handle.o: IO.h
io_handle.o: edStructs.h fort.h
io_handle.o: edUtils.h qual.h
io_handle.o: gap-dbstruct.h
io_handle.o: gap-if.h
io_handle.o: io_handle.h
io_handle.o: io_utils.h
io_handle.o: primlib.h
io_handle.o: tagDefs.h fortran.h
io_handle.o: template.h list.h
io_handle.o: tkEdNames.h io-reg.h
io_handle.o: tkEditor.h
io_handle.o: tman_display.h
io_handle.o: undo.h tagUtils.h
io_utils.o: $(SRCROOT)/Misc/FtoC.h
io_utils.o: $(SRCROOT)/Misc/array.h
io_utils.o: $(SRCROOT)/Misc/bitmap.h
io_utils.o: $(SRCROOT)/Misc/misc.h
io_utils.o: $(SRCROOT)/Misc/os.h
io_utils.o: $(SRCROOT)/Misc/xalloc.h
io_utils.o: $(SRCROOT)/Misc/xerror.h
io_utils.o: $(SRCROOT)/g/freetree.h
io_utils.o: $(SRCROOT)/g/g-defs.h
io_utils.o: $(SRCROOT)/g/g-error.h
io_utils.o: $(SRCROOT)/g/g-filedefs.h
io_utils.o: $(SRCROOT)/g/g-os.h
io_utils.o: $(SRCROOT)/g/g-struct.h
io_utils.o: $(PWD)/staden_config.h
io_utils.o: IO.h
io_utils.o: gap-dbstruct.h
io_utils.o: gap-if.h
io_utils.o: io_handle.h
io_utils.o: io_utils.h
join.o: $(SRCROOT)/Misc/array.h
join.o: $(SRCROOT)/Misc/bitmap.h
join.o: $(SRCROOT)/Misc/dstring.h
join.o: $(SRCROOT)/Misc/misc.h
join.o: $(SRCROOT)/Misc/os.h
join.o: $(SRCROOT)/Misc/xalloc.h
join.o: $(SRCROOT)/Misc/xerror.h
join.o: $(SRCROOT)/g/freetree.h
join.o: $(SRCROOT)/g/g-defs.h
join.o: $(SRCROOT)/g/g-error.h
join.o: $(SRCROOT)/g/g-filedefs.h
join.o: $(SRCROOT)/g/g-os.h
join.o: $(SRCROOT)/g/g-struct.h
join.o: $(SRCROOT)/primer3/src/dpal.h
join.o: $(SRCROOT)/primer3/src/primer3.h
join.o: $(SRCROOT)/seq_utils/align.h
join.o: $(SRCROOT)/seq_utils/align_lib.h
join.o: $(SRCROOT)/seq_utils/align_lib_old.h
join.o: $(SRCROOT)/seq_utils/dna_utils.h
join.o: $(PWD)/staden_config.h
join.o: $(SRCROOT)/tk_utils/cli_arg.h
join.o: $(SRCROOT)/tk_utils/postscript.h
join.o: $(SRCROOT)/tk_utils/sheet.h
join.o: $(SRCROOT)/tk_utils/tcl_utils.h
join.o: $(SRCROOT)/tk_utils/tkSheet.h
join.o: $(SRCROOT)/tk_utils/tkSheet_common.h
join.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
join.o: $(SRCROOT)/tk_utils/tkTrace.h
join.o: IO.h
join.o: IO2.h contigEditor.h
join.o: edStructs.h
join.o: edUtils.h qual.h
join.o: fort.h io_utils.h
join.o: fortran.h
join.o: gap-dbstruct.h
join.o: gap-if.h
join.o: gap_globals.h
join.o: hash_lib.h notes.h
join.o: io_handle.h tkEditor.h
join.o: primlib.h
join.o: tagUtils.h tagDefs.h
join.o: template.h list.h undo.h
join.o: tkEdNames.h io-reg.h
join.o: tman_display.h
legacy_f2c.o: f2c.h
list.o: $(SRCROOT)/Misc/xalloc.h
list.o: list.h
list_proc.o: $(SRCROOT)/Misc/FtoC.h
list_proc.o: $(SRCROOT)/Misc/array.h
list_proc.o: $(SRCROOT)/Misc/bitmap.h
list_proc.o: $(SRCROOT)/Misc/misc.h
list_proc.o: $(SRCROOT)/Misc/os.h
list_proc.o: $(SRCROOT)/Misc/xalloc.h
list_proc.o: $(SRCROOT)/Misc/xerror.h
list_proc.o: $(SRCROOT)/g/freetree.h
list_proc.o: $(SRCROOT)/g/g-defs.h
list_proc.o: $(SRCROOT)/g/g-error.h
list_proc.o: $(SRCROOT)/g/g-filedefs.h
list_proc.o: $(SRCROOT)/g/g-os.h
list_proc.o: $(SRCROOT)/g/g-struct.h
list_proc.o: $(PWD)/staden_config.h
list_proc.o: $(SRCROOT)/tk_utils/tcl_utils.h
list_proc.o: IO.h
list_proc.o: gap-dbstruct.h
list_proc.o: gap-if.h
list_proc.o: io_handle.h
list_proc.o: io_utils.h
list_proc.o: list.h gap_globals.h
list_proc.o: newgap_cmds.h
mess.o: $(SRCROOT)/tk_utils/text_output.h
newgap_cmds.o: $(SRCROOT)/Misc/array.h
newgap_cmds.o: $(SRCROOT)/Misc/bitmap.h
newgap_cmds.o: $(SRCROOT)/Misc/dstring.h
newgap_cmds.o: $(SRCROOT)/Misc/locks.h
newgap_cmds.o: $(SRCROOT)/Misc/misc.h
newgap_cmds.o: $(SRCROOT)/Misc/os.h
newgap_cmds.o: $(SRCROOT)/Misc/xalloc.h
newgap_cmds.o: $(SRCROOT)/Misc/xerror.h
newgap_cmds.o: $(SRCROOT)/g/freetree.h
newgap_cmds.o: $(SRCROOT)/g/g-defs.h
newgap_cmds.o: $(SRCROOT)/g/g-error.h
newgap_cmds.o: $(SRCROOT)/g/g-filedefs.h
newgap_cmds.o: $(SRCROOT)/g/g-os.h
newgap_cmds.o: $(SRCROOT)/g/g-struct.h
newgap_cmds.o: $(SRCROOT)/primer3/src/dpal.h
newgap_cmds.o: $(SRCROOT)/primer3/src/primer3.h
newgap_cmds.o: $(SRCROOT)/seq_utils/align.h
newgap_cmds.o: $(SRCROOT)/seq_utils/align_lib.h
newgap_cmds.o: $(SRCROOT)/seq_utils/align_lib_old.h
newgap_cmds.o: $(SRCROOT)/seq_utils/genetic_code.h
newgap_cmds.o: $(SRCROOT)/seq_utils/read_matrix.h
newgap_cmds.o: $(SRCROOT)/seq_utils/renz_utils.h
newgap_cmds.o: $(SRCROOT)/seq_utils/sequence_formats.h
newgap_cmds.o: $(PWD)/staden_config.h
newgap_cmds.o: $(SRCROOT)/tk_utils/canvas_box.h
newgap_cmds.o: $(SRCROOT)/tk_utils/cli_arg.h
newgap_cmds.o: $(SRCROOT)/tk_utils/postscript.h
newgap_cmds.o: $(SRCROOT)/tk_utils/sheet.h
newgap_cmds.o: $(SRCROOT)/tk_utils/tclXkeylist.h
newgap_cmds.o: $(SRCROOT)/tk_utils/tcl_utils.h
newgap_cmds.o: $(SRCROOT)/tk_utils/text_output.h
newgap_cmds.o: $(SRCROOT)/tk_utils/tkSheet.h
newgap_cmds.o: $(SRCROOT)/tk_utils/tkSheet_common.h
newgap_cmds.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
newgap_cmds.o: $(SRCROOT)/tk_utils/tkTrace.h
newgap_cmds.o: IO.h
newgap_cmds.o: active_tags.h
newgap_cmds.o: alter_rel.h
newgap_cmds.o: assemble_direct.h
newgap_cmds.o: auto_assemble.h
newgap_cmds.o: auto_break.h
newgap_cmds.o: break_contig.h
newgap_cmds.o: check_assembly.h
newgap_cmds.o: clip.h notes.h
newgap_cmds.o: complement.h
newgap_cmds.o: confidence_graph.h
newgap_cmds.o: consen.h
newgap_cmds.o: consistency_display.h
newgap_cmds.o: contig_order.h
newgap_cmds.o: contig_selector.h
newgap_cmds.o: copy_db.h
newgap_cmds.o: cs-object.h
newgap_cmds.o: dbcheck.h
newgap_cmds.o: dis_readings.h
newgap_cmds.o: dstrand.h
newgap_cmds.o: edStructs.h
newgap_cmds.o: extract.h preass.h
newgap_cmds.o: find_oligo.h
newgap_cmds.o: find_repeats.h
newgap_cmds.o: fort.h tkEditor.h
newgap_cmds.o: fortran.h
newgap_cmds.o: gap-create.h
newgap_cmds.o: gap-dbstruct.h
newgap_cmds.o: gap-if.h
newgap_cmds.o: gap_array.h fij.h
newgap_cmds.o: gap_canvas_box.h
newgap_cmds.o: gap_cli_arg.h
newgap_cmds.o: gap_globals.h
newgap_cmds.o: hash.h
newgap_cmds.o: hash_lib.h
newgap_cmds.o: io-reg.h
newgap_cmds.o: io_handle.h
newgap_cmds.o: io_utils.h
newgap_cmds.o: list_proc.h
newgap_cmds.o: newgap_cmds.h
newgap_cmds.o: newgap_structs.h
newgap_cmds.o: oligo_sel.h
newgap_cmds.o: plot_quality.h
newgap_cmds.o: primlib.h
newgap_cmds.o: probe.h
newgap_cmds.o: qual.h qualIO.h
newgap_cmds.o: quality_plot.h
newgap_cmds.o: reactions.h
newgap_cmds.o: reading_coverage.h
newgap_cmds.o: readpair.h
newgap_cmds.o: readpair_coverage.h
newgap_cmds.o: restriction_enzymes.h
newgap_cmds.o: ruler_display.h
newgap_cmds.o: seqInfo.h
newgap_cmds.o: show_relationships.h
newgap_cmds.o: shuffle_pads.h
newgap_cmds.o: stop_codon.h
newgap_cmds.o: strand_coverage.h
newgap_cmds.o: tagDefs.h
newgap_cmds.o: tagUtils.h
newgap_cmds.o: tagdb.h
newgap_cmds.o: template.h list.h
newgap_cmds.o: template_display.h
newgap_cmds.o: tk-io-reg.h
newgap_cmds.o: tkEdNames.h
newgap_cmds.o: tman_display.h
newgap_cmds.o: undo.h edUtils.h
notedb.o: $(SRCROOT)/Misc/misc.h
notedb.o: $(SRCROOT)/Misc/os.h
notedb.o: $(SRCROOT)/Misc/parse_db.h
notedb.o: $(SRCROOT)/Misc/xalloc.h
notedb.o: $(PWD)/staden_config.h
notedb.o: notedb.h
notes.o: $(SRCROOT)/Misc/array.h
notes.o: $(SRCROOT)/Misc/bitmap.h
notes.o: $(SRCROOT)/Misc/dstring.h
notes.o: $(SRCROOT)/Misc/os.h
notes.o: $(SRCROOT)/Misc/xalloc.h
notes.o: $(SRCROOT)/Misc/xerror.h
notes.o: $(SRCROOT)/g/freetree.h
notes.o: $(SRCROOT)/g/g-defs.h
notes.o: $(SRCROOT)/g/g-error.h
notes.o: $(SRCROOT)/g/g-filedefs.h
notes.o: $(SRCROOT)/g/g-os.h
notes.o: $(SRCROOT)/g/g-struct.h
notes.o: $(SRCROOT)/primer3/src/dpal.h
notes.o: $(SRCROOT)/primer3/src/primer3.h
notes.o: $(PWD)/staden_config.h
notes.o: $(SRCROOT)/tk_utils/cli_arg.h
notes.o: $(SRCROOT)/tk_utils/postscript.h
notes.o: $(SRCROOT)/tk_utils/sheet.h
notes.o: $(SRCROOT)/tk_utils/tcl_utils.h
notes.o: $(SRCROOT)/tk_utils/text_output.h
notes.o: $(SRCROOT)/tk_utils/tkSheet.h
notes.o: $(SRCROOT)/tk_utils/tkSheet_common.h
notes.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
notes.o: $(SRCROOT)/tk_utils/tkTrace.h
notes.o: IO.h
notes.o: edUtils.h qual.h
notes.o: fort.h tkEditor.h
notes.o: gap-dbstruct.h
notes.o: gap-if.h
notes.o: gap_cli_arg.h notes.h
notes.o: gap_globals.h
notes.o: io_handle.h io_utils.h
notes.o: primlib.h
notes.o: tagDefs.h fortran.h
notes.o: tagUtils.h edStructs.h
notes.o: template.h list.h undo.h
notes.o: tkEdNames.h io-reg.h
notes.o: tman_display.h
oligo.o: $(SRCROOT)/Misc/array.h
oligo.o: $(SRCROOT)/Misc/bitmap.h
oligo.o: $(SRCROOT)/Misc/dstring.h
oligo.o: $(SRCROOT)/Misc/misc.h
oligo.o: $(SRCROOT)/Misc/os.h
oligo.o: $(SRCROOT)/Misc/reg_exp.h
oligo.o: $(SRCROOT)/Misc/xalloc.h
oligo.o: $(SRCROOT)/Misc/xerror.h
oligo.o: $(SRCROOT)/g/freetree.h
oligo.o: $(SRCROOT)/g/g-defs.h
oligo.o: $(SRCROOT)/g/g-error.h
oligo.o: $(SRCROOT)/g/g-filedefs.h
oligo.o: $(SRCROOT)/g/g-os.h
oligo.o: $(SRCROOT)/g/g-struct.h
oligo.o: $(SRCROOT)/primer3/src/dpal.h
oligo.o: $(SRCROOT)/primer3/src/primer3.h
oligo.o: $(SRCROOT)/seq_utils/dna_utils.h
oligo.o: $(PWD)/staden_config.h
oligo.o: $(SRCROOT)/tk_utils/cli_arg.h
oligo.o: $(SRCROOT)/tk_utils/postscript.h
oligo.o: $(SRCROOT)/tk_utils/sheet.h
oligo.o: $(SRCROOT)/tk_utils/tkSheet.h
oligo.o: $(SRCROOT)/tk_utils/tkSheet_common.h
oligo.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
oligo.o: $(SRCROOT)/tk_utils/tkTrace.h
oligo.o: IO.h
oligo.o: edStructs.h
oligo.o: edUtils.h qual.h
oligo.o: fort.h io_utils.h
oligo.o: fortran.h
oligo.o: gap-dbstruct.h
oligo.o: gap-if.h
oligo.o: io_handle.h tkEditor.h
oligo.o: oligo.h
oligo.o: primlib.h
oligo.o: tagUtils.h tagDefs.h
oligo.o: template.h list.h undo.h
oligo.o: tkEdNames.h io-reg.h
oligo.o: tman_display.h
oligo_sel.o: $(SRCROOT)/Misc/FtoC.h
oligo_sel.o: $(SRCROOT)/Misc/array.h
oligo_sel.o: $(SRCROOT)/Misc/bitmap.h
oligo_sel.o: $(SRCROOT)/Misc/dstring.h
oligo_sel.o: $(SRCROOT)/Misc/misc.h
oligo_sel.o: $(SRCROOT)/Misc/os.h
oligo_sel.o: $(SRCROOT)/Misc/xalloc.h
oligo_sel.o: $(SRCROOT)/Misc/xerror.h
oligo_sel.o: $(SRCROOT)/g/freetree.h
oligo_sel.o: $(SRCROOT)/g/g-defs.h
oligo_sel.o: $(SRCROOT)/g/g-error.h
oligo_sel.o: $(SRCROOT)/g/g-filedefs.h
oligo_sel.o: $(SRCROOT)/g/g-os.h
oligo_sel.o: $(SRCROOT)/g/g-struct.h
oligo_sel.o: $(SRCROOT)/primer3/src/dpal.h
oligo_sel.o: $(SRCROOT)/primer3/src/primer3.h
oligo_sel.o: $(PWD)/staden_config.h
oligo_sel.o: $(SRCROOT)/tk_utils/cli_arg.h
oligo_sel.o: $(SRCROOT)/tk_utils/postscript.h
oligo_sel.o: $(SRCROOT)/tk_utils/sheet.h
oligo_sel.o: $(SRCROOT)/tk_utils/tcl_utils.h
oligo_sel.o: $(SRCROOT)/tk_utils/text_output.h
oligo_sel.o: $(SRCROOT)/tk_utils/tkSheet.h
oligo_sel.o: $(SRCROOT)/tk_utils/tkSheet_common.h
oligo_sel.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
oligo_sel.o: $(SRCROOT)/tk_utils/tkTrace.h
oligo_sel.o: IO.h
oligo_sel.o: dstrand.h
oligo_sel.o: edStructs.h
oligo_sel.o: fort.h tkEditor.h
oligo_sel.o: gap-dbstruct.h
oligo_sel.o: gap-if.h
oligo_sel.o: gap_globals.h
oligo_sel.o: io_handle.h qual.h
oligo_sel.o: io_utils.h
oligo_sel.o: list_proc.h
oligo_sel.o: oligo_sel.h
oligo_sel.o: primlib.h
oligo_sel.o: tagDefs.h fortran.h
oligo_sel.o: tagUtils.h
oligo_sel.o: template.h list.h
oligo_sel.o: tkEdNames.h io-reg.h
oligo_sel.o: tman_display.h
oligo_sel.o: undo.h edUtils.h
parse_ft.o: $(SRCROOT)/Misc/xalloc.h
parse_ft.o: parse_ft.h
preass.o: $(SRCROOT)/Misc/array.h
preass.o: $(SRCROOT)/Misc/bitmap.h
preass.o: $(SRCROOT)/Misc/dstring.h
preass.o: $(SRCROOT)/Misc/misc.h
preass.o: $(SRCROOT)/Misc/os.h
preass.o: $(SRCROOT)/Misc/xalloc.h
preass.o: $(SRCROOT)/Misc/xerror.h
preass.o: $(SRCROOT)/g/freetree.h
preass.o: $(SRCROOT)/g/g-defs.h
preass.o: $(SRCROOT)/g/g-error.h
preass.o: $(SRCROOT)/g/g-filedefs.h
preass.o: $(SRCROOT)/g/g-os.h
preass.o: $(SRCROOT)/g/g-struct.h
preass.o: $(SRCROOT)/primer3/src/dpal.h
preass.o: $(SRCROOT)/primer3/src/primer3.h
preass.o: $(PWD)/staden_config.h
preass.o: $(SRCROOT)/tk_utils/cli_arg.h
preass.o: $(SRCROOT)/tk_utils/postscript.h
preass.o: $(SRCROOT)/tk_utils/sheet.h
preass.o: $(SRCROOT)/tk_utils/text_output.h
preass.o: $(SRCROOT)/tk_utils/tkSheet.h
preass.o: $(SRCROOT)/tk_utils/tkSheet_common.h
preass.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
preass.o: $(SRCROOT)/tk_utils/tkTrace.h
preass.o: IO.h
preass.o: clones.h dbcheck.h
preass.o: fort.h tkEditor.h
preass.o: gap-dbstruct.h
preass.o: gap-if.h
preass.o: io_handle.h io_utils.h
preass.o: preass.h
preass.o: primlib.h
preass.o: seqInfo.h
preass.o: tagDefs.h fortran.h
preass.o: tagUtils.h edStructs.h
preass.o: template.h list.h
preass.o: tkEdNames.h io-reg.h
preass.o: tman_display.h
preass.o: undo.h edUtils.h qual.h
primlib.o: $(SRCROOT)/Misc/os.h
primlib.o: $(SRCROOT)/primer3/src/dpal.h
primlib.o: $(SRCROOT)/primer3/src/primer3.h
primlib.o: $(PWD)/staden_config.h
primlib.o: $(SRCROOT)/tk_utils/tcl_utils.h
primlib.o: gap_globals.h
primlib.o: primlib.h
probe.o: $(SRCROOT)/Misc/array.h
probe.o: $(SRCROOT)/Misc/bitmap.h
probe.o: $(SRCROOT)/Misc/misc.h
probe.o: $(SRCROOT)/Misc/os.h
probe.o: $(SRCROOT)/Misc/xalloc.h
probe.o: $(SRCROOT)/Misc/xerror.h
probe.o: $(SRCROOT)/g/freetree.h
probe.o: $(SRCROOT)/g/g-defs.h
probe.o: $(SRCROOT)/g/g-error.h
probe.o: $(SRCROOT)/g/g-filedefs.h
probe.o: $(SRCROOT)/g/g-os.h
probe.o: $(SRCROOT)/g/g-struct.h
probe.o: $(SRCROOT)/primer3/src/dpal.h
probe.o: $(SRCROOT)/primer3/src/primer3.h
probe.o: $(SRCROOT)/seq_utils/align.h
probe.o: $(SRCROOT)/seq_utils/align_lib_old.h
probe.o: $(SRCROOT)/seq_utils/dna_utils.h
probe.o: $(SRCROOT)/seq_utils/sequence_formats.h
probe.o: $(PWD)/staden_config.h
probe.o: $(SRCROOT)/tk_utils/tcl_utils.h
probe.o: $(SRCROOT)/tk_utils/text_output.h
probe.o: IO.h
probe.o: assemble_direct.h
probe.o: fort.h primlib.h
probe.o: gap-dbstruct.h
probe.o: gap-if.h
probe.o: gap_globals.h
probe.o: io_handle.h io_utils.h
probe.o: probe.h
probe.o: seqInfo.h
qual.o: $(SRCROOT)/Misc/misc.h
qual.o: $(SRCROOT)/Misc/os.h
qual.o: $(SRCROOT)/Misc/xalloc.h
qual.o: $(PWD)/staden_config.h
qual.o: $(SRCROOT)/tk_utils/tcl_utils.h
qual.o: gap_globals.h
qual.o: qual.h
qual.o: qualP.h
qualIO.o: $(SRCROOT)/Misc/array.h
qualIO.o: $(SRCROOT)/Misc/bitmap.h
qualIO.o: $(SRCROOT)/Misc/misc.h
qualIO.o: $(SRCROOT)/Misc/os.h
qualIO.o: $(SRCROOT)/Misc/xalloc.h
qualIO.o: $(SRCROOT)/Misc/xerror.h
qualIO.o: $(SRCROOT)/g/freetree.h
qualIO.o: $(SRCROOT)/g/g-defs.h
qualIO.o: $(SRCROOT)/g/g-error.h
qualIO.o: $(SRCROOT)/g/g-filedefs.h
qualIO.o: $(SRCROOT)/g/g-os.h
qualIO.o: $(SRCROOT)/g/g-struct.h
qualIO.o: $(PWD)/staden_config.h
qualIO.o: $(SRCROOT)/tk_utils/tcl_utils.h
qualIO.o: IO.h
qualIO.o: IO2.h qual.h qualP.h
qualIO.o: gap-dbstruct.h
qualIO.o: gap-if.h
qualIO.o: gap_globals.h
qualIO.o: io_handle.h io_utils.h
quality_plot.o: $(SRCROOT)/Misc/array.h
quality_plot.o: $(SRCROOT)/Misc/bitmap.h
quality_plot.o: $(SRCROOT)/Misc/dstring.h
quality_plot.o: $(SRCROOT)/Misc/misc.h
quality_plot.o: $(SRCROOT)/Misc/os.h
quality_plot.o: $(SRCROOT)/Misc/xalloc.h
quality_plot.o: $(SRCROOT)/Misc/xerror.h
quality_plot.o: $(SRCROOT)/g/freetree.h
quality_plot.o: $(SRCROOT)/g/g-defs.h
quality_plot.o: $(SRCROOT)/g/g-error.h
quality_plot.o: $(SRCROOT)/g/g-filedefs.h
quality_plot.o: $(SRCROOT)/g/g-os.h
quality_plot.o: $(SRCROOT)/g/g-struct.h
quality_plot.o: $(SRCROOT)/primer3/src/dpal.h
quality_plot.o: $(SRCROOT)/primer3/src/primer3.h
quality_plot.o: $(PWD)/staden_config.h
quality_plot.o: $(SRCROOT)/tk_utils/canvas_box.h
quality_plot.o: $(SRCROOT)/tk_utils/cli_arg.h
quality_plot.o: $(SRCROOT)/tk_utils/postscript.h
quality_plot.o: $(SRCROOT)/tk_utils/sheet.h
quality_plot.o: $(SRCROOT)/tk_utils/tclXkeylist.h
quality_plot.o: $(SRCROOT)/tk_utils/tcl_utils.h
quality_plot.o: $(SRCROOT)/tk_utils/text_output.h
quality_plot.o: $(SRCROOT)/tk_utils/tkSheet.h
quality_plot.o: $(SRCROOT)/tk_utils/tkSheet_common.h
quality_plot.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
quality_plot.o: $(SRCROOT)/tk_utils/tkTrace.h
quality_plot.o: IO.h
quality_plot.o: edStructs.h
quality_plot.o: edUtils.h
quality_plot.o: fort.h
quality_plot.o: fortran.h
quality_plot.o: gap-dbstruct.h
quality_plot.o: gap-if.h
quality_plot.o: gap_canvas_box.h
quality_plot.o: gap_globals.h
quality_plot.o: hash.h template.h
quality_plot.o: io-reg.h
quality_plot.o: io_handle.h
quality_plot.o: io_utils.h
quality_plot.o: list.h
quality_plot.o: newgap_cmds.h
quality_plot.o: primlib.h
quality_plot.o: qual.h
quality_plot.o: quality_plot.h
quality_plot.o: ruler_display.h
quality_plot.o: tagDefs.h
quality_plot.o: template_display.h
quality_plot.o: tkEdNames.h
quality_plot.o: tkEditor.h
quality_plot.o: tman_display.h
quality_plot.o: undo.h tagUtils.h
reactions.o: $(SRCROOT)/Misc/array.h
reactions.o: $(SRCROOT)/Misc/bitmap.h
reactions.o: $(SRCROOT)/Misc/dstring.h
reactions.o: $(SRCROOT)/Misc/misc.h
reactions.o: $(SRCROOT)/Misc/os.h
reactions.o: $(SRCROOT)/Misc/xalloc.h
reactions.o: $(SRCROOT)/Misc/xerror.h
reactions.o: $(SRCROOT)/g/freetree.h
reactions.o: $(SRCROOT)/g/g-defs.h
reactions.o: $(SRCROOT)/g/g-error.h
reactions.o: $(SRCROOT)/g/g-filedefs.h
reactions.o: $(SRCROOT)/g/g-os.h
reactions.o: $(SRCROOT)/g/g-struct.h
reactions.o: $(SRCROOT)/primer3/src/dpal.h
reactions.o: $(SRCROOT)/primer3/src/primer3.h
reactions.o: $(PWD)/staden_config.h
reactions.o: $(SRCROOT)/tk_utils/cli_arg.h
reactions.o: $(SRCROOT)/tk_utils/postscript.h
reactions.o: $(SRCROOT)/tk_utils/sheet.h
reactions.o: $(SRCROOT)/tk_utils/tcl_utils.h
reactions.o: $(SRCROOT)/tk_utils/tkSheet.h
reactions.o: $(SRCROOT)/tk_utils/tkSheet_common.h
reactions.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
reactions.o: $(SRCROOT)/tk_utils/tkTrace.h
reactions.o: IO.h
reactions.o: edStructs.h
reactions.o: gap-dbstruct.h
reactions.o: gap-if.h
reactions.o: gap_globals.h
reactions.o: io_handle.h fort.h
reactions.o: io_utils.h
reactions.o: primlib.h
reactions.o: qual.h qualP.h
reactions.o: tagDefs.h fortran.h
reactions.o: tagUtils.h
reactions.o: template.h list.h
reactions.o: tkEdNames.h io-reg.h
reactions.o: tkEditor.h
reactions.o: tman_display.h
reactions.o: undo.h edUtils.h
reading_coverage.o: $(SRCROOT)/Misc/array.h
reading_coverage.o: $(SRCROOT)/Misc/bitmap.h
reading_coverage.o: $(SRCROOT)/Misc/misc.h
reading_coverage.o: $(SRCROOT)/Misc/os.h
reading_coverage.o: $(SRCROOT)/Misc/xalloc.h
reading_coverage.o: $(SRCROOT)/Misc/xerror.h
reading_coverage.o: $(SRCROOT)/g/freetree.h
reading_coverage.o: $(SRCROOT)/g/g-defs.h
reading_coverage.o: $(SRCROOT)/g/g-error.h
reading_coverage.o: $(SRCROOT)/g/g-filedefs.h
reading_coverage.o: $(SRCROOT)/g/g-os.h
reading_coverage.o: $(SRCROOT)/g/g-struct.h
reading_coverage.o: $(PWD)/staden_config.h
reading_coverage.o: $(SRCROOT)/tk_utils/canvas_box.h
reading_coverage.o: $(SRCROOT)/tk_utils/tclXkeylist.h
reading_coverage.o: $(SRCROOT)/tk_utils/tcl_utils.h
reading_coverage.o: $(SRCROOT)/tk_utils/text_output.h
reading_coverage.o: consen.h
reading_coverage.o: consistency_display.h
reading_coverage.o: gap-dbstruct.h
reading_coverage.o: gap-if.h
reading_coverage.o: gap_canvas_box.h
reading_coverage.o: gap_globals.h
reading_coverage.o: hash.h list.h
reading_coverage.o: io-reg.h IO.h
reading_coverage.o: io_handle.h
reading_coverage.o: io_utils.h
reading_coverage.o: newgap_cmds.h
reading_coverage.o: qual.h
reading_coverage.o: reading_coverage.h
reading_coverage.o: ruler_display.h
reading_coverage.o: template.h
reading_coverage.o: template_display.h
readpair.o: $(SRCROOT)/Misc/array.h
readpair.o: $(SRCROOT)/Misc/bitmap.h
readpair.o: $(SRCROOT)/Misc/dstring.h
readpair.o: $(SRCROOT)/Misc/misc.h
readpair.o: $(SRCROOT)/Misc/os.h
readpair.o: $(SRCROOT)/Misc/xalloc.h
readpair.o: $(SRCROOT)/Misc/xerror.h
readpair.o: $(SRCROOT)/g/freetree.h
readpair.o: $(SRCROOT)/g/g-defs.h
readpair.o: $(SRCROOT)/g/g-error.h
readpair.o: $(SRCROOT)/g/g-filedefs.h
readpair.o: $(SRCROOT)/g/g-os.h
readpair.o: $(SRCROOT)/g/g-struct.h
readpair.o: $(SRCROOT)/primer3/src/dpal.h
readpair.o: $(SRCROOT)/primer3/src/primer3.h
readpair.o: $(PWD)/staden_config.h
readpair.o: $(SRCROOT)/tk_utils/canvas_box.h
readpair.o: $(SRCROOT)/tk_utils/cli_arg.h
readpair.o: $(SRCROOT)/tk_utils/postscript.h
readpair.o: $(SRCROOT)/tk_utils/sheet.h
readpair.o: $(SRCROOT)/tk_utils/tclXkeylist.h
readpair.o: $(SRCROOT)/tk_utils/tcl_utils.h
readpair.o: $(SRCROOT)/tk_utils/text_output.h
readpair.o: $(SRCROOT)/tk_utils/tkSheet.h
readpair.o: $(SRCROOT)/tk_utils/tkSheet_common.h
readpair.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
readpair.o: $(SRCROOT)/tk_utils/tkTrace.h
readpair.o: IO.h
readpair.o: complement.h
readpair.o: contigEditor.h
readpair.o: contig_selector.h
readpair.o: edStructs.h
readpair.o: edUtils.h qual.h
readpair.o: fortran.h fort.h
readpair.o: gap-dbstruct.h
readpair.o: gap-if.h
readpair.o: io-reg.h cs-object.h
readpair.o: io_handle.h
readpair.o: io_utils.h
readpair.o: list.h gap_globals.h
readpair.o: newgap_cmds.h
readpair.o: readpair.h
readpair.o: tagUtils.h tagDefs.h
readpair.o: template.h undo.h
readpair.o: tkEdNames.h primlib.h
readpair.o: tkEditor.h
readpair.o: tman_display.h
readpair_coverage.o: $(SRCROOT)/Misc/array.h
readpair_coverage.o: $(SRCROOT)/Misc/bitmap.h
readpair_coverage.o: $(SRCROOT)/Misc/misc.h
readpair_coverage.o: $(SRCROOT)/Misc/os.h
readpair_coverage.o: $(SRCROOT)/Misc/xalloc.h
readpair_coverage.o: $(SRCROOT)/Misc/xerror.h
readpair_coverage.o: $(SRCROOT)/g/freetree.h
readpair_coverage.o: $(SRCROOT)/g/g-defs.h
readpair_coverage.o: $(SRCROOT)/g/g-error.h
readpair_coverage.o: $(SRCROOT)/g/g-filedefs.h
readpair_coverage.o: $(SRCROOT)/g/g-os.h
readpair_coverage.o: $(SRCROOT)/g/g-struct.h
readpair_coverage.o: $(PWD)/staden_config.h
readpair_coverage.o: $(SRCROOT)/tk_utils/canvas_box.h
readpair_coverage.o: $(SRCROOT)/tk_utils/tclXkeylist.h
readpair_coverage.o: $(SRCROOT)/tk_utils/tcl_utils.h
readpair_coverage.o: $(SRCROOT)/tk_utils/text_output.h
readpair_coverage.o: IO.h
readpair_coverage.o: consen.h
readpair_coverage.o: consistency_display.h
readpair_coverage.o: gap-dbstruct.h
readpair_coverage.o: gap-if.h
readpair_coverage.o: gap_canvas_box.h
readpair_coverage.o: gap_globals.h
readpair_coverage.o: hash.h
readpair_coverage.o: io-reg.h
readpair_coverage.o: io_handle.h
readpair_coverage.o: io_utils.h
readpair_coverage.o: list.h
readpair_coverage.o: newgap_cmds.h
readpair_coverage.o: qual.h
readpair_coverage.o: reading_coverage.h
readpair_coverage.o: readpair_coverage.h
readpair_coverage.o: ruler_display.h
readpair_coverage.o: template.h
readpair_coverage.o: template_display.h
restriction_enzymes.o: $(SRCROOT)/Misc/array.h
restriction_enzymes.o: $(SRCROOT)/Misc/bitmap.h
restriction_enzymes.o: $(SRCROOT)/Misc/dstring.h
restriction_enzymes.o: $(SRCROOT)/Misc/getfile.h
restriction_enzymes.o: $(SRCROOT)/Misc/misc.h
restriction_enzymes.o: $(SRCROOT)/Misc/os.h
restriction_enzymes.o: $(SRCROOT)/Misc/xalloc.h
restriction_enzymes.o: $(SRCROOT)/Misc/xerror.h
restriction_enzymes.o: $(SRCROOT)/g/freetree.h
restriction_enzymes.o: $(SRCROOT)/g/g-defs.h
restriction_enzymes.o: $(SRCROOT)/g/g-error.h
restriction_enzymes.o: $(SRCROOT)/g/g-filedefs.h
restriction_enzymes.o: $(SRCROOT)/g/g-os.h
restriction_enzymes.o: $(SRCROOT)/g/g-struct.h
restriction_enzymes.o: $(SRCROOT)/primer3/src/dpal.h
restriction_enzymes.o: $(SRCROOT)/primer3/src/primer3.h
restriction_enzymes.o: $(SRCROOT)/seq_utils/dna_utils.h
restriction_enzymes.o: $(SRCROOT)/seq_utils/renz_utils.h
restriction_enzymes.o: $(PWD)/staden_config.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/canvas_box.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/cli_arg.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/postscript.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/restriction_enzyme_map.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/sheet.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/tclXkeylist.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/tcl_utils.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/text_output.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/tkSheet.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/tkSheet_common.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
restriction_enzymes.o: $(SRCROOT)/tk_utils/tkTrace.h
restriction_enzymes.o: IO.h
restriction_enzymes.o: edStructs.h
restriction_enzymes.o: edUtils.h
restriction_enzymes.o: fort.h
restriction_enzymes.o: fortran.h
restriction_enzymes.o: gap-dbstruct.h
restriction_enzymes.o: gap-if.h
restriction_enzymes.o: gap_canvas_box.h
restriction_enzymes.o: gap_globals.h
restriction_enzymes.o: hash.h
restriction_enzymes.o: io-reg.h
restriction_enzymes.o: io_handle.h
restriction_enzymes.o: io_utils.h
restriction_enzymes.o: list.h
restriction_enzymes.o: list_proc.h
restriction_enzymes.o: newgap_cmds.h
restriction_enzymes.o: primlib.h
restriction_enzymes.o: qual.h
restriction_enzymes.o: restriction_enzymes.h
restriction_enzymes.o: ruler_display.h
restriction_enzymes.o: tagDefs.h
restriction_enzymes.o: tagUtils.h
restriction_enzymes.o: template.h
restriction_enzymes.o: template_display.h
restriction_enzymes.o: tkEdNames.h
restriction_enzymes.o: tkEditor.h
restriction_enzymes.o: tman_display.h
restriction_enzymes.o: undo.h
ruler_display.o: $(SRCROOT)/Misc/array.h
ruler_display.o: $(SRCROOT)/Misc/bitmap.h
ruler_display.o: $(SRCROOT)/Misc/misc.h
ruler_display.o: $(SRCROOT)/Misc/os.h
ruler_display.o: $(SRCROOT)/Misc/xalloc.h
ruler_display.o: $(SRCROOT)/Misc/xerror.h
ruler_display.o: $(SRCROOT)/g/freetree.h
ruler_display.o: $(SRCROOT)/g/g-defs.h
ruler_display.o: $(SRCROOT)/g/g-error.h
ruler_display.o: $(SRCROOT)/g/g-filedefs.h
ruler_display.o: $(SRCROOT)/g/g-os.h
ruler_display.o: $(SRCROOT)/g/g-struct.h
ruler_display.o: $(PWD)/staden_config.h
ruler_display.o: $(SRCROOT)/tk_utils/canvas_box.h
ruler_display.o: $(SRCROOT)/tk_utils/ruler_tick.h
ruler_display.o: $(SRCROOT)/tk_utils/tclXkeylist.h
ruler_display.o: $(SRCROOT)/tk_utils/tcl_utils.h
ruler_display.o: IO.h
ruler_display.o: gap-dbstruct.h
ruler_display.o: gap-if.h
ruler_display.o: gap_globals.h
ruler_display.o: io-reg.h hash.h
ruler_display.o: io_handle.h
ruler_display.o: io_utils.h
ruler_display.o: list.h
ruler_display.o: ruler_display.h
ruler_display.o: template.h
ruler_display.o: template_display.h
searchUtils.o: $(SRCROOT)/Misc/array.h
searchUtils.o: $(SRCROOT)/Misc/bitmap.h
searchUtils.o: $(SRCROOT)/Misc/dstring.h
searchUtils.o: $(SRCROOT)/Misc/misc.h
searchUtils.o: $(SRCROOT)/Misc/os.h
searchUtils.o: $(SRCROOT)/Misc/reg_exp.h
searchUtils.o: $(SRCROOT)/Misc/xalloc.h
searchUtils.o: $(SRCROOT)/Misc/xerror.h
searchUtils.o: $(SRCROOT)/g/freetree.h
searchUtils.o: $(SRCROOT)/g/g-defs.h
searchUtils.o: $(SRCROOT)/g/g-error.h
searchUtils.o: $(SRCROOT)/g/g-filedefs.h
searchUtils.o: $(SRCROOT)/g/g-os.h
searchUtils.o: $(SRCROOT)/g/g-struct.h
searchUtils.o: $(SRCROOT)/primer3/src/dpal.h
searchUtils.o: $(SRCROOT)/primer3/src/primer3.h
searchUtils.o: $(SRCROOT)/seq_utils/dna_utils.h
searchUtils.o: $(SRCROOT)/seq_utils/search_utils.h
searchUtils.o: $(PWD)/staden_config.h
searchUtils.o: $(SRCROOT)/tk_utils/cli_arg.h
searchUtils.o: $(SRCROOT)/tk_utils/postscript.h
searchUtils.o: $(SRCROOT)/tk_utils/sheet.h
searchUtils.o: $(SRCROOT)/tk_utils/text_output.h
searchUtils.o: $(SRCROOT)/tk_utils/tkSheet.h
searchUtils.o: $(SRCROOT)/tk_utils/tkSheet_common.h
searchUtils.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
searchUtils.o: $(SRCROOT)/tk_utils/tkTrace.h
searchUtils.o: IO.h
searchUtils.o: active_tags.h
searchUtils.o: contigEditor.h
searchUtils.o: edStructs.h
searchUtils.o: edUtils.h
searchUtils.o: fort.h io_utils.h
searchUtils.o: fortran.h
searchUtils.o: gap-dbstruct.h
searchUtils.o: gap-if.h
searchUtils.o: io-reg.h primlib.h
searchUtils.o: io_handle.h
searchUtils.o: qual.h
searchUtils.o: select.h
searchUtils.o: tagDefs.h
searchUtils.o: tagdb.h
searchUtils.o: template.h list.h
searchUtils.o: tkEdNames.h
searchUtils.o: tkEditor.h
searchUtils.o: tman_display.h
searchUtils.o: tman_interface.h
searchUtils.o: undo.h tagUtils.h
seqInfo.o: $(SRCROOT)/Misc/FtoC.h
seqInfo.o: $(SRCROOT)/Misc/array.h
seqInfo.o: $(SRCROOT)/Misc/bitmap.h
seqInfo.o: $(SRCROOT)/Misc/dstring.h
seqInfo.o: $(SRCROOT)/Misc/misc.h
seqInfo.o: $(SRCROOT)/Misc/os.h
seqInfo.o: $(SRCROOT)/Misc/xalloc.h
seqInfo.o: $(SRCROOT)/Misc/xerror.h
seqInfo.o: $(SRCROOT)/g/freetree.h
seqInfo.o: $(SRCROOT)/g/g-defs.h
seqInfo.o: $(SRCROOT)/g/g-error.h
seqInfo.o: $(SRCROOT)/g/g-filedefs.h
seqInfo.o: $(SRCROOT)/g/g-os.h
seqInfo.o: $(SRCROOT)/g/g-struct.h
seqInfo.o: $(SRCROOT)/primer3/src/dpal.h
seqInfo.o: $(SRCROOT)/primer3/src/primer3.h
seqInfo.o: $(PWD)/staden_config.h
seqInfo.o: $(SRCROOT)/tk_utils/cli_arg.h
seqInfo.o: $(SRCROOT)/tk_utils/postscript.h
seqInfo.o: $(SRCROOT)/tk_utils/sheet.h
seqInfo.o: $(SRCROOT)/tk_utils/tkSheet.h
seqInfo.o: $(SRCROOT)/tk_utils/tkSheet_common.h
seqInfo.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
seqInfo.o: $(SRCROOT)/tk_utils/tkTrace.h
seqInfo.o: IO.h
seqInfo.o: edStructs.h fort.h
seqInfo.o: gap-dbstruct.h
seqInfo.o: gap-if.h
seqInfo.o: io_handle.h io_utils.h
seqInfo.o: parse_ft.h tagUtils.h
seqInfo.o: primlib.h
seqInfo.o: qual.h tman_display.h
seqInfo.o: seqInfo.h
seqInfo.o: tagDefs.h fortran.h
seqInfo.o: tagdb.h
seqInfo.o: template.h list.h
seqInfo.o: tkEdNames.h io-reg.h
seqInfo.o: tkEditor.h
seqInfo.o: undo.h edUtils.h
show_relationships.o: $(SRCROOT)/Misc/array.h
show_relationships.o: $(SRCROOT)/Misc/bitmap.h
show_relationships.o: $(SRCROOT)/Misc/misc.h
show_relationships.o: $(SRCROOT)/Misc/os.h
show_relationships.o: $(SRCROOT)/Misc/xalloc.h
show_relationships.o: $(SRCROOT)/Misc/xerror.h
show_relationships.o: $(SRCROOT)/g/freetree.h
show_relationships.o: $(SRCROOT)/g/g-defs.h
show_relationships.o: $(SRCROOT)/g/g-error.h
show_relationships.o: $(SRCROOT)/g/g-filedefs.h
show_relationships.o: $(SRCROOT)/g/g-os.h
show_relationships.o: $(SRCROOT)/g/g-struct.h
show_relationships.o: $(PWD)/staden_config.h
show_relationships.o: $(SRCROOT)/tk_utils/cli_arg.h
show_relationships.o: $(SRCROOT)/tk_utils/text_output.h
show_relationships.o: IO.h
show_relationships.o: gap-dbstruct.h
show_relationships.o: gap-if.h
show_relationships.o: io_handle.h
show_relationships.o: io_utils.h
show_relationships.o: newgap_structs.h
shuffle_pads.o: $(SRCROOT)/Misc/array.h
shuffle_pads.o: $(SRCROOT)/Misc/bitmap.h
shuffle_pads.o: $(SRCROOT)/Misc/dstring.h
shuffle_pads.o: $(SRCROOT)/Misc/misc.h
shuffle_pads.o: $(SRCROOT)/Misc/os.h
shuffle_pads.o: $(SRCROOT)/Misc/xalloc.h
shuffle_pads.o: $(SRCROOT)/Misc/xerror.h
shuffle_pads.o: $(SRCROOT)/g/freetree.h
shuffle_pads.o: $(SRCROOT)/g/g-defs.h
shuffle_pads.o: $(SRCROOT)/g/g-error.h
shuffle_pads.o: $(SRCROOT)/g/g-filedefs.h
shuffle_pads.o: $(SRCROOT)/g/g-os.h
shuffle_pads.o: $(SRCROOT)/g/g-struct.h
shuffle_pads.o: $(SRCROOT)/primer3/src/dpal.h
shuffle_pads.o: $(SRCROOT)/primer3/src/primer3.h
shuffle_pads.o: $(SRCROOT)/seq_utils/align.h
shuffle_pads.o: $(SRCROOT)/seq_utils/align_lib.h
shuffle_pads.o: $(SRCROOT)/seq_utils/align_lib_old.h
shuffle_pads.o: $(SRCROOT)/seq_utils/dna_utils.h
shuffle_pads.o: $(PWD)/staden_config.h
shuffle_pads.o: $(SRCROOT)/tk_utils/cli_arg.h
shuffle_pads.o: $(SRCROOT)/tk_utils/postscript.h
shuffle_pads.o: $(SRCROOT)/tk_utils/sheet.h
shuffle_pads.o: $(SRCROOT)/tk_utils/text_output.h
shuffle_pads.o: $(SRCROOT)/tk_utils/tkSheet.h
shuffle_pads.o: $(SRCROOT)/tk_utils/tkSheet_common.h
shuffle_pads.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
shuffle_pads.o: $(SRCROOT)/tk_utils/tkTrace.h
shuffle_pads.o: IO.h
shuffle_pads.o: dis_readings.h
shuffle_pads.o: edStructs.h
shuffle_pads.o: fort.h tkEditor.h
shuffle_pads.o: fortran.h
shuffle_pads.o: gap-dbstruct.h
shuffle_pads.o: gap-if.h
shuffle_pads.o: io-reg.h
shuffle_pads.o: io_handle.h
shuffle_pads.o: io_utils.h
shuffle_pads.o: primlib.h
shuffle_pads.o: qual.h
shuffle_pads.o: shuffle_pads.h
shuffle_pads.o: tagDefs.h
shuffle_pads.o: tagUtils.h
shuffle_pads.o: template.h list.h
shuffle_pads.o: tkEdNames.h
shuffle_pads.o: tman_display.h
shuffle_pads.o: undo.h edUtils.h
stack_dump.o: $(SRCROOT)/Misc/array.h
stack_dump.o: $(SRCROOT)/Misc/bitmap.h
stack_dump.o: $(SRCROOT)/Misc/misc.h
stack_dump.o: $(SRCROOT)/Misc/os.h
stack_dump.o: $(SRCROOT)/Misc/xalloc.h
stack_dump.o: $(SRCROOT)/Misc/xerror.h
stack_dump.o: $(SRCROOT)/g/g-error.h
stack_dump.o: $(PWD)/staden_config.h
stack_dump.o: gap-error.h
stack_dump.o: stack_dump.h
stop_codon.o: $(SRCROOT)/Misc/array.h
stop_codon.o: $(SRCROOT)/Misc/bitmap.h
stop_codon.o: $(SRCROOT)/Misc/dstring.h
stop_codon.o: $(SRCROOT)/Misc/misc.h
stop_codon.o: $(SRCROOT)/Misc/os.h
stop_codon.o: $(SRCROOT)/Misc/xalloc.h
stop_codon.o: $(SRCROOT)/Misc/xerror.h
stop_codon.o: $(SRCROOT)/g/freetree.h
stop_codon.o: $(SRCROOT)/g/g-defs.h
stop_codon.o: $(SRCROOT)/g/g-error.h
stop_codon.o: $(SRCROOT)/g/g-filedefs.h
stop_codon.o: $(SRCROOT)/g/g-os.h
stop_codon.o: $(SRCROOT)/g/g-struct.h
stop_codon.o: $(SRCROOT)/primer3/src/dpal.h
stop_codon.o: $(SRCROOT)/primer3/src/primer3.h
stop_codon.o: $(SRCROOT)/seq_utils/dna_utils.h
stop_codon.o: $(SRCROOT)/seq_utils/renz_utils.h
stop_codon.o: $(PWD)/staden_config.h
stop_codon.o: $(SRCROOT)/tk_utils/canvas_box.h
stop_codon.o: $(SRCROOT)/tk_utils/cli_arg.h
stop_codon.o: $(SRCROOT)/tk_utils/postscript.h
stop_codon.o: $(SRCROOT)/tk_utils/restriction_enzyme_map.h
stop_codon.o: $(SRCROOT)/tk_utils/sheet.h
stop_codon.o: $(SRCROOT)/tk_utils/tclXkeylist.h
stop_codon.o: $(SRCROOT)/tk_utils/tcl_utils.h
stop_codon.o: $(SRCROOT)/tk_utils/tkSheet.h
stop_codon.o: $(SRCROOT)/tk_utils/tkSheet_common.h
stop_codon.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
stop_codon.o: $(SRCROOT)/tk_utils/tkTrace.h
stop_codon.o: edStructs.h fort.h
stop_codon.o: edUtils.h
stop_codon.o: gap-dbstruct.h
stop_codon.o: gap-if.h
stop_codon.o: gap_canvas_box.h
stop_codon.o: gap_globals.h
stop_codon.o: hash.h list.h
stop_codon.o: io-reg.h IO.h
stop_codon.o: io_handle.h
stop_codon.o: io_utils.h
stop_codon.o: primlib.h
stop_codon.o: qual.h
stop_codon.o: restriction_enzymes.h
stop_codon.o: ruler_display.h
stop_codon.o: stop_codon.h
stop_codon.o: tagDefs.h fortran.h
stop_codon.o: template.h
stop_codon.o: template_display.h
stop_codon.o: tkEdNames.h
stop_codon.o: tkEditor.h
stop_codon.o: tman_display.h
stop_codon.o: undo.h tagUtils.h
strand_coverage.o: $(SRCROOT)/Misc/array.h
strand_coverage.o: $(SRCROOT)/Misc/bitmap.h
strand_coverage.o: $(SRCROOT)/Misc/misc.h
strand_coverage.o: $(SRCROOT)/Misc/os.h
strand_coverage.o: $(SRCROOT)/Misc/xalloc.h
strand_coverage.o: $(SRCROOT)/Misc/xerror.h
strand_coverage.o: $(SRCROOT)/g/freetree.h
strand_coverage.o: $(SRCROOT)/g/g-defs.h
strand_coverage.o: $(SRCROOT)/g/g-error.h
strand_coverage.o: $(SRCROOT)/g/g-filedefs.h
strand_coverage.o: $(SRCROOT)/g/g-os.h
strand_coverage.o: $(SRCROOT)/g/g-struct.h
strand_coverage.o: $(PWD)/staden_config.h
strand_coverage.o: $(SRCROOT)/tk_utils/canvas_box.h
strand_coverage.o: $(SRCROOT)/tk_utils/tclXkeylist.h
strand_coverage.o: $(SRCROOT)/tk_utils/tcl_utils.h
strand_coverage.o: $(SRCROOT)/tk_utils/text_output.h
strand_coverage.o: consen.h
strand_coverage.o: consistency_display.h
strand_coverage.o: gap-dbstruct.h
strand_coverage.o: gap-if.h
strand_coverage.o: gap_canvas_box.h
strand_coverage.o: gap_globals.h
strand_coverage.o: hash.h list.h
strand_coverage.o: io-reg.h IO.h
strand_coverage.o: io_handle.h
strand_coverage.o: io_utils.h
strand_coverage.o: newgap_cmds.h
strand_coverage.o: qual.h
strand_coverage.o: ruler_display.h
strand_coverage.o: strand_coverage.h
strand_coverage.o: template.h
strand_coverage.o: template_display.h
tagEditor.o: $(SRCROOT)/Misc/array.h
tagEditor.o: $(SRCROOT)/Misc/bitmap.h
tagEditor.o: $(SRCROOT)/Misc/dstring.h
tagEditor.o: $(SRCROOT)/Misc/misc.h
tagEditor.o: $(SRCROOT)/Misc/os.h
tagEditor.o: $(SRCROOT)/Misc/xalloc.h
tagEditor.o: $(SRCROOT)/Misc/xerror.h
tagEditor.o: $(SRCROOT)/g/freetree.h
tagEditor.o: $(SRCROOT)/g/g-defs.h
tagEditor.o: $(SRCROOT)/g/g-error.h
tagEditor.o: $(SRCROOT)/g/g-filedefs.h
tagEditor.o: $(SRCROOT)/g/g-os.h
tagEditor.o: $(SRCROOT)/g/g-struct.h
tagEditor.o: $(SRCROOT)/primer3/src/dpal.h
tagEditor.o: $(SRCROOT)/primer3/src/primer3.h
tagEditor.o: $(PWD)/staden_config.h
tagEditor.o: $(SRCROOT)/tk_utils/cli_arg.h
tagEditor.o: $(SRCROOT)/tk_utils/postscript.h
tagEditor.o: $(SRCROOT)/tk_utils/sheet.h
tagEditor.o: $(SRCROOT)/tk_utils/tcl_utils.h
tagEditor.o: $(SRCROOT)/tk_utils/tkSheet.h
tagEditor.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tagEditor.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tagEditor.o: $(SRCROOT)/tk_utils/tkTrace.h
tagEditor.o: IO.h
tagEditor.o: edStructs.h fort.h
tagEditor.o: edUtils.h qual.h
tagEditor.o: gap-dbstruct.h
tagEditor.o: gap-if.h
tagEditor.o: gap_globals.h
tagEditor.o: io_handle.h
tagEditor.o: io_utils.h
tagEditor.o: primlib.h
tagEditor.o: select.h
tagEditor.o: tagDefs.h fortran.h
tagEditor.o: template.h list.h
tagEditor.o: tkEdNames.h io-reg.h
tagEditor.o: tkEditor.h
tagEditor.o: tman_display.h
tagEditor.o: undo.h tagUtils.h
tagU1.o: $(SRCROOT)/Misc/array.h
tagU1.o: $(SRCROOT)/Misc/bitmap.h
tagU1.o: $(SRCROOT)/Misc/dstring.h
tagU1.o: $(SRCROOT)/Misc/misc.h
tagU1.o: $(SRCROOT)/Misc/os.h
tagU1.o: $(SRCROOT)/Misc/xalloc.h
tagU1.o: $(SRCROOT)/Misc/xerror.h
tagU1.o: $(SRCROOT)/g/freetree.h
tagU1.o: $(SRCROOT)/g/g-defs.h
tagU1.o: $(SRCROOT)/g/g-error.h
tagU1.o: $(SRCROOT)/g/g-filedefs.h
tagU1.o: $(SRCROOT)/g/g-os.h
tagU1.o: $(SRCROOT)/g/g-struct.h
tagU1.o: $(SRCROOT)/primer3/src/dpal.h
tagU1.o: $(SRCROOT)/primer3/src/primer3.h
tagU1.o: $(SRCROOT)/seq_utils/dna_utils.h
tagU1.o: $(PWD)/staden_config.h
tagU1.o: $(SRCROOT)/tk_utils/cli_arg.h
tagU1.o: $(SRCROOT)/tk_utils/postscript.h
tagU1.o: $(SRCROOT)/tk_utils/sheet.h
tagU1.o: $(SRCROOT)/tk_utils/tkSheet.h
tagU1.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tagU1.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tagU1.o: $(SRCROOT)/tk_utils/tkTrace.h
tagU1.o: IO.h
tagU1.o: contigEditor.h tagdb.h
tagU1.o: edStructs.h
tagU1.o: edUtils.h qual.h
tagU1.o: fort.h io_utils.h
tagU1.o: gap-dbstruct.h
tagU1.o: gap-if.h
tagU1.o: io_handle.h tkEditor.h
tagU1.o: notedb.h active_tags.h
tagU1.o: primlib.h
tagU1.o: select.h
tagU1.o: tagDefs.h undo.h
tagU1.o: tagUtils.h fortran.h
tagU1.o: template.h list.h
tagU1.o: tkEdNames.h io-reg.h
tagU1.o: tman_display.h
tagU2.o: $(SRCROOT)/Misc/FtoC.h
tagU2.o: $(SRCROOT)/Misc/array.h
tagU2.o: $(SRCROOT)/Misc/bitmap.h
tagU2.o: $(SRCROOT)/Misc/dstring.h
tagU2.o: $(SRCROOT)/Misc/misc.h
tagU2.o: $(SRCROOT)/Misc/os.h
tagU2.o: $(SRCROOT)/Misc/xalloc.h
tagU2.o: $(SRCROOT)/Misc/xerror.h
tagU2.o: $(SRCROOT)/g/freetree.h
tagU2.o: $(SRCROOT)/g/g-defs.h
tagU2.o: $(SRCROOT)/g/g-error.h
tagU2.o: $(SRCROOT)/g/g-filedefs.h
tagU2.o: $(SRCROOT)/g/g-os.h
tagU2.o: $(SRCROOT)/g/g-struct.h
tagU2.o: $(SRCROOT)/primer3/src/dpal.h
tagU2.o: $(SRCROOT)/primer3/src/primer3.h
tagU2.o: $(PWD)/staden_config.h
tagU2.o: $(SRCROOT)/tk_utils/cli_arg.h
tagU2.o: $(SRCROOT)/tk_utils/postscript.h
tagU2.o: $(SRCROOT)/tk_utils/sheet.h
tagU2.o: $(SRCROOT)/tk_utils/tcl_utils.h
tagU2.o: $(SRCROOT)/tk_utils/text_output.h
tagU2.o: $(SRCROOT)/tk_utils/tkSheet.h
tagU2.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tagU2.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tagU2.o: $(SRCROOT)/tk_utils/tkTrace.h
tagU2.o: IO.h
tagU2.o: contigEditor.h
tagU2.o: dbcheck.h active_tags.h
tagU2.o: edUtils.h qual.h
tagU2.o: fort.h io_utils.h
tagU2.o: gap-dbstruct.h
tagU2.o: gap-if.h
tagU2.o: io_handle.h tkEditor.h
tagU2.o: primlib.h
tagU2.o: qualIO.h gap_globals.h
tagU2.o: tagDefs.h fortran.h
tagU2.o: tagUtils.h edStructs.h
tagU2.o: tagdb.h
tagU2.o: template.h list.h undo.h
tagU2.o: tkEdNames.h io-reg.h
tagU2.o: tman_display.h
tagdb.o: $(SRCROOT)/Misc/misc.h
tagdb.o: $(SRCROOT)/Misc/os.h
tagdb.o: $(SRCROOT)/Misc/parse_db.h
tagdb.o: $(SRCROOT)/Misc/xalloc.h
tagdb.o: $(PWD)/staden_config.h
tagdb.o: tagdb.h
template.o: $(SRCROOT)/Misc/array.h
template.o: $(SRCROOT)/Misc/bitmap.h
template.o: $(SRCROOT)/Misc/dstring.h
template.o: $(SRCROOT)/Misc/misc.h
template.o: $(SRCROOT)/Misc/os.h
template.o: $(SRCROOT)/Misc/xalloc.h
template.o: $(SRCROOT)/Misc/xerror.h
template.o: $(SRCROOT)/g/freetree.h
template.o: $(SRCROOT)/g/g-defs.h
template.o: $(SRCROOT)/g/g-error.h
template.o: $(SRCROOT)/g/g-filedefs.h
template.o: $(SRCROOT)/g/g-os.h
template.o: $(SRCROOT)/g/g-struct.h
template.o: $(SRCROOT)/primer3/src/dpal.h
template.o: $(SRCROOT)/primer3/src/primer3.h
template.o: $(PWD)/staden_config.h
template.o: $(SRCROOT)/tk_utils/cli_arg.h
template.o: $(SRCROOT)/tk_utils/postscript.h
template.o: $(SRCROOT)/tk_utils/sheet.h
template.o: $(SRCROOT)/tk_utils/tcl_utils.h
template.o: $(SRCROOT)/tk_utils/tkSheet.h
template.o: $(SRCROOT)/tk_utils/tkSheet_common.h
template.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
template.o: $(SRCROOT)/tk_utils/tkTrace.h
template.o: IO.h
template.o: edStructs.h fort.h
template.o: gap-dbstruct.h
template.o: gap-if.h
template.o: gap_globals.h
template.o: io_handle.h
template.o: io_utils.h
template.o: primlib.h
template.o: qual.h tman_display.h
template.o: tagDefs.h fortran.h
template.o: tagUtils.h
template.o: template.h list.h
template.o: tkEdNames.h io-reg.h
template.o: tkEditor.h
template.o: undo.h edUtils.h
template_display.o: $(SRCROOT)/Misc/FtoC.h
template_display.o: $(SRCROOT)/Misc/array.h
template_display.o: $(SRCROOT)/Misc/bitmap.h
template_display.o: $(SRCROOT)/Misc/dstring.h
template_display.o: $(SRCROOT)/Misc/misc.h
template_display.o: $(SRCROOT)/Misc/os.h
template_display.o: $(SRCROOT)/Misc/vlen.h
template_display.o: $(SRCROOT)/Misc/xalloc.h
template_display.o: $(SRCROOT)/Misc/xerror.h
template_display.o: $(SRCROOT)/g/freetree.h
template_display.o: $(SRCROOT)/g/g-defs.h
template_display.o: $(SRCROOT)/g/g-error.h
template_display.o: $(SRCROOT)/g/g-filedefs.h
template_display.o: $(SRCROOT)/g/g-os.h
template_display.o: $(SRCROOT)/g/g-struct.h
template_display.o: $(SRCROOT)/primer3/src/dpal.h
template_display.o: $(SRCROOT)/primer3/src/primer3.h
template_display.o: $(PWD)/staden_config.h
template_display.o: $(SRCROOT)/tk_utils/canvas_box.h
template_display.o: $(SRCROOT)/tk_utils/cli_arg.h
template_display.o: $(SRCROOT)/tk_utils/postscript.h
template_display.o: $(SRCROOT)/tk_utils/ruler_tick.h
template_display.o: $(SRCROOT)/tk_utils/sheet.h
template_display.o: $(SRCROOT)/tk_utils/tclXkeylist.h
template_display.o: $(SRCROOT)/tk_utils/tcl_utils.h
template_display.o: $(SRCROOT)/tk_utils/text_output.h
template_display.o: $(SRCROOT)/tk_utils/tkSheet.h
template_display.o: $(SRCROOT)/tk_utils/tkSheet_common.h
template_display.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
template_display.o: $(SRCROOT)/tk_utils/tkTrace.h
template_display.o: active_tags.h
template_display.o: contig_selector.h
template_display.o: cs-object.h
template_display.o: edStructs.h
template_display.o: edUtils.h
template_display.o: fort.h
template_display.o: fortran.h
template_display.o: gap-dbstruct.h
template_display.o: gap-if.h
template_display.o: gap_canvas_box.h
template_display.o: gap_globals.h
template_display.o: hash.h list.h
template_display.o: io-reg.h IO.h
template_display.o: io_handle.h
template_display.o: io_utils.h
template_display.o: primlib.h
template_display.o: qual.h
template_display.o: ruler_display.h
template_display.o: tagDefs.h
template_display.o: tagUtils.h
template_display.o: tagdb.h
template_display.o: template.h
template_display.o: template_display.h
template_display.o: tkEdNames.h
template_display.o: tkEditor.h
template_display.o: tman_display.h
template_display.o: undo.h
text-io-reg.o: $(SRCROOT)/Misc/array.h
text-io-reg.o: $(SRCROOT)/Misc/bitmap.h
text-io-reg.o: $(SRCROOT)/Misc/os.h
text-io-reg.o: $(SRCROOT)/Misc/xerror.h
text-io-reg.o: $(SRCROOT)/g/freetree.h
text-io-reg.o: $(SRCROOT)/g/g-defs.h
text-io-reg.o: $(SRCROOT)/g/g-error.h
text-io-reg.o: $(SRCROOT)/g/g-filedefs.h
text-io-reg.o: $(SRCROOT)/g/g-os.h
text-io-reg.o: $(SRCROOT)/g/g-struct.h
text-io-reg.o: $(PWD)/staden_config.h
text-io-reg.o: IO.h
text-io-reg.o: gap-dbstruct.h
text-io-reg.o: gap-if.h
text-io-reg.o: io_handle.h
text-io-reg.o: io_utils.h
tk-io-reg.o: $(SRCROOT)/Misc/array.h
tk-io-reg.o: $(SRCROOT)/Misc/bitmap.h
tk-io-reg.o: $(SRCROOT)/Misc/misc.h
tk-io-reg.o: $(SRCROOT)/Misc/os.h
tk-io-reg.o: $(SRCROOT)/Misc/xalloc.h
tk-io-reg.o: $(SRCROOT)/Misc/xerror.h
tk-io-reg.o: $(SRCROOT)/g/freetree.h
tk-io-reg.o: $(SRCROOT)/g/g-defs.h
tk-io-reg.o: $(SRCROOT)/g/g-error.h
tk-io-reg.o: $(SRCROOT)/g/g-filedefs.h
tk-io-reg.o: $(SRCROOT)/g/g-os.h
tk-io-reg.o: $(SRCROOT)/g/g-struct.h
tk-io-reg.o: $(PWD)/staden_config.h
tk-io-reg.o: $(SRCROOT)/tk_utils/canvas_box.h
tk-io-reg.o: $(SRCROOT)/tk_utils/cli_arg.h
tk-io-reg.o: $(SRCROOT)/tk_utils/tclXkeylist.h
tk-io-reg.o: $(SRCROOT)/tk_utils/tcl_utils.h
tk-io-reg.o: IO.h
tk-io-reg.o: consistency_display.h
tk-io-reg.o: contig_selector.h
tk-io-reg.o: cs-object.h
tk-io-reg.o: gap-dbstruct.h
tk-io-reg.o: gap-if.h
tk-io-reg.o: gap_cli_arg.h
tk-io-reg.o: gap_globals.h
tk-io-reg.o: hash.h template.h
tk-io-reg.o: io_handle.h
tk-io-reg.o: io_utils.h
tk-io-reg.o: newgap_cmds.h list.h
tk-io-reg.o: newgap_structs.h
tk-io-reg.o: ruler_display.h
tk-io-reg.o: template_display.h
tk-io-reg.o: tk-io-reg.h io-reg.h
tkEdNames.o: $(SRCROOT)/Misc/array.h
tkEdNames.o: $(SRCROOT)/Misc/bitmap.h
tkEdNames.o: $(SRCROOT)/Misc/dstring.h
tkEdNames.o: $(SRCROOT)/Misc/os.h
tkEdNames.o: $(SRCROOT)/Misc/xalloc.h
tkEdNames.o: $(SRCROOT)/Misc/xerror.h
tkEdNames.o: $(SRCROOT)/g/freetree.h
tkEdNames.o: $(SRCROOT)/g/g-defs.h
tkEdNames.o: $(SRCROOT)/g/g-error.h
tkEdNames.o: $(SRCROOT)/g/g-filedefs.h
tkEdNames.o: $(SRCROOT)/g/g-os.h
tkEdNames.o: $(SRCROOT)/g/g-struct.h
tkEdNames.o: $(SRCROOT)/primer3/src/dpal.h
tkEdNames.o: $(SRCROOT)/primer3/src/primer3.h
tkEdNames.o: $(PWD)/staden_config.h
tkEdNames.o: $(SRCROOT)/tk_utils/cli_arg.h
tkEdNames.o: $(SRCROOT)/tk_utils/postscript.h
tkEdNames.o: $(SRCROOT)/tk_utils/sheet.h
tkEdNames.o: $(SRCROOT)/tk_utils/tcl_utils.h
tkEdNames.o: $(SRCROOT)/tk_utils/tkSheet.h
tkEdNames.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tkEdNames.o: $(SRCROOT)/tk_utils/tkSheet_config.h
tkEdNames.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tkEdNames.o: $(SRCROOT)/tk_utils/tkTrace.h
tkEdNames.o: $(SRCROOT)/tk_utils/tk_defs.h
tkEdNames.o: IO.h
tkEdNames.o: edStructs.h
tkEdNames.o: edUtils.h qual.h
tkEdNames.o: fort.h io_utils.h
tkEdNames.o: gap-dbstruct.h
tkEdNames.o: gap-if.h
tkEdNames.o: io_handle.h
tkEdNames.o: notes.h
tkEdNames.o: primlib.h
tkEdNames.o: tagDefs.h fortran.h
tkEdNames.o: template.h list.h
tkEdNames.o: tkEdNames.h
tkEdNames.o: tkEditor.h io-reg.h
tkEdNames.o: tman_display.h
tkEdNames.o: undo.h tagUtils.h
tkEdUtils.o: $(SRCROOT)/Misc/array.h
tkEdUtils.o: $(SRCROOT)/Misc/bitmap.h
tkEdUtils.o: $(SRCROOT)/Misc/dstring.h
tkEdUtils.o: $(SRCROOT)/Misc/misc.h
tkEdUtils.o: $(SRCROOT)/Misc/os.h
tkEdUtils.o: $(SRCROOT)/Misc/xalloc.h
tkEdUtils.o: $(SRCROOT)/Misc/xerror.h
tkEdUtils.o: $(SRCROOT)/g/freetree.h
tkEdUtils.o: $(SRCROOT)/g/g-defs.h
tkEdUtils.o: $(SRCROOT)/g/g-error.h
tkEdUtils.o: $(SRCROOT)/g/g-filedefs.h
tkEdUtils.o: $(SRCROOT)/g/g-os.h
tkEdUtils.o: $(SRCROOT)/g/g-struct.h
tkEdUtils.o: $(SRCROOT)/primer3/src/dpal.h
tkEdUtils.o: $(SRCROOT)/primer3/src/primer3.h
tkEdUtils.o: $(SRCROOT)/seq_utils/dna_utils.h
tkEdUtils.o: $(SRCROOT)/seq_utils/genetic_code.h
tkEdUtils.o: $(PWD)/staden_config.h
tkEdUtils.o: $(SRCROOT)/tk_utils/cli_arg.h
tkEdUtils.o: $(SRCROOT)/tk_utils/postscript.h
tkEdUtils.o: $(SRCROOT)/tk_utils/sheet.h
tkEdUtils.o: $(SRCROOT)/tk_utils/tcl_utils.h
tkEdUtils.o: $(SRCROOT)/tk_utils/tkSheet.h
tkEdUtils.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tkEdUtils.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tkEdUtils.o: $(SRCROOT)/tk_utils/tkTrace.h
tkEdUtils.o: IO.h
tkEdUtils.o: edStructs.h
tkEdUtils.o: edUtils.h qual.h
tkEdUtils.o: fort.h io_utils.h
tkEdUtils.o: gap-dbstruct.h
tkEdUtils.o: gap-if.h
tkEdUtils.o: gap_globals.h
tkEdUtils.o: io_handle.h
tkEdUtils.o: notedb.h
tkEdUtils.o: primlib.h
tkEdUtils.o: tagDefs.h fortran.h
tkEdUtils.o: template.h list.h
tkEdUtils.o: tkEdNames.h io-reg.h
tkEdUtils.o: tkEditor.h
tkEdUtils.o: tman_display.h
tkEdUtils.o: tman_interface.h
tkEdUtils.o: undo.h tagUtils.h
tkEditor.o: $(SRCROOT)/Misc/array.h
tkEditor.o: $(SRCROOT)/Misc/bitmap.h
tkEditor.o: $(SRCROOT)/Misc/dstring.h
tkEditor.o: $(SRCROOT)/Misc/os.h
tkEditor.o: $(SRCROOT)/Misc/xalloc.h
tkEditor.o: $(SRCROOT)/Misc/xerror.h
tkEditor.o: $(SRCROOT)/g/freetree.h
tkEditor.o: $(SRCROOT)/g/g-defs.h
tkEditor.o: $(SRCROOT)/g/g-error.h
tkEditor.o: $(SRCROOT)/g/g-filedefs.h
tkEditor.o: $(SRCROOT)/g/g-os.h
tkEditor.o: $(SRCROOT)/g/g-struct.h
tkEditor.o: $(SRCROOT)/primer3/src/dpal.h
tkEditor.o: $(SRCROOT)/primer3/src/primer3.h
tkEditor.o: $(PWD)/staden_config.h
tkEditor.o: $(SRCROOT)/tk_utils/cli_arg.h
tkEditor.o: $(SRCROOT)/tk_utils/postscript.h
tkEditor.o: $(SRCROOT)/tk_utils/sheet.h
tkEditor.o: $(SRCROOT)/tk_utils/tcl_utils.h
tkEditor.o: $(SRCROOT)/tk_utils/tkSheet.h
tkEditor.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tkEditor.o: $(SRCROOT)/tk_utils/tkSheet_config.h
tkEditor.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tkEditor.o: $(SRCROOT)/tk_utils/tkTrace.h
tkEditor.o: $(SRCROOT)/tk_utils/tk_defs.h
tkEditor.o: IO.h
tkEditor.o: contigEditor.h
tkEditor.o: edCommands.h
tkEditor.o: edStructs.h
tkEditor.o: edUtils.h qual.h
tkEditor.o: fort.h io_utils.h
tkEditor.o: gap-dbstruct.h
tkEditor.o: gap-if.h
tkEditor.o: gap_globals.h
tkEditor.o: io_handle.h
tkEditor.o: oligo.h
tkEditor.o: primlib.h
tkEditor.o: tagDefs.h fortran.h
tkEditor.o: template.h list.h
tkEditor.o: tkEdNames.h io-reg.h
tkEditor.o: tkEditor.h
tkEditor.o: tman_display.h
tkEditor.o: tman_interface.h
tkEditor.o: undo.h tagUtils.h
tman_cons.o: $(SRCROOT)/Misc/array.h
tman_cons.o: $(SRCROOT)/Misc/bitmap.h
tman_cons.o: $(SRCROOT)/Misc/dstring.h
tman_cons.o: $(SRCROOT)/Misc/misc.h
tman_cons.o: $(SRCROOT)/Misc/os.h
tman_cons.o: $(SRCROOT)/Misc/xalloc.h
tman_cons.o: $(SRCROOT)/Misc/xerror.h
tman_cons.o: $(SRCROOT)/g/freetree.h
tman_cons.o: $(SRCROOT)/g/g-defs.h
tman_cons.o: $(SRCROOT)/g/g-error.h
tman_cons.o: $(SRCROOT)/g/g-filedefs.h
tman_cons.o: $(SRCROOT)/g/g-os.h
tman_cons.o: $(SRCROOT)/g/g-struct.h
tman_cons.o: $(SRCROOT)/primer3/src/dpal.h
tman_cons.o: $(SRCROOT)/primer3/src/primer3.h
tman_cons.o: $(PWD)/staden_config.h
tman_cons.o: $(SRCROOT)/tk_utils/cli_arg.h
tman_cons.o: $(SRCROOT)/tk_utils/postscript.h
tman_cons.o: $(SRCROOT)/tk_utils/sheet.h
tman_cons.o: $(SRCROOT)/tk_utils/tclXkeylist.h
tman_cons.o: $(SRCROOT)/tk_utils/tcl_utils.h
tman_cons.o: $(SRCROOT)/tk_utils/tkSheet.h
tman_cons.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tman_cons.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tman_cons.o: $(SRCROOT)/tk_utils/tkTrace.h
tman_cons.o: $(SRCROOT)/tk_utils/tkTraceIO.h
tman_cons.o: IO.h
tman_cons.o: edStructs.h
tman_cons.o: edUtils.h qual.h
tman_cons.o: fort.h io_utils.h
tman_cons.o: gap-dbstruct.h
tman_cons.o: gap-if.h
tman_cons.o: gap_globals.h
tman_cons.o: io_handle.h
tman_cons.o: primlib.h
tman_cons.o: tagDefs.h fortran.h
tman_cons.o: template.h list.h
tman_cons.o: tkEdNames.h io-reg.h
tman_cons.o: tkEditor.h
tman_cons.o: tman_display.h
tman_cons.o: tman_interface.h
tman_cons.o: undo.h tagUtils.h
tman_diff.o: $(SRCROOT)/Misc/array.h
tman_diff.o: $(SRCROOT)/Misc/bitmap.h
tman_diff.o: $(SRCROOT)/Misc/dstring.h
tman_diff.o: $(SRCROOT)/Misc/misc.h
tman_diff.o: $(SRCROOT)/Misc/os.h
tman_diff.o: $(SRCROOT)/Misc/xalloc.h
tman_diff.o: $(SRCROOT)/Misc/xerror.h
tman_diff.o: $(SRCROOT)/g/freetree.h
tman_diff.o: $(SRCROOT)/g/g-defs.h
tman_diff.o: $(SRCROOT)/g/g-error.h
tman_diff.o: $(SRCROOT)/g/g-filedefs.h
tman_diff.o: $(SRCROOT)/g/g-os.h
tman_diff.o: $(SRCROOT)/g/g-struct.h
tman_diff.o: $(SRCROOT)/mutlib/mutlib.h
tman_diff.o: $(SRCROOT)/mutlib/staden.h
tman_diff.o: $(SRCROOT)/primer3/src/dpal.h
tman_diff.o: $(SRCROOT)/primer3/src/primer3.h
tman_diff.o: $(SRCROOT)/seq_utils/dna_utils.h
tman_diff.o: $(PWD)/staden_config.h
tman_diff.o: $(SRCROOT)/tk_utils/cli_arg.h
tman_diff.o: $(SRCROOT)/tk_utils/postscript.h
tman_diff.o: $(SRCROOT)/tk_utils/sheet.h
tman_diff.o: $(SRCROOT)/tk_utils/tclXkeylist.h
tman_diff.o: $(SRCROOT)/tk_utils/tcl_utils.h
tman_diff.o: $(SRCROOT)/tk_utils/tkSheet.h
tman_diff.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tman_diff.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tman_diff.o: $(SRCROOT)/tk_utils/tkTrace.h
tman_diff.o: $(SRCROOT)/tk_utils/tkTraceIO.h
tman_diff.o: IO.h
tman_diff.o: edStructs.h
tman_diff.o: edUtils.h qual.h
tman_diff.o: fort.h io_utils.h
tman_diff.o: gap-dbstruct.h
tman_diff.o: gap-if.h
tman_diff.o: gap_globals.h
tman_diff.o: io_handle.h
tman_diff.o: primlib.h
tman_diff.o: tagDefs.h fortran.h
tman_diff.o: template.h list.h
tman_diff.o: tkEdNames.h io-reg.h
tman_diff.o: tkEditor.h
tman_diff.o: tman_display.h
tman_diff.o: tman_interface.h
tman_diff.o: undo.h tagUtils.h
tman_display.o: $(SRCROOT)/Misc/array.h
tman_display.o: $(SRCROOT)/Misc/bitmap.h
tman_display.o: $(SRCROOT)/Misc/dstring.h
tman_display.o: $(SRCROOT)/Misc/misc.h
tman_display.o: $(SRCROOT)/Misc/os.h
tman_display.o: $(SRCROOT)/Misc/xalloc.h
tman_display.o: $(SRCROOT)/Misc/xerror.h
tman_display.o: $(SRCROOT)/g/freetree.h
tman_display.o: $(SRCROOT)/g/g-defs.h
tman_display.o: $(SRCROOT)/g/g-error.h
tman_display.o: $(SRCROOT)/g/g-filedefs.h
tman_display.o: $(SRCROOT)/g/g-os.h
tman_display.o: $(SRCROOT)/g/g-struct.h
tman_display.o: $(SRCROOT)/primer3/src/dpal.h
tman_display.o: $(SRCROOT)/primer3/src/primer3.h
tman_display.o: $(PWD)/staden_config.h
tman_display.o: $(SRCROOT)/tk_utils/cli_arg.h
tman_display.o: $(SRCROOT)/tk_utils/postscript.h
tman_display.o: $(SRCROOT)/tk_utils/sheet.h
tman_display.o: $(SRCROOT)/tk_utils/tclXkeylist.h
tman_display.o: $(SRCROOT)/tk_utils/tcl_utils.h
tman_display.o: $(SRCROOT)/tk_utils/tkSheet.h
tman_display.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tman_display.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tman_display.o: $(SRCROOT)/tk_utils/tkTrace.h
tman_display.o: IO.h
tman_display.o: edStructs.h
tman_display.o: fort.h io_utils.h
tman_display.o: fortran.h
tman_display.o: gap-dbstruct.h
tman_display.o: gap-if.h
tman_display.o: gap_globals.h
tman_display.o: io-reg.h
tman_display.o: io_handle.h
tman_display.o: primlib.h
tman_display.o: qual.h tagUtils.h
tman_display.o: tagDefs.h
tman_display.o: template.h list.h
tman_display.o: tkEdNames.h
tman_display.o: tkEditor.h
tman_display.o: tman_display.h
tman_display.o: tman_interface.h
tman_display.o: undo.h edUtils.h
tman_interface.o: $(SRCROOT)/Misc/array.h
tman_interface.o: $(SRCROOT)/Misc/bitmap.h
tman_interface.o: $(SRCROOT)/Misc/dstring.h
tman_interface.o: $(SRCROOT)/Misc/misc.h
tman_interface.o: $(SRCROOT)/Misc/os.h
tman_interface.o: $(SRCROOT)/Misc/xalloc.h
tman_interface.o: $(SRCROOT)/Misc/xerror.h
tman_interface.o: $(SRCROOT)/g/freetree.h
tman_interface.o: $(SRCROOT)/g/g-defs.h
tman_interface.o: $(SRCROOT)/g/g-error.h
tman_interface.o: $(SRCROOT)/g/g-filedefs.h
tman_interface.o: $(SRCROOT)/g/g-os.h
tman_interface.o: $(SRCROOT)/g/g-struct.h
tman_interface.o: $(SRCROOT)/primer3/src/dpal.h
tman_interface.o: $(SRCROOT)/primer3/src/primer3.h
tman_interface.o: $(PWD)/staden_config.h
tman_interface.o: $(SRCROOT)/tk_utils/cli_arg.h
tman_interface.o: $(SRCROOT)/tk_utils/postscript.h
tman_interface.o: $(SRCROOT)/tk_utils/sheet.h
tman_interface.o: $(SRCROOT)/tk_utils/tkSheet.h
tman_interface.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tman_interface.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tman_interface.o: $(SRCROOT)/tk_utils/tkTrace.h
tman_interface.o: IO.h
tman_interface.o: contigEditor.h
tman_interface.o: edStructs.h
tman_interface.o: edUtils.h
tman_interface.o: fort.h
tman_interface.o: fortran.h
tman_interface.o: gap-dbstruct.h
tman_interface.o: gap-if.h
tman_interface.o: io-reg.h
tman_interface.o: io_handle.h
tman_interface.o: io_utils.h
tman_interface.o: list.h undo.h
tman_interface.o: primlib.h
tman_interface.o: qual.h
tman_interface.o: tagDefs.h
tman_interface.o: tagUtils.h
tman_interface.o: template.h
tman_interface.o: tkEdNames.h
tman_interface.o: tkEditor.h
tman_interface.o: tman_display.h
tman_interface.o: tman_interface.h
undo.o: $(SRCROOT)/Misc/array.h
undo.o: $(SRCROOT)/Misc/bitmap.h
undo.o: $(SRCROOT)/Misc/dstring.h
undo.o: $(SRCROOT)/Misc/os.h
undo.o: $(SRCROOT)/Misc/xalloc.h
undo.o: $(SRCROOT)/Misc/xerror.h
undo.o: $(SRCROOT)/g/freetree.h
undo.o: $(SRCROOT)/g/g-defs.h
undo.o: $(SRCROOT)/g/g-error.h
undo.o: $(SRCROOT)/g/g-filedefs.h
undo.o: $(SRCROOT)/g/g-os.h
undo.o: $(SRCROOT)/g/g-struct.h
undo.o: $(SRCROOT)/primer3/src/dpal.h
undo.o: $(SRCROOT)/primer3/src/primer3.h
undo.o: $(PWD)/staden_config.h
undo.o: $(SRCROOT)/tk_utils/cli_arg.h
undo.o: $(SRCROOT)/tk_utils/postscript.h
undo.o: $(SRCROOT)/tk_utils/sheet.h
undo.o: $(SRCROOT)/tk_utils/text_output.h
undo.o: $(SRCROOT)/tk_utils/tkSheet.h
undo.o: $(SRCROOT)/tk_utils/tkSheet_common.h
undo.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
undo.o: $(SRCROOT)/tk_utils/tkTrace.h
undo.o: IO.h
undo.o: edStructs.h
undo.o: extend.h contigEditor.h
undo.o: fort.h io_utils.h
undo.o: fortran.h
undo.o: gap-dbstruct.h
undo.o: gap-if.h
undo.o: io_handle.h tkEditor.h
undo.o: primlib.h
undo.o: tagUtils.h tagDefs.h
undo.o: template.h list.h
undo.o: tkEdNames.h io-reg.h
undo.o: tman_display.h
undo.o: undo.h edUtils.h qual.h
vseqs.o: $(SRCROOT)/Misc/array.h
vseqs.o: $(SRCROOT)/Misc/bitmap.h
vseqs.o: $(SRCROOT)/Misc/misc.h
vseqs.o: $(SRCROOT)/Misc/os.h
vseqs.o: $(SRCROOT)/Misc/xalloc.h
vseqs.o: $(SRCROOT)/Misc/xerror.h
vseqs.o: $(SRCROOT)/g/freetree.h
vseqs.o: $(SRCROOT)/g/g-defs.h
vseqs.o: $(SRCROOT)/g/g-error.h
vseqs.o: $(SRCROOT)/g/g-filedefs.h
vseqs.o: $(SRCROOT)/g/g-os.h
vseqs.o: $(SRCROOT)/g/g-struct.h
vseqs.o: $(PWD)/staden_config.h
vseqs.o: $(SRCROOT)/tk_utils/text_output.h
vseqs.o: IO.h
vseqs.o: gap-dbstruct.h
vseqs.o: gap-if.h
vseqs.o: io_handle.h io_utils.h
vseqs.o: vseqs.h qual.h
