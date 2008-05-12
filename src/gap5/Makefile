LIBS = gap
ifeq ($(MACHINE),windows)
PROGS = $(LIBS) $(O)/copy_db 
else
PROGS = $(LIBS) $(O)/gap4sh $(O)/copy_db 
endif
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
#SRCROOT=/home5/pubseq/share/src
#SRCROOT=/nfs/arran/home5/pubseq/share/src
#REMOTESRC=$(GAP4SRC)

include	 $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

GAP4SRC=$(SRCROOT)/gap5

# Local comment: Comment out next line for remote compilation
GAP4BIN=$(O)

G_INC=-I/nfs/team71/psg/jkb/work/tgap
G_DEP=-L/nfs/team71/psg/jkb/work/tgap -ltgap
INCLUDES_E += -I$(HOME)/work/tgap \
	      -I$(GAP4SRC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(MUT_INC) \
	      $(F77_INC) $(P3_INC)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)


TGAP_LIB=-L/nfs/team71/psg/jkb/work/tgap -ltgap


#
# Building the programs
# This should be just a linking phase because all of the object
# files and library files are generated using implicit rules.
# We use the fortran compiler to do linking.
#
#	$(GAP4BIN)/contig_editor.o
CEDITOR = \
	$(GAP4BIN)/editor_view.o \
	$(GAP4BIN)/tkEdNames.o\
	$(GAP4BIN)/tkEditor.o\
	$(GAP4BIN)/tman_interface.o\
	$(GAP4BIN)/tman_display.o

CEDITOR_SCRAP=\
	$(GAP4BIN)/contigEditor.o\
	$(GAP4BIN)/tkEditor.o\
	$(GAP4BIN)/tkEdUtils.o\
	$(GAP4BIN)/edInterface.o\
	$(GAP4BIN)/edUtils2.o \
	$(GAP4BIN)/tagU1.o\
	$(GAP4BIN)/undo.o\
	$(GAP4BIN)/edExtend.o\
	$(GAP4BIN)/edCommands.o\
	$(GAP4BIN)/edMutations.o\
	$(GAP4BIN)/tagEditor.o\
	$(GAP4BIN)/searchUtils.o\
	$(GAP4BIN)/tman_interface.o\
	$(GAP4BIN)/tman_display.o\
	$(GAP4BIN)/tman_cons.o\
	$(GAP4BIN)/tman_diff.o\
	$(GAP4BIN)/join.o\
	$(GAP4BIN)/oligo.o

GAP5=\
	$(GAP4BIN)/newgap5_cmds.o \
	$(GAP4BIN)/gap4_compat.o \
        $(GAP4BIN)/actf.o \
	$(GAP4BIN)/gap_hash.o\
	$(GAP4BIN)/hash_lib.o \
	$(GAP4BIN)/qual.o \
	$(GAP4BIN)/qualIO.o \
	$(GAP4BIN)/gap_globals.o \
	$(GAP4BIN)/tagdb.o\
	$(GAP4BIN)/init.o \
	$(GAP4BIN)/notedb.o\
	$(GAP4BIN)/active_tags.o\
	$(GAP4BIN)/list_proc.o\
	$(GAP4BIN)/gap-error.o\
	$(GAP4BIN)/stack_dump.o\
	$(GAP4BIN)/tk-io-reg.o\
	$(GAP4BIN)/tkAppInit.o \
	$(GAP4BIN)/consen.o\
	$(GAP4BIN)/contig_selector.o\
	$(GAP4BIN)/cs-object.o\
	$(GAP4BIN)/find_repeats.o\
	$(CEDITOR)

GAP5SH_LIBS=\
	$(IOLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(MISC_LIB) \
	$(P3_LIB) \
	$(TGAP_LIB)


$(O)/gap5sh:	$(GAPOBJS) $(GAP5)
	$(CXX) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(GAP5) $(GAP5SH_LIBS) $(CXX_DEP) $(F77_DEP)



#
# Trace manager objects
#

#
# Here we define the necessary object files needed for an external program
# to link with the gap4 I/O routines
#
#
# The object file list are in  $(GAPDB_EXT_OBJS).
# The linker libraries are in  $(GAPDB_EXT_LIBS).
# The link dependencies are in $(GAPDB_EXT_DEPS).
# To use these, simply add "include $(SRCROOT)/mk/gap4_defs.mk" to your
# Makefile after the global and machine includes.
#
# For an example of usage, see the "convert" program Makefile.

GAPDB_LOW=\
	$(GAP4BIN)/gap-local.o\
	$(GAP4BIN)/gap-remote.o\
	$(GAP4BIN)/gap-if.o\
	$(GAP4BIN)/gap-init.o\
	$(GAP4BIN)/gap-dbstruct.o\
	$(GAP4BIN)/gap-create.o\
	$(GAP4BIN)/gap-error.o\
	$(GAP4BIN)/stack_dump.o\
	$(GAP4BIN)/gap-io.o

GAPDB_MID=\
        $(GAP4BIN)/io-reg.o \
        $(GAP4BIN)/actf.o

#        $(GAP4BIN)/IO.o \

GAPDB_UTILS=\
	$(GAP4BIN)/io_handle.o \
	$(GAP4BIN)/io_utils.o

# GAPDB_EXT_OBJS is basically the low and mid level files plus the
# text-io-reg.o object. This is simply stub routines to enable safe linking
# and isn't needed by gap4 itself.

GAPDB_EXT_OBJS=\
	$(GAPDB_LOW) \
	$(GAPDB_MID) \
	$(GAPDB_UTILS) \
	$(GAP4BIN)/text-io-reg.o

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
	$(TCL_LIB) \
	$(TGAP_LIB)

GAPDB=\
	$(GAPDB_LOW)\
	$(GAPDB_MID)\
	$(GAP4BIN)/seqInfo.o\
	$(GAP4BIN)/parse_ft.o\
	$(GAP4BIN)/tk-io-reg.o\
	$(GAP4BIN)/io_handle.o \
	$(GAP4BIN)/io_utils.o \
	$(GAP4BIN)/gap4_compat.o

#	$(GAP4BIN)/gap-tcl.o \
#	$(GAP4BIN)/IO2.o\
#	$(GAP4BIN)/IO3.o \

TCLBITS=\
	$(GAP4BIN)/gap_cli_arg.o \
	$(GAP4BIN)/gap_globals.o \
	$(GAP4BIN)/newgap_cmds.o \
	$(GAP4BIN)/init.o

GAPFUNCS = \
	$(GAP4BIN)/bubbl3.o\
	$(GAP4BIN)/tagdb.o\
	$(GAP4BIN)/notedb.o\
	$(GAP4BIN)/active_tags.o\
	$(GAP4BIN)/dbcheck.o\
	$(GAP4BIN)/clones.o\
	$(GAP4BIN)/extract.o\
	$(GAP4BIN)/preass.o \
	$(GAP4BIN)/list.o \
	$(GAP4BIN)/reactions.o \
	$(GAP4BIN)/probe.o \
	$(GAP4BIN)/template.o \
	$(GAP4BIN)/template_display.o \
	$(GAP4BIN)/ruler_display.o \
	$(GAP4BIN)/gap_canvas_box.o \
	$(GAP4BIN)/hash.o \
	$(GAP4BIN)/gap_array.o \
	$(GAP4BIN)/show_relationships.o \
	$(GAP4BIN)/fij.o \
	$(GAP4BIN)/do_fij.o \
	$(GAP4BIN)/auto_assemble.o \
	$(GAP4BIN)/dis_readings.o \
	$(GAP4BIN)/break_contig.o \
	$(GAP4BIN)/quality_plot.o \
	$(GAP4BIN)/readpair.o \
	$(GAP4BIN)/contig_selector.o \
	$(GAP4BIN)/complement.o \
	$(GAP4BIN)/list_proc.o\
	$(GAP4BIN)/dstrand.o\
	$(GAP4BIN)/oligo_sel.o\
	$(GAP4BIN)/primlib.o\
	$(GAP4BIN)/alter_rel.o\
	$(GAP4BIN)/restriction_enzymes.o \
	$(GAP4BIN)/stop_codon.o \
	$(GAP4BIN)/assemble_direct.o\
	$(GAP4BIN)/check_assembly.o\
	$(GAP4BIN)/tagU2.o\
	$(GAP4BIN)/mess.o\
	$(GAP4BIN)/find_oligo.o\
	$(GAP4BIN)/copy_db.o \
	$(GAP4BIN)/contig_order.o\
	$(GAP4BIN)/clip.o\
	$(GAP4BIN)/notes.o\
	$(GAP4BIN)/consistency_display.o\
	$(GAP4BIN)/consistency_canvas_box.o\
	$(GAP4BIN)/confidence_graph.o\
	$(GAP4BIN)/reading_coverage.o\
	$(GAP4BIN)/readpair_coverage.o\
	$(GAP4BIN)/strand_coverage.o\
	$(GAP4BIN)/find_fragments.o\
	$(GAP4BIN)/vseqs.o\
	$(GAP4BIN)/$(GAP4_LEGACY)\
	$(GAP4BIN)/f2c.o \
	$(GAP4BIN)/shuffle_pads.o \
	$(GAP4BIN)/auto_break.o

GAPOBJS=\
	$(COMMONOBJS)\
	$(GAPDB)\
	$(TCLBITS)\
	$(CONSEN)\
	$(GAPFUNCS)\
	$(CEDITOR)

GAPOBJS2=\
	$(GAP4BIN)/tkAppInit.o

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
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(MISC_LIB) \
	$(MUT_LIB) \
	$(P3_LIB) \
	$(TGAP_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(GAPOBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LINK_PATHFLAG)$(L) $(GAPOBJS) $(GAP4_LIBS) $(CXX_DEP) $(SHLIB_DEP) $(F77_DEP)

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(GAPOBJS)
	$(MKDEFL) $@ $(GAPOBJS)

$(O)/gap4sh:	$(GAPOBJS) $(GAPOBJS2)
	$(CXX) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(GAPOBJS) $(GAPOBJS2) $(GAP4SH_LIBS) $(CXX_DEP) $(F77_DEP)

COPYDBOBJS=\
        $(GAPDB_LOW) \
        $(GAPDB_MID) \
        $(GAP4BIN)/copy_db.o \
        $(GAP4BIN)/copy_db_main.o \
        $(GAP4BIN)/text-io-reg.o\
	$(GAP4BIN)/io_utils.o \
	$(GAP4BIN)/io_handle.o

CPLIB=\
	$(G_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB) \
	$(TCL_LIB)

$(O)/copy_db:   $(COPYDBOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(COPYDBOBJS) $(CPLIB) $(LIBSC)

THRASHOBJS=\
        $(GAPDB_LOW) \
        $(GAPDB_MID) \
        $(GAP4BIN)/gap-thrash3.o \
        $(GAP4BIN)/text-io-reg.o

THRASHBOBJS=\
        $(GAPDB_LOW) \
        $(GAPDB_MID) \
        $(GAP4BIN)/gap-thrash2bug.o \
        $(GAP4BIN)/text-io-reg.o

TLIB=\
	$(G_LIB) \
	$(TEXTUTILS_LIB) \
	./libmisc.a

$(O)/thrash:   $(THRASHOBJS)
	$(CLD) -o $@ $(THRASHOBJS) $(TLIB) $(LIBSC) -ldmalloc
$(O)/thrashb:   $(THRASHBOBJS)
	$(CLD) -o $@ $(THRASHBOBJS) $(TLIB) $(LIBSC) -ldmalloc

DEPEND_OBJ = $(GAPOBJS) $(COPYDBOBJS)

install:
	cp $(O)/copy_db $(INSTALLBIN)
	-cp gap4 $(INSTALLSCRIPT)
	-mkdir $(INSTALLLIB)/gap
	cp *.tcl tclIndex $(INSTALLLIB)/gap
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

distsrc: distsrc_dirs
	-cp -R *.[chf] *.tcl tclIndex Makefile gap4 gap4viewer gap4.bat \
		dependencies $(DIRNAME)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
