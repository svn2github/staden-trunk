LIBS = gap
ifeq ($(MACHINE),windows)
PROGS = $(O)/copy_db $(LIBS)
else
PROGS = $(O)/gap4sh $(O)/copy_db $(LIBS)
endif
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
#SRCROOT=/home5/pubseq/share/src
#SRCROOT=/nfs/arran/home5/pubseq/share/src
#REMOTESRC=$(GAP4SRC)

include	 $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

# Local comment: Comment out next line for remote compilation
GAP4BIN=$(O)

INCLUDES_E += -I$(GAP4SRC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

#
# Common objects, needed by every program
# 7/1/99 johnt - these objects no longer used by anything
#
#COMMONOBJS=\
#	$(STADENBIN)/seeme-$(MACHINE).o\
#	$(STADENBIN)/Cme-$(MACHINE).o

#
# Building the programs
# This should be just a linking phase because all of the object
# files and library files are generated using implicit rules.
# We use the fortran compiler to do linking.
#

#
# Trace manager objects
#

include $(SRCROOT)/mk/gap4_defs.mk

GAPDB=\
	$(GAPDB_LOW)\
	$(GAPDB_MID)\
	$(GAP4BIN)/IO2.o\
	$(GAP4BIN)/seqInfo.o\
	$(GAP4BIN)/IO3.o \
	$(GAP4BIN)/io_utils.o \
	$(GAP4BIN)/io_handle.o \
	$(GAP4BIN)/gap-tcl.o \
	$(GAP4BIN)/tk-io-reg.o

CONSEN=\
	$(GAP4BIN)/consen.o\
	$(GAP4BIN)/gap_hash.o\
	$(GAP4BIN)/qual.o\
	$(GAP4BIN)/qualIO.o

OSP=\
	$(GAP4BIN)/osp_analysis.o\
	$(GAP4BIN)/oligocom.o

TCLBITS=\
	$(GAP4BIN)/gap_cli_arg.o \
	$(GAP4BIN)/gap_globals.o \
	$(GAP4BIN)/newgap_cmds.o \
	$(GAP4BIN)/init.o

EDITOR = \
	$(GAP4BIN)/tkEditor.o\
	$(GAP4BIN)/tkEdNames.o\
	$(GAP4BIN)/tkEdUtils.o\
	$(GAP4BIN)/edInterface.o\
	$(GAP4BIN)/edUtils2.o \
	$(GAP4BIN)/tagU1.o\
	$(GAP4BIN)/undo.o\
	$(GAP4BIN)/edExtend.o\
	$(GAP4BIN)/edCommands.o\
	$(GAP4BIN)/tagEditor.o\
	$(GAP4BIN)/searchUtils.o\
	$(GAP4BIN)/tman_interface.o\
	$(GAP4BIN)/tman_display.o\
	$(GAP4BIN)/tman_cons.o\
	$(GAP4BIN)/tman_diff.o\
	$(GAP4BIN)/contigEditor.o\
	$(GAP4BIN)/join.o\
	$(GAP4BIN)/oligo.o

GAPFUNCS = \
	$(GAP4BIN)/legacy.o\
	$(GAP4BIN)/bubbl3.o\
	$(GAP4BIN)/parse_db.o\
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
	$(GAP4BIN)/hash_lib.o \
	$(GAP4BIN)/do_fij.o \
	$(GAP4BIN)/auto_assemble.o \
	$(GAP4BIN)/dis_readings.o \
	$(GAP4BIN)/find_repeats.o \
	$(GAP4BIN)/break_contig.o \
	$(GAP4BIN)/quality_plot.o \
	$(GAP4BIN)/readpair.o \
	$(GAP4BIN)/contig_selector.o \
	$(GAP4BIN)/complement.o \
	$(GAP4BIN)/cs-object.o\
	$(GAP4BIN)/list_proc.o\
	$(GAP4BIN)/dstrand.o\
	$(GAP4BIN)/oligo_sel.o\
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
	$(GAP4BIN)/vseqs.o

GAPOBJS=\
	$(COMMONOBJS)\
	$(GAPDB)\
	$(TCLBITS)\
	$(CONSEN)\
	$(GAPFUNCS)\
	$(OSP)\
	$(EDITOR)

GAPOBJS2=\
	$(GAP4BIN)/tkAppInit.o

#GAPDEP=\
#	$(IOLIBBIN)/libexp.a \
#	$(IOLIBBIN)/libread.a \
#	$(IOLIBBIN)/libio-utils.a \
#	$(SEQUTILSBIN)/libseq_utils.a \
#	$(TCLUTILSBIN)/libtcl_utils.a \
#	$(MISCBIN)/libmisc.a \
#	$(GBIN)/libg.a
GAPDEP=\
	$(EXP_DEP) \
	$(IOLIB_DEP) \
	$(SEQUTILS_DEP) \
	$(G_DEP) \
	$(TKUTILS_DEP) \
	$(MISC_DEP) \
	$(F77_DEP)

GAP4_LIBS=\
	$(IOLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(MISC_LIB) \
	$(TK_LIB) \
	$(BIOLIMS_LIB) \
	$(BIOGAP_LIB) \
	$(G_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(GAPOBJS) $(DEF_FILE)
	-$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LINK_PATHFLAG)$(L) $(GAPOBJS) $(GAP4_LIBS) $(LIBSF) $(GAPDEP) $(SHLIB_DEP)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for WINNT
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(GAPOBJS)
	$(MKDEFL) $@ $(GAPOBJS)

$(O)/gap4sh:	$(GAPOBJS) $(GAPOBJS2)
	$(CXX) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(GAPOBJS) $(GAPOBJS2) $(GAP4_LIBS) $(LIBSF) $(GAPDEP)

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

CPDEP=\
	$(G_DEP) \
	$(TEXTUTILS_DEP) \
	$(MISC_DEP) \
	$(TCL_DEP)

$(O)/copy_db:   $(COPYDBOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(COPYDBOBJS) $(CPLIB) $(LIBSC) $(CPDEP)

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
	$(CLD) -o $@ $(THRASHOBJS) $(TLIB) $(LIBSC) $(CPDEP) -ldmalloc
$(O)/thrashb:   $(THRASHBOBJS)
	$(CLD) -o $@ $(THRASHBOBJS) $(TLIB) $(LIBSC) $(CPDEP) -ldmalloc

DEPEND_OBJ = $(GAPOBJS) $(COPYDBOBJS)

install:
	cp $(O)/gap4sh $(O)/copy_db $(INSTALLBIN)
	cp gap4 $(INSTALLSCRIPT)
	-mkdir $(INSTALLLIB)/gap
	cp *.tcl tclIndex $(INSTALLLIB)/gap
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

distsrc: distsrc_dirs
	-cp -R *.[chf] *.tcl tclIndex Makefile gap4 gap4.bat dependencies \
		$(DIRNAME)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
