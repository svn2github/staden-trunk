LIBS = gap5
PROGS = $(LIBS) $(O)/gap5sh 
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
#SRCROOT=/home5/pubseq/share/src
#SRCROOT=/nfs/arran/home5/pubseq/share/src
#REMOTESRC=$(GAP5SRC)

include	 $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

GAP5SRC=$(SRCROOT)/gap5

# Local comment: Comment out next line for remote compilation
GAP5BIN=$(O)

TGAP_INC=-I$(SRCROOT)/tgap
TGAP_DEP=-ltgap

INCLUDES_E += -I$(GAP5SRC) $(TGAP_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(MUT_INC) \
	      $(F77_INC) $(P3_INC)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)


TGAP_LIB=-ltgap


#
# Building the programs
# This should be just a linking phase because all of the object
# files and library files are generated using implicit rules.
# We use the fortran compiler to do linking.
#
#	$(GAP5BIN)/contig_editor.o
CEDITOR = \
	$(GAP5BIN)/editor_view.o \
	$(GAP5BIN)/editor_join.o \
	$(GAP5BIN)/tkEdNames.o\
	$(GAP5BIN)/tkEditor.o\
	$(GAP5BIN)/tman_interface.o\
	$(GAP5BIN)/tman_display.o

CEDITOR_SCRAP=\
	$(GAP5BIN)/contigEditor.o\
	$(GAP5BIN)/tkEditor.o\
	$(GAP5BIN)/tkEdUtils.o\
	$(GAP5BIN)/edInterface.o\
	$(GAP5BIN)/edUtils2.o \
	$(GAP5BIN)/tagU1.o\
	$(GAP5BIN)/undo.o\
	$(GAP5BIN)/edExtend.o\
	$(GAP5BIN)/edCommands.o\
	$(GAP5BIN)/edMutations.o\
	$(GAP5BIN)/tagEditor.o\
	$(GAP5BIN)/searchUtils.o\
	$(GAP5BIN)/tman_interface.o\
	$(GAP5BIN)/tman_display.o\
	$(GAP5BIN)/tman_cons.o\
	$(GAP5BIN)/tman_diff.o\
	$(GAP5BIN)/join.o\
	$(GAP5BIN)/oligo.o

GAP5=\
	$(GAP5BIN)/newgap5_cmds.o \
	$(GAP5BIN)/gap4_compat.o \
        $(GAP5BIN)/actf.o \
	$(GAP5BIN)/gap_hash.o\
	$(GAP5BIN)/hash_lib.o \
	$(GAP5BIN)/qual.o \
	$(GAP5BIN)/qualIO.o \
	$(GAP5BIN)/gap_globals.o \
	$(GAP5BIN)/tagdb.o\
	$(GAP5BIN)/init.o \
	$(GAP5BIN)/notedb.o\
	$(GAP5BIN)/active_tags.o\
	$(GAP5BIN)/list_proc.o\
	$(GAP5BIN)/gap-error.o\
	$(GAP5BIN)/stack_dump.o\
	$(GAP5BIN)/tk-io-reg.o\
	$(GAP5BIN)/tkAppInit.o \
	$(GAP5BIN)/consen.o\
	$(GAP5BIN)/contig_selector.o\
	$(GAP5BIN)/cs-object.o\
	$(GAP5BIN)/find_repeats.o\
	$(GAP5BIN)/consensus.o\
	$(GAP5BIN)/read_depth.o\
	$(GAP5BIN)/fij.o\
	$(GAP5BIN)/do_fij.o\
	$(CEDITOR)

GAP5SH_LIBS=\
	$(IOLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(MISC_LIB) \
	$(P3_LIB) \
	$(TGAP_LIB)


$(O)/gap5sh:	$(GAP5)
	$(CXX) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(GAP5) $(GAP5SH_LIBS) $(CXX_DEP) $(F77_DEP)


GAP5_LIBS=\
	$(IOLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(BIOLIMS_LIB) \
	$(BIOGAP_LIB) \
	$(TGAP_DEP) \
	$(MUT_LIB) \
	$(MISC_LIB) \
	$(P3_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(GAP5) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LINK_PATHFLAG)$(L) $(GAP5) $(GAP5_LIBS) $(CXX_DEP) $(SHLIB_DEP) $(F77_DEP)


$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(GAP5)
	$(MKDEFL) $@ $(GAP5)


DEPEND_OBJ = $(GAP5)

install:
	-cp gap5 $(INSTALLSCRIPT)
	-mkdir $(INSTALLLIB)/gap5
	cp *.tcl tclIndex $(INSTALLLIB)/gap5
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

distsrc: distsrc_dirs
	-cp -R *.[chf] *.tcl tclIndex Makefile gap5 gap5viewer gap5.bat \
		dependencies $(DIRNAME)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
