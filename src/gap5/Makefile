LIBS = gap5
PROGS = $(LIBS) $(O)/gap5sh 
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

TGAP_INC=-I$(SRCROOT)/tgap
TGAP_DEP=-ltgap

INCLUDES_E += -I$(GAP4SRC) $(TGAP_INC) $(TK_INC) \
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
	$(GAP4BIN)/consensus.o\
	$(GAP4BIN)/read_depth.o\
	$(GAP4BIN)/fij.o\
	$(GAP4BIN)/do_fij.o\
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
