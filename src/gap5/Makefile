LIBS = gap5
PROGS = $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

TGAP_INC=-I$(SRCROOT)/tgap
TGAP_DEP=-ltgap

INCLUDES_E += $(MISC_INC) $(TGAP_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC)
	      
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)


TGAP_LIB=-ltgap


#
# Building the programs
# This should be just a linking phase because all of the object
# files and library files are generated using implicit rules.
# We use the fortran compiler to do linking.
CEDITOR = \
	editor_view.o \
	editor_join.o \
	tkEdNames.o\
	tkEditor.o\
	tman_interface.o\
	tman_display.o

CEDITOR_SCRAP=\
	contigEditor.o\
	tkEditor.o\
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
	join.o\
	oligo.o

GAP5=\
	newgap5_cmds.o \
	gap4_compat.o \
        actf.o \
	gap_hash.o\
	hash_lib.o \
	qual.o \
	qualIO.o \
	gap_globals.o \
	tagdb.o\
	init.o \
	notedb.o\
	active_tags.o\
	list_proc.o\
	gap-error.o\
	stack_dump.o\
	tk-io-reg.o\
	tkAppInit.o \
	consen.o\
	contig_selector.o\
	cs-object.o\
	find_repeats.o\
	consensus.o\
	read_depth.o\
	template_display.o \
	fij.o\
	do_fij.o\
	break_contig.o \
	export_contigs.o \
	$(CEDITOR)

GAP5_LIBS=\
	$(IOLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(TGAP_DEP) \
	$(MISC_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(GAP5) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LINK_PATHFLAG)$(L) $(GAP5) $(GAP5_LIBS) $(SHLIB_DEP)


$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(GAP5)
	$(MKDEFL) $@ $(GAP5)


DEPEND_OBJ = $(GAP5)

S=$(SRCROOT)/gap5
install:
	-cp $(S)/gap5 $(INSTALLBIN)
	-mkdir $(INSTALLTCL)/gap5
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/gap5
	cp $(S)/gap5rc $(S)/gap5rc_menu_full $(INSTALLETC)
	cp $(PROGLIBS) $(INSTALLLIB)

distsrc: distsrc_dirs
	-cp -R *.[chf] *.tcl tclIndex Makefile gap5 gap5viewer gap5.bat \
		dependencies $(DIRNAME)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
