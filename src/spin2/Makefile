#
# Makefile for spin2 routines
#
SPIN2_LIB     	= $(SPIN2_LIB_S) $(LINK_LIBFLAG)spin2$(LIB_EXT) $(SPIN2_LIB_E)
SPIN2_INC	= -I$(SPIN2SRC)
SPIN2SRC	= $(SRCROOT)/spin2
SPIN2BIN	= $(SPIN2SRC)/$(O)

LIBS = spin2
PROGS = lib$(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include	$(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(TK_INC) $(ITK_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(IOLIB_INC)
CFLAGS += $(SHLIB_CFLAGS)

SPIN2BIN=$(O)

OBJS = \
	$(SPIN2BIN)/seq_reg.o \
	$(SPIN2BIN)/seq_reg_cmds.o \
	$(SPIN2BIN)/sequtils.o \
	$(SPIN2BIN)/seq_results.o \
	$(SPIN2BIN)/seq_sendto.o \
	$(SPIN2BIN)/spin_globals.o \
	$(SPIN2BIN)/seqed.o\
	$(SPIN2BIN)/tkSeqed.o \
	$(SPIN2BIN)/tkSeqedUtils.o \
	$(SPIN2BIN)/tkSeqedNames.o \
	$(SPIN2BIN)/seqed_write.o \
	$(SPIN2BIN)/seqedInterface.o\
	$(SPIN2BIN)/seqed_translate.o\
	$(SPIN2BIN)/seqed_restriction_enzymes.o\
	$(SPIN2BIN)/seqed_search.o\
	$(SPIN2BIN)/sequence_pair_display.o\
	$(SPIN2BIN)/emboss_input_funcs.o \
	$(SPIN2BIN)/readpam.o\
	$(SPIN2BIN)/compare_spans.o \
	$(SPIN2BIN)/sip_hash.o \
	$(SPIN2BIN)/sip_similar_spans.o\
	$(SPIN2BIN)/sip_find_identity.o\
	$(SPIN2BIN)/sip_quick_scan.o\
	$(SPIN2BIN)/sip_align.o \
	$(SPIN2BIN)/sip_sim.o\
	$(SPIN2BIN)/sip_globals.o \
	$(SPIN2BIN)/sip_cmds.o \
	$(SPIN2BIN)/sip_results.o \
	$(SPIN2BIN)/probs.o \
	$(SPIN2BIN)/sim.o \
	$(SPIN2BIN)/init.o\
	$(SPIN2BIN)/nip_globals.o \
	$(SPIN2BIN)/nip_cmds.o \
	$(SPIN2BIN)/nip_base_comp.o \
	$(SPIN2BIN)/nip_gene_search.o \
	$(SPIN2BIN)/nip_trna_search.o\
	$(SPIN2BIN)/codon_content.o \
	$(SPIN2BIN)/nip_stop_codon.o \
	$(SPIN2BIN)/trna_search.o \
	$(SPIN2BIN)/splice_search.o \
	$(SPIN2BIN)/nip_canvas_box.o \
	$(SPIN2BIN)/nip_splice_search.o \
	$(SPIN2BIN)/nip_string_search.o \
	$(SPIN2BIN)/nip_wtmatrix_search.o \
	$(SPIN2BIN)/nip_restriction_enzymes.o\
	$(SPIN2BIN)/dinuc_freqs.o\
	$(SPIN2BIN)/ft_viewer.o\
	$(SPIN2BIN)/feature_selector.o\
	$(SPIN2BIN)/seq_element.o\
	$(SPIN2BIN)/seq_element_cmds.o\
	$(SPIN2BIN)/spin_cli_arg.o

SPIN2OBJS2=\
	$(TKUTILSBIN)/stash.o\
	$(TKUTILSBIN)/tkMain.o

#
# Main dependency
#
lib$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

DEPS = $(TCL_LIB)\
       $(TK_LIB) \
	$(ITCL_LIB)\
       $(ITK_LIB) \
       $(IOUTILS_LIB) \
       $(TKUTILS_LIB) \
       $(SEQUTILS_LIB) \
       $(MISC_LIB)


$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(DEPS)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for windows
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl *.itk tclIndex Makefile spin2 spin2.bat dependencies\
		$(DIRNAME)

install:
	-cp spin2 $(INSTALLSCRIPT)
	-mkdir $(INSTALLLIB)/spin2
	cp *.tcl *.itk tclIndex $(INSTALLLIB)/spin2
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
# DO NOT DELETE

