#
# Makefile for spin routines
#

LIBS = spin
PROGS = lib$(LIBS) 
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(TK_INC)  $(TKUTILS_INC) $(SEQUTILS_INC) $(IOLIB_INC)
CFLAGS += $(SHLIB_CFLAGS)

SPINBIN=$(O)

OBJS = \
	$(SPINBIN)/seq_reg.o \
	$(SPINBIN)/seq_reg_cmds.o \
	$(SPINBIN)/raster_cmds.o \
	$(SPINBIN)/seq_raster.o \
	$(SPINBIN)/sequtils.o \
	$(SPINBIN)/seq_results.o \
	$(SPINBIN)/seq_sendto.o \
	$(SPINBIN)/spin_globals.o \
	$(SPINBIN)/seq_plot_funcs.o\
	$(SPINBIN)/seqed.o\
	$(SPINBIN)/tkSeqed.o \
	$(SPINBIN)/tkSeqedUtils.o \
	$(SPINBIN)/tkSeqedNames.o \
	$(SPINBIN)/seqed_write.o \
	$(SPINBIN)/seqedInterface.o\
	$(SPINBIN)/seqed_translate.o\
	$(SPINBIN)/seqed_restriction_enzymes.o\
	$(SPINBIN)/seqed_search.o\
	$(SPINBIN)/sequence_pair_display.o\
	$(SPINBIN)/emboss_input_funcs.o \
	$(SPINBIN)/readpam.o\
	$(SPINBIN)/compare_spans.o \
	$(SPINBIN)/sip_hash.o \
	$(SPINBIN)/sip_similar_spans.o\
	$(SPINBIN)/sip_find_identity.o\
	$(SPINBIN)/sip_quick_scan.o\
	$(SPINBIN)/sip_align.o \
	$(SPINBIN)/sip_globals.o \
	$(SPINBIN)/sip_cmds.o \
	$(SPINBIN)/sip_results.o \
	$(SPINBIN)/rescan_matches.o \
	$(SPINBIN)/probs.o \
	$(SPINBIN)/sim.o \
	$(SPINBIN)/sip_sim.o\
	$(SPINBIN)/init.o\
	$(SPINBIN)/nip_globals.o \
	$(SPINBIN)/nip_cmds.o \
	$(SPINBIN)/nip_base_comp.o \
	$(SPINBIN)/nip_gene_search.o \
	$(SPINBIN)/nip_trna_search.o\
	$(SPINBIN)/codon_content.o \
	$(SPINBIN)/nip_stop_codon.o \
	$(SPINBIN)/trna_search.o \
	$(SPINBIN)/splice_search.o \
	$(SPINBIN)/nip_splice_search.o \
	$(SPINBIN)/nip_string_search.o \
	$(SPINBIN)/nip_wtmatrix_search.o \
	$(SPINBIN)/nip_restriction_enzymes.o\
	$(SPINBIN)/nip_canvas_box.o\
	$(SPINBIN)/dinuc_freqs.o

#
# Main dependency
#
lib$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

DEPS = $(TCL_LIB)\
       $(TK_LIB) \
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
	-cp -R *.[ch] *.tcl *.itk tclIndex Makefile spin spin.bat dependencies\
		$(DIRNAME)

install:
	-cp spin $(INSTALLSCRIPT)
	-mkdir $(INSTALLLIB)/spin
	cp *.tcl *.itk tclIndex $(INSTALLLIB)/spin
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
# DO NOT DELETE

