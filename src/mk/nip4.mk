#LIBS = nip
#PROGS = $(LIBS)
#PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include	$(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

NIP4BIN=$(O)

INCLUDES_E += -I$(GAP4SRC) -I$(SIP4SRC) $(TK_INC) $(IOLIB_INC) $(TCLUTILS_INC) \
	      $(TKUTILS_INC) $(SEQUTILS_INC) $(MISC_INC) $(SEQLIB_INC) \
	$(SPIN_INC) $(G_INC)

NIPOBJS = \
	$(NIP4BIN)/init.o\
	$(NIP4BIN)/nip_globals.o \
	$(NIP4BIN)/nip_cmds.o \
	$(NIP4BIN)/nip_base_comp.o \
	$(NIP4BIN)/nip_gene_search.o \
	$(NIP4BIN)/array_arith.o \
	$(NIP4BIN)/codon_content.o \
	$(NIP4BIN)/nip_stop_codon.o \
	$(NIP4BIN)/trna_search.o \
	$(NIP4BIN)/nip_trna_search.o \
	$(NIP4BIN)/splice_search.o \
	$(NIP4BIN)/nip_splice_search.o \
	$(NIP4BIN)/nip_string_search.o \
	$(NIP4BIN)/nip_wtmatrix_search.o \
	$(NIP4BIN)/nip_restriction_enzymes.o\
	$(NIP4BIN)/sequtils_cmds.o\
	$(NIP4BIN)/interp_global.o\
	$(NIP4BIN)/nip_canvas_box.o\
	$(NIP4BIN)/dinuc_freqs.o

NIPOBJS2=\
	$(NIP4BIN)/tkAppInit.o

#	$(GAP4BIN)/tkMain.o

# 12/1/99 johnt - added additional dependencies for WINNT support
DEPS=\
	$(IOLIB_DEP) \
	$(SEQLIB_DEP) \
	$(MISC_DEP) \
	$(SEQUTILS_DEP) \
	$(SPIN_DEP) \
	$(TKUTILS_DEP) \
	$(TK_DEP) \
	$(TCL_DEP)

NIP4_LIBS=\
	$(SEQLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(SPIN_LIB) \
	$(TCLUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB)



$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@


$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(NIPOBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(NIPOBJS) $(DEPS)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for WINNT
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(NIPOBJS)
	$(MKDEFL) $@ $(NIPOBJS)


DEPEND_OBJ = $(NIP4OBJS)

distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl tclIndex Makefile dependencies \
		$(DIRNAME)
	-rm $(DIRNAME)/resource.h

install:
	-mkdir $(INSTALLLIB)/nip
	cp *.tcl tclIndex $(INSTALLLIB)/nip
#	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.


