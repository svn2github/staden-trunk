LIBS = prefinish
PROGS = x$(LIBS) $(O)/insert_size_dist $(O)/read_size_dist $(O)/read_conf_dist
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..

include	 $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

# Local comment: Comment out next line for remote compilation
FINBIN=$(O)

#L=$(O)

INCLUDES_E := $(GAP4_INC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(P3_INC) \
	      $(INCLUDES_E)
#INCLUDES += -Iprimer3_0_9_test/src/
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	$(FINBIN)/finish.o \
	$(FINBIN)/finish_hash.o \
	$(FINBIN)/finish_long.o \
	$(FINBIN)/finish_main.o \
	$(FINBIN)/finish_utils.o \
	$(FINBIN)/finish_walk.o \
	$(FINBIN)/finish_reverse.o \
	$(FINBIN)/finish_filter.o \
	$(FINBIN)/dust.o

FIN_DEP=\
	$(TK_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(MISC_LIB) \
	$(GAP_LIB) \
	$(P3_LIB)

x$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(FIN_DEP) $(SHLIB_DEP)

include $(SRCROOT)/mk/gap4_defs.mk
ISLIBS=\
	$(G_LIB) \
	$(GAP_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB) \
	$(TCL_LIB)

ISOBJS=\
	$(O)/insert_size_dist.o

RSOBJS=\
	$(O)/read_size_dist.o

RCOBJS=\
	$(O)/read_conf_dist.o

$(O)/insert_size_dist: $(ISOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(ISOBJS) $(ISLIBS) $(LIBSC)

$(O)/read_size_dist: $(RSOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(RSOBJS) $(ISLIBS) $(LIBSC)

$(O)/read_conf_dist: $(RCOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(RCOBJS) $(ISLIBS) $(LIBSC)

install:
	-mkdir $(INSTALLLIB)/prefinish
	cp *.tcl *.itcl *.itk help args.template tclIndex $(INSTALLLIB)/prefinish
	cp config.prefinish $(INSTALLLIB)/prefinish
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)
	cp prefinish finish_sanger finish_cDNA finish_cDNA_ends_only \
	    $(INSTALLBIN)
	-mkdir $(INSTALLLIB)/prefinish/rules
	cp rules/*.rule $(INSTALLLIB)/prefinish/rules

distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl *.itcl *.itk help args.template Makefile \
	    dependencies CHEMISTRY METHODS README prefinish rules \
	    config.prefinish $(DIRNAME)

DEPEND_OBJ = $(OBJS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
