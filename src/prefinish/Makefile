LIBS = prefinish
PROGS = x$(LIBS) insert_size_dist read_size_dist read_conf_dist
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E := $(MISC_INC) $(GAP4_INC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(P3_INC) \
	      $(INCLUDES_E)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	finish.o \
	finish_hash.o \
	finish_long.o \
	finish_main.o \
	finish_utils.o \
	finish_walk.o \
	finish_reverse.o \
	finish_filter.o \
	finish_pcr.o \
	dust.o

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

ISLIBS=\
	$(G_LIB) \
	$(GAP_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB) \
	$(TCL_LIB)

ISOBJS= insert_size_dist.o

RSOBJS=	read_size_dist.o

RCOBJS=	read_conf_dist.o

insert_size_dist: $(ISOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(ISOBJS) $(ISLIBS) $(LIBSC)

read_size_dist: $(RSOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(RSOBJS) $(ISLIBS) $(LIBSC)

read_conf_dist: $(RCOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(RCOBJS) $(ISLIBS) $(LIBSC)

S=$(SRCROOT)/prefinish
install:
	-mkdir $(INSTALLTCL)/prefinish
	cp $(S)/*.tcl $(S)/*.itcl $(S)/*.itk $(S)/help $(S)/args.template $(S)/tclIndex $(INSTALLTCL)/prefinish
	cp $(S)/config.prefinish $(INSTALLTCL)/prefinish
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/prefinish $(INSTALLBIN)
	cp $(S)/finish_sanger $(S)/finish_cDNA $(S)/finish_cDNA_ends_only \
	    $(INSTALLSCRIPTS)
	-mkdir $(INSTALLTCL)/prefinish/rules
	cp $(S)/rules/*.rule $(INSTALLTCL)/prefinish/rules
	cp $(S)/config.prefinish $(S)/prefinishrc $(INSTALLETC)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/*.itcl $(S)/*.itk $(S)/help \
	   $(S)/args.template $(S)/Makefile $(S)/prefinish $(S)/tclIndex \
	   $(S)/config.prefinish $(S)/finish_sanger $(S)/finish_cDNA \
	   $(S)/finish_cDNA_ends_only $(S)/prefinishrc \
	   $(DIRNAME)
	-mkdir $(DIRNAME)/rules
	cp $(S)/rules/* $(DIRNAME)/rules
	-mkdir $(DIRNAME)/docs
	cp $(S)/docs/* $(DIRNAME)/docs

DEPEND_OBJ = $(OBJS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
