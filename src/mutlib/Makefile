#
# Makefile for Mutation Library
#
LIBS = mutlib
PROGS= $(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

CFLAGS += $(SHLIB_CXXFLAGS) $(MISC_DLL)

INCLUDES_E += $(MUT_INC) $(SEQUTILS_INC)

DEFINES += -DMUTLIB_EXPORTS -DDYNMAT

MUTBIN=$(O)

#
# Objects
#
OBJS = \
	$(MUTBIN)/align.o \
	$(MUTBIN)/basecall.o \
	$(MUTBIN)/caller.o \
	$(MUTBIN)/caller_snr.o \
	$(MUTBIN)/caller_base.o \
	$(MUTBIN)/caller_level.o \
	$(MUTBIN)/dnatable.o \
	$(MUTBIN)/muttag.o \
	$(MUTBIN)/mutationtag.o \
	$(MUTBIN)/mutationtag_utils.o \
	$(MUTBIN)/pathutil.o \
	$(MUTBIN)/peakcall.o \
	$(MUTBIN)/stringlist.o \
	$(MUTBIN)/tagarray.o \
	$(MUTBIN)/trace.o \
	$(MUTBIN)/validate.o \
	$(MUTBIN)/read_matrix.o \
	$(MUTBIN)/sp_align_lib.o \
	$(MUTBIN)/sp_alignment.o \
	$(MUTBIN)/sp_hash_lib.o \
	$(MUTBIN)/sp_matrix.o \
	$(MUTBIN)/tracediff.o \
	$(MUTBIN)/tracediff_helper.o \
	$(MUTBIN)/tracediff_mutations.o \
	$(MUTBIN)/tracediff_parameters.o \
	$(MUTBIN)/tracediff_validate.o \
	$(MUTBIN)/tracealign.o \
	$(MUTBIN)/tracealign_cache.o \
	$(MUTBIN)/tracealign_helper.o \
	$(MUTBIN)/tracealign_interpolate.o \
	$(MUTBIN)/tracealign_preprocess.o \
	$(MUTBIN)/tracealign_validate.o \
	$(MUTBIN)/mutscan.o \
	$(MUTBIN)/mutscan_helper.o \
	$(MUTBIN)/mutscan_validate.o \
	$(MUTBIN)/mutscan_analysis.o \
	$(MUTBIN)/mutscan_preprocess.o \
	$(MUTBIN)/mutscan_parameters.o


MUT_DEP=\
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB) \
	$(MATH_LIB)

ifeq ($(MACHINE),macosx)
MUT_DEP += $(TK_LIB)
endif

#ifeq ($(MACHINE),alpha)
#MUT_DEP += -Bstatic -lgcc
#endif

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS)
	-cp $(MUTSRC)/$(LIBS).def $(DEF_FILE) 2>&-
	$(SHLIB_LDXX) $(SHLIB_LDXXFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(MUT_DEP)

DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	-cp -R *.[ch] *.[ch]pp Makefile dependencies $(DIRNAME)
	mkdir $(DIRNAME)/docs
	cp docs/* $(DIRNAME)/docs

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
