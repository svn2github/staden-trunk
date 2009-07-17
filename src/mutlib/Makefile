#
# Makefile for Mutation Library
#
LIBS = mutlib
PROGS= $(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

CXXFLAGS += $(SHLIB_CXXFLAGS) $(MISC_DLL)

INCLUDES_E += $(MISC_INC) $(MUT_INC) $(SEQUTILS_INC) $(IOLIB_INC)

DEFINES += -DMUTLIB_EXPORTS -DDYNMAT

#
# Objects
#
OBJS = \
	align.o \
	basecall.o \
	caller.o \
	caller_snr.o \
	caller_base.o \
	caller_level.o \
	dnatable.o \
	muttag.o \
	mutationtag.o \
	mutationtag_utils.o \
	pathutil.o \
	peakcall.o \
	stringlist.o \
	tagarray.o \
	trace.o \
	validate.o \
	read_matrix.o \
	sp_align_lib.o \
	sp_alignment.o \
	sp_hash_lib.o \
	sp_matrix.o \
	tracediff.o \
	tracediff_helper.o \
	tracediff_mutations.o \
	tracediff_parameters.o \
	tracediff_validate.o \
	tracealign.o \
	tracealign_cache.o \
	tracealign_helper.o \
	tracealign_interpolate.o \
	tracealign_preprocess.o \
	tracealign_validate.o \
	mutscan.o \
	mutscan_helper.o \
	mutscan_validate.o \
	mutscan_analysis.o \
	mutscan_preprocess.o \
	mutscan_parameters.o


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

install:
	cp $(PROGLIBS) $(INSTALLLIB)

include dependencies
