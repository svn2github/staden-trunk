#
# Makefile for seq_utils routines
#

LIBS 	= seq_utils
PROGS	= $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

CFLAGS += $(SHLIB_CFLAGS) $(SEQ_UTILS_DLL)

INCLUDES_E += $(MISC_INC) $(TEXTUTILS_INC)

#
# Objects
#
OBJS = \
	align.o \
	align_lib_old.o \
	align_ss.o \
	align_ss2.o \
	align_sv.o \
	dna_utils.o \
	genetic_code.o \
	renz_utils.o \
	sequence_formats.o \
	scramble.o \
	base_comp.o \
	open_reading_frames.o\
	edge.o\
	search_utils.o\
	align_lib.o\
	read_matrix.o\
	filter_words.o


SU_LIBS = \
	$(TKUTILS_LIB) \
	$(MISC_LIB)

ifeq ($(MACHINE),macosx)
SU_LIBS += $(TK_LIB)
endif

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@


$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(SU_LIBS)


$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.gbl $(S)/Makefile $(DIRNAME)
	cp $(S)/align_lib_nuc_matrix $(S)/nuc_matrix $(DIRNAME)

S=$(SRCROOT)/seq_utils
install:
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/align_lib_nuc_matrix $(S)/nuc_matrix $(INSTALLETC)

include dependencies
# DO NOT DELETE
