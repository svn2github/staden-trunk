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

# DO NOT DELETE THIS LINE -- make depend depends on it.

align.o: $(SRCROOT)/Misc/misc.h
align.o: $(SRCROOT)/Misc/os.h
align.o: $(SRCROOT)/Misc/xalloc.h
align.o: $(PWD)/staden_config.h
align.o: align.h
align.o: align_lib_old.h
align_lib.o: $(SRCROOT)/Misc/misc.h
align_lib.o: $(SRCROOT)/Misc/os.h
align_lib.o: $(SRCROOT)/Misc/xalloc.h
align_lib.o: $(PWD)/staden_config.h
align_lib.o: align_lib.h
align_lib.o: dna_utils.h
align_lib.o: read_matrix.h
align_lib_old.o: $(SRCROOT)/Misc/misc.h
align_lib_old.o: $(SRCROOT)/Misc/os.h
align_lib_old.o: $(SRCROOT)/Misc/xalloc.h
align_lib_old.o: $(PWD)/staden_config.h
align_lib_old.o: align.h
align_lib_old.o: align_lib_old.h
align_ss.o: $(SRCROOT)/Misc/misc.h
align_ss.o: $(SRCROOT)/Misc/os.h
align_ss.o: $(SRCROOT)/Misc/xalloc.h
align_ss.o: $(PWD)/staden_config.h
align_ss.o: align.h
align_ss.o: align_lib_old.h
align_ss2.o: $(SRCROOT)/Misc/misc.h
align_ss2.o: $(SRCROOT)/Misc/os.h
align_ss2.o: $(SRCROOT)/Misc/xalloc.h
align_ss2.o: $(PWD)/staden_config.h
align_ss2.o: align.h
align_ss2.o: align_lib_old.h
align_sv.o: $(SRCROOT)/Misc/misc.h
align_sv.o: $(SRCROOT)/Misc/os.h
align_sv.o: $(SRCROOT)/Misc/xalloc.h
align_sv.o: $(PWD)/staden_config.h
align_sv.o: align.h
align_sv.o: align_lib_old.h
align_sv.o: uascii.gbl
base_comp.o: $(SRCROOT)/Misc/misc.h
base_comp.o: $(SRCROOT)/Misc/os.h
base_comp.o: $(SRCROOT)/Misc/xalloc.h
base_comp.o: $(PWD)/staden_config.h
base_comp.o: base_comp.h
base_comp.o: dna_utils.h
base_comp.o: edge.h
dna_utils.o: $(SRCROOT)/Misc/FtoC.h
dna_utils.o: $(SRCROOT)/Misc/misc.h
dna_utils.o: $(SRCROOT)/Misc/os.h
dna_utils.o: $(SRCROOT)/Misc/xalloc.h
dna_utils.o: $(PWD)/staden_config.h
dna_utils.o: dna_utils.h
edge.o: $(SRCROOT)/Misc/misc.h
edge.o: $(SRCROOT)/Misc/os.h
edge.o: $(SRCROOT)/Misc/xalloc.h
edge.o: $(PWD)/staden_config.h
edge.o: dna_utils.h
filter_words.o: dna_utils.h
filter_words.o: filter_words.h
genetic_code.o: $(SRCROOT)/Misc/array_arith.h
genetic_code.o: $(SRCROOT)/Misc/misc.h
genetic_code.o: $(SRCROOT)/Misc/os.h
genetic_code.o: $(SRCROOT)/Misc/xalloc.h
genetic_code.o: $(SRCROOT)/text_utils/text_output.h
genetic_code.o: $(PWD)/staden_config.h
genetic_code.o: dna_utils.h
genetic_code.o: genetic_code.h
open_reading_frames.o: $(SRCROOT)/Misc/array_arith.h
open_reading_frames.o: $(SRCROOT)/Misc/misc.h
open_reading_frames.o: $(SRCROOT)/Misc/os.h
open_reading_frames.o: $(SRCROOT)/Misc/xalloc.h
open_reading_frames.o: $(SRCROOT)/text_utils/text_output.h
open_reading_frames.o: $(PWD)/staden_config.h
open_reading_frames.o: dna_utils.h
open_reading_frames.o: genetic_code.h
read_matrix.o: $(SRCROOT)/Misc/misc.h
read_matrix.o: $(SRCROOT)/Misc/os.h
read_matrix.o: $(SRCROOT)/Misc/xalloc.h
read_matrix.o: $(PWD)/staden_config.h
read_matrix.o: read_matrix.h
renz_utils.o: $(SRCROOT)/Misc/getfile.h
renz_utils.o: $(SRCROOT)/Misc/misc.h
renz_utils.o: $(SRCROOT)/Misc/os.h
renz_utils.o: $(SRCROOT)/Misc/xalloc.h
renz_utils.o: $(PWD)/staden_config.h
renz_utils.o: dna_utils.h
renz_utils.o: renz_utils.h
scramble.o: sequence_formats.h
search_utils.o: search_utils.h
sequence_formats.o: $(SRCROOT)/Misc/array.h
sequence_formats.o: $(SRCROOT)/Misc/getfile.h
sequence_formats.o: $(SRCROOT)/Misc/misc.h
sequence_formats.o: $(SRCROOT)/Misc/os.h
sequence_formats.o: $(SRCROOT)/Misc/xalloc.h
sequence_formats.o: $(SRCROOT)/Misc/xerror.h
sequence_formats.o: $(PWD)/staden_config.h
sequence_formats.o: sequence_formats.h
