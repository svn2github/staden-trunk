#
# Makefile for seq_utils routines
#

LIBS 	= seq_utils
PROGS	= $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

# 7/1/990 johnt - define SEQ_UTILS_DLL to cause export of global variables with Visual C++
CFLAGS += $(SHLIB_CFLAGS) $(SEQ_UTILS_DLL)
SEQUTILSBIN=$(O)

INCLUDES_E += $(TEXTUTILS_INC)

#
# Objects
#
OBJS = \
	$(SEQUTILSBIN)/align.o \
	$(SEQUTILSBIN)/align_lib_old.o \
	$(SEQUTILSBIN)/align_ss.o \
	$(SEQUTILSBIN)/align_ss2.o \
	$(SEQUTILSBIN)/align_sv.o \
	$(SEQUTILSBIN)/dna_utils.o \
	$(SEQUTILSBIN)/genetic_code.o \
	$(SEQUTILSBIN)/renz_utils.o \
	$(SEQUTILSBIN)/sequence_formats.o \
	$(SEQUTILSBIN)/scramble.o \
	$(SEQUTILSBIN)/base_comp.o \
	$(SEQUTILSBIN)/open_reading_frames.o\
	$(SEQUTILSBIN)/edge.o\
	$(SEQUTILSBIN)/search_utils.o\
	$(SEQUTILSBIN)/align_lib.o\
	$(SEQUTILSBIN)/read_matrix.o

# 7/1/99 johnt - added TKUTILS_DEP (for Visual C++)
DEPS = \
	$(TEXTUTILS_DEP) \
	$(TKUTILS_DEP) \
	$(MISC_DEP)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(DEPS)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for WINNT
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	-cp -R *.[ch] *.gbl Makefile dependencies $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
# DO NOT DELETE
