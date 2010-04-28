#
# Makefile for vector_clip
#

PROGS = vector_clip

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC) $(SEQUTILS_INC) $(MISC_INC) $(TEXTUTILS_INC)

#EXTRA_LIBS += -ldmalloc

OBJS	= \
	vector_clip.o

HLIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(IOUTILS_LIB) \
	$(MISC_LIB) \
	$(IOLIB_LIB) \
	$(MATH_LIB)


vector_clip: $(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(HLIBS) $(LIBSC) 


DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE

vector_clip.o: $(SRCROOT)/Misc/getfile.h
vector_clip.o: $(SRCROOT)/Misc/misc.h
vector_clip.o: $(SRCROOT)/Misc/os.h
vector_clip.o: $(SRCROOT)/Misc/xalloc.h
vector_clip.o: $(PWD)/staden_config.h
vector_clip.o: $(SRCROOT)/seq_utils/dna_utils.h
