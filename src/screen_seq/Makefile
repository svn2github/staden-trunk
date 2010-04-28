#
# Makefile for screen_seq
#

PROGS = screen_seq

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(TEXTUTILS_INC) $(IOLIB_INC) $(SEQUTILS_INC) $(MISC_INC)

OBJS	= screen_seq.o

HLIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(IOUTILS_LIB) \
	$(MISC_LIB) \
	$(IOLIB_LIB) \
	$(MATH_LIB)


screen_seq: $(OBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(HLIBS) $(LIBSC)


DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it.

screen_seq.o: $(SRCROOT)/Misc/getfile.h
screen_seq.o: $(SRCROOT)/Misc/misc.h
screen_seq.o: $(SRCROOT)/Misc/os.h
screen_seq.o: $(SRCROOT)/Misc/xalloc.h
screen_seq.o: $(PWD)/staden_config.h
screen_seq.o: $(SRCROOT)/seq_utils/dna_utils.h
