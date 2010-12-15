#
# Makefile for make_weights
#

PROGS = make_weights

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(SEQUTILS_INC) $(MISC_INC) $(TEXTUTILS_INC)

#EXTRA_LIBS += -ldmalloc

OBJS	= make_weights.o

HLIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB) \
	$(MATH_LIB) \
	$(IOLIB_LIB)
	


make_weights: $(OBJS) $(DEPS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJS) $(HLIBS) $(LIBSC) 
#	purify $(CLD) -o $@ $(OBJS) $(HLIBS) $(LIBSC) 



DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it

make_weights.o: $(SRCROOT)/Misc/misc.h
make_weights.o: $(SRCROOT)/Misc/os.h
make_weights.o: $(SRCROOT)/Misc/xalloc.h
make_weights.o: $(PWD)/staden_config.h
make_weights.o: $(SRCROOT)/seq_utils/dna_utils.h
