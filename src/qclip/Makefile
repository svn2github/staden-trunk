#
# Makefile for clip
#

PROGS= 	qclip

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC) $(TKUTILS_INC) $(MISC_INC)

OBJ=\
	qclip.o\
	consen.o \
	seqInfo.o

QCLIP_LIBS=\
	$(TEXTUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB)

qclip: $(OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJ) $(QCLIP_LIBS) $(LIBSC)

DEPEND_OBJ = $(OBJ)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it.

consen.o: $(SRCROOT)/Misc/misc.h
consen.o: $(SRCROOT)/Misc/os.h
consen.o: $(SRCROOT)/Misc/xalloc.h
consen.o: $(PWD)/staden_config.h
consen.o: $(SRCROOT)/qclip/consen.h
qclip.o: $(SRCROOT)/Misc/misc.h
qclip.o: $(SRCROOT)/Misc/os.h
qclip.o: $(SRCROOT)/Misc/xalloc.h
qclip.o: $(PWD)/staden_config.h
qclip.o: $(SRCROOT)/qclip/consen.h
qclip.o: $(SRCROOT)/qclip/seqInfo.h
seqInfo.o: $(SRCROOT)/Misc/misc.h
seqInfo.o: $(SRCROOT)/Misc/os.h
seqInfo.o: $(SRCROOT)/Misc/xalloc.h
seqInfo.o: $(PWD)/staden_config.h
seqInfo.o: $(SRCROOT)/qclip/seqInfo.h
