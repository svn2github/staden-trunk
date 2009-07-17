#
# Makefile for clip
#

PROGS= 	qclip

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC) $(TKUTILS_INC)

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

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
