#
# Makefile for clip
#

PROGS= 	$(O)/qclip

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC) $(TKUTILS_INC)

QCLIPSRC=$(O)

all: $(O)/qclip

OBJ=\
	$(QCLIPSRC)/qclip.o\
	$(QCLIPSRC)/consen.o \
	$(QCLIPSRC)/seqInfo.o

QCLIP_LIBS=\
	$(TEXTUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB)

$(O)/qclip: $(OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJ) $(QCLIP_LIBS) $(LIBSC)

DEPEND_OBJ = $(OBJ)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
