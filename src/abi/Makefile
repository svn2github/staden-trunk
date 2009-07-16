#
# Makefile for getABISampleName
#

PROGS = $(O)/getABISampleName $(O)/getABIcomment $(O)/getABIdate \
	$(O)/getABIfield

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC)
ABIBIN=$(O)

OBJSSN = $(ABIBIN)/abiIO.o $(ABIBIN)/getABISampleName.o
OBJSD  = $(ABIBIN)/abiIO.o $(ABIBIN)/getABIdate.o
OBJSFI = $(ABIBIN)/abiIO.o $(ABIBIN)/getABIfield.o
OBJSC  = $(ABIBIN)/abiIO.o $(ABIBIN)/getABIcomment.o

LIBS=\
	$(IOLIB_LIB) \
	$(MISC_LIB)

$(O)/getABISampleName: $(OBJSSN)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSSN) $(LIBS) $(LIBSC)

$(O)/getABIdate: $(OBJSD)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(IOUTILS_LIB) $(MISC_LIB) $(LIBSC)

$(O)/getABIcomment: $(OBJSC)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSC) $(LIBS) $(LIBSC)

$(O)/getABIfield: $(OBJSFI)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSFI) $(LIBS) $(LIBSC)

DEPEND_OBJ = $(OBJSSN) $(OBJSST) $(OBJSHX) $(OBJSRW) $(OBJSD) $(OBJSC) $(OBJSF)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
