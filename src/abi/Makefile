#
# Makefile for getABISampleName
#

PROGS = $(O)/getABISampleName $(O)/getABIcomment $(O)/getABIdate \
	$(O)/getABIstring $(O)/getABIhex $(O)/getABIraw

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC)
ABIBIN=$(O)

OBJSSN = $(ABIBIN)/getABISampleName.o
OBJSD  = $(ABIBIN)/getABIdate.o
OBJSHX = $(ABIBIN)/getABIhex.o
OBJSRW = $(ABIBIN)/getABIraw.o
OBJSST = $(ABIBIN)/getABIstring.o
OBJSC  = $(ABIBIN)/getABIcomment.o

LIBS=\
	$(IOLIB_LIB) \
	$(MISC_LIB)

$(O)/getABISampleName: $(OBJSSN)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSSN) $(LIBS) $(LIBSC)

$(O)/getABIdate: $(OBJSD)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(IOUTILS_LIB) $(MISC_LIB) $(LIBSC)

$(O)/getABIcomment: $(OBJSC)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSC) $(LIBS) $(LIBSC)

$(O)/getABIstring: $(OBJSST)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSST) $(LIBS) $(LIBSC)

$(O)/getABIhex: $(OBJSHX)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSHX) $(LIBS) $(LIBSC)

$(O)/getABIraw: $(OBJSRW)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSRW) $(LIBS) $(LIBSC)

DEPEND_OBJ = $(OBJSSN) $(OBJSST) $(OBJSHX) $(OBJSRW) $(OBJSD) $(OBJSC)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
