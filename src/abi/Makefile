#
# Makefile for getABISampleName
#

PROGS = getABISampleName getABIcomment getABIdate getABIfield

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC)

OBJSSN = abiIO.o getABISampleName.o
OBJSD  = abiIO.o getABIdate.o
OBJSFI = abiIO.o getABIfield.o
OBJSC  = abiIO.o getABIcomment.o

LIBS=\
	$(IOLIB_LIB) \
	$(MISC_LIB)

getABISampleName: $(OBJSSN)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSSN) $(LIBS) $(LIBSC)

getABIdate: $(OBJSD)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(IOUTILS_LIB) $(MISC_LIB) $(LIBSC)

getABIcomment: $(OBJSC)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSC) $(LIBS) $(LIBSC)

getABIfield: $(OBJSFI)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSFI) $(LIBS) $(LIBSC)

DEPEND_OBJ = $(OBJSSN) $(OBJSD) $(OBJSC) $(OBJSFI)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
