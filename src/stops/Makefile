#
# Makefile for stops
#

PROGS= 	stops

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC) $(TKUTILS_INC)

OBJ=	stops.o

STOPS_LIBS=\
	$(TEXTUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB)

stops: $(OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJ) $(STOPS_LIBS) $(LIBSC) $(MATH_LIB)

DEPEND_OBJ = $(OBJ)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
