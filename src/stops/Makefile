#
# Makefile for stops
#

PROGS= 	$(O)/stops

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC) $(TKUTILS_INC)

STOPSSRC=$(O)
COPTDEBUG=$(COPT)

all: $(O)/stops

OBJ=\
	$(STOPSSRC)/stops.o

STOPS_LIBS=\
	$(TEXTUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB)

$(O)/stops: $(OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJ) $(STOPS_LIBS) $(LIBSC) $(MATH_LIB)

DEPEND_OBJ = $(OBJ)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
