#
# Makefile for clip
#

PROGS= 	$(O)/qclip

SRCROOT=..

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC) $(TKUTILS_INC)

QCLIPSRC=$(O)

all: $(O)/qclip

OBJ=\
	$(QCLIPSRC)/qclip.o\
	$(QCLIPSRC)/consen.o

QCLIP_LIBS=\
	$(TEXTUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB)

DEPS=\
	$(IOLIB_DEP) \
	$(MISC_DEP)


$(O)/qclip: $(OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJ) $(QCLIP_LIBS) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJ)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
