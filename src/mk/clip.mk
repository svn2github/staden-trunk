#
# Makefile for clip
#

PROGS= 	$(O)/clip

SRCROOT=..

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC) $(TEXTUTILS_INC)

CLIPSRC=$(O)

all: $(O)/clip

OBJ=\
	$(CLIPSRC)/clip.o\
	$(CLIPSRC)/consen.o

CLIP_LIBS=\
	$(IOLIB_LIB) \
	$(MISC_LIB)

DEPS=\
	$(TEXTUTILS_LIB) \
	$(IOLIB_DEP) \
	$(MISC_DEP)


$(O)/clip: $(OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJ) $(CLIP_LIBS) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJ)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
