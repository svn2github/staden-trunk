#
# Makefile for get_scf_field
#

PROGS= 	$(O)/get_scf_field

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(IOLIB_INC)

SCFFIELDSRC=$(O)

all: $(O)/get_scf_field

OBJS	= $(SCFFIELDSRC)/get_scf_field.o

GSF_LIBS=\
	$(IOLIB_LIB) \
	$(IOUTILS_LIB) \
	$(MISC_LIB)

$(O)/get_scf_field:	$(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(IOLIB_LIB) $(LIBSC) $(GSF_LIBS)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
