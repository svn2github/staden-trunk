#
# Makefile for get_scf_field
#

PROGS= 	get_scf_field

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(IOLIB_INC)

OBJS	= get_scf_field.o

GSF_LIBS=\
	$(IOLIB_LIB) \
	$(IOUTILS_LIB) \
	$(MISC_LIB)

get_scf_field:	$(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(IOLIB_LIB) $(LIBSC) $(GSF_LIBS)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
