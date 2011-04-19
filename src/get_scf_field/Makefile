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

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it.

get_scf_field.o: $(SRCROOT)/Misc/misc.h
get_scf_field.o: $(SRCROOT)/Misc/os.h
get_scf_field.o: $(SRCROOT)/Misc/xalloc.h
get_scf_field.o: $(PWD)/staden_config.h
