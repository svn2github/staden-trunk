#
# Makefile for convert (Alpha version)
#

PROGS = convert

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(GAPDB_EXT_INC) $(TK_INC) $(IOLIB_INC) $(TKUTILS_INC) $(P3_INC) $(TCL_INC) $(GAP4_INC) $(G_INC)

OBJS = \
	list.o \
	main.o \
	dapIO.o \
	dapDB.o \
	bapIO.o \
	bapDB.o \
	gapDB.o \
	process.o \
	flat_sd.o

DEPS = $(GAP_LIB)

convert: $(OBJS)
	$(CLD) -o $@ $(OBJS) $(LIBS) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJS)

include dependencies

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(S)/convert.doc $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)
