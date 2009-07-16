#
# Makefile for convert (Alpha version)
#

PROGS = $(O)/convert

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include $(SRCROOT)/gap4/gap4_defs.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(GAPDB_EXT_INC) $(TK_INC) $(IOLIB_INC) $(TKUTILS_INC) $(P3_INC) $(TCL_INC)
CONVERTBIN=$(O)

OBJS = \
	$(CONVERTBIN)/list.o \
	$(CONVERTBIN)/main.o \
	$(CONVERTBIN)/dapIO.o \
	$(CONVERTBIN)/dapDB.o \
	$(CONVERTBIN)/bapIO.o \
	$(CONVERTBIN)/bapDB.o \
	$(CONVERTBIN)/gapDB.o \
	$(CONVERTBIN)/process.o \
	$(CONVERTBIN)/flat_sd.o \
	$(GAPDB_EXT_OBJS)

DEPS = $(GAPDB_EXT_DEPS)

LIBS = $(GAPDB_EXT_LIBS)

$(O)/convert: $(OBJS)
	$(CLD) -o $@ $(OBJS) $(LIBS) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJS)

include dependencies

distsrc: distsrc_dirs
	-cp -R *.[ch] Makefile dependencies convert.doc $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)
