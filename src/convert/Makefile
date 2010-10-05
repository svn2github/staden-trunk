#
# Makefile for convert (Alpha version)
#

PROGS = staden_convert

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

DEPS =  $(GAP_LIB) \
	$(G_LIB) \
	$(P3_LIB) \
	$(MUT_LIB) \
	$(TKUTILS_LIB) \
	$(SEQUTILS_LIB) \
	$(MISC_LIB) \
	$(IOLIB_LIB)

staden_convert: $(OBJS)
	$(CLD) -o $@ $(OBJS) $(LIBS) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(S)/convert.doc $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it.

bapDB.o: $(SRCROOT)/Misc/misc.h
bapDB.o: $(SRCROOT)/Misc/os.h
bapDB.o: $(SRCROOT)/Misc/xalloc.h
bapDB.o: $(PWD)/staden_config.h
bapDB.o: $(SRCROOT)/convert/bapIO.h
bapDB.o: $(SRCROOT)/convert/list.h
bapDB.o: $(SRCROOT)/convert/newtypes.h
bapDB.o: $(SRCROOT)/convert/process.h
bapIO.o: $(SRCROOT)/Misc/misc.h
bapIO.o: $(SRCROOT)/Misc/os.h
bapIO.o: $(SRCROOT)/Misc/xalloc.h
bapIO.o: $(PWD)/staden_config.h
bapIO.o: $(SRCROOT)/convert/bapIO.h
bapIO.o: $(SRCROOT)/convert/newtypes.h
dapDB.o: $(SRCROOT)/Misc/misc.h
dapDB.o: $(SRCROOT)/Misc/os.h
dapDB.o: $(SRCROOT)/Misc/xalloc.h
dapDB.o: $(PWD)/staden_config.h
dapDB.o: $(SRCROOT)/convert/dapDB.h
dapDB.o: $(SRCROOT)/convert/dapIO.h
dapDB.o: $(SRCROOT)/convert/list.h
dapDB.o: $(SRCROOT)/convert/newtypes.h
dapDB.o: $(SRCROOT)/convert/process.h
dapIO.o: $(SRCROOT)/Misc/misc.h
dapIO.o: $(SRCROOT)/Misc/os.h
dapIO.o: $(SRCROOT)/Misc/xalloc.h
dapIO.o: $(PWD)/staden_config.h
dapIO.o: $(SRCROOT)/convert/dapIO.h
dapIO.o: $(SRCROOT)/convert/newtypes.h
flat_sd.o: $(SRCROOT)/Misc/misc.h
flat_sd.o: $(SRCROOT)/Misc/os.h
flat_sd.o: $(SRCROOT)/Misc/xalloc.h
flat_sd.o: $(PWD)/staden_config.h
flat_sd.o: $(SRCROOT)/convert/list.h
flat_sd.o: $(SRCROOT)/convert/process.h
gapDB.o: $(SRCROOT)/Misc/array.h
gapDB.o: $(SRCROOT)/Misc/bitmap.h
gapDB.o: $(SRCROOT)/Misc/misc.h
gapDB.o: $(SRCROOT)/Misc/os.h
gapDB.o: $(SRCROOT)/Misc/xalloc.h
gapDB.o: $(SRCROOT)/Misc/xerror.h
gapDB.o: $(PWD)/staden_config.h
gapDB.o: $(SRCROOT)/convert/gapDB.h
gapDB.o: $(SRCROOT)/convert/list.h
gapDB.o: $(SRCROOT)/convert/process.h
gapDB.o: $(SRCROOT)/g/freetree.h
gapDB.o: $(SRCROOT)/g/g-defs.h
gapDB.o: $(SRCROOT)/g/g-error.h
gapDB.o: $(SRCROOT)/g/g-filedefs.h
gapDB.o: $(SRCROOT)/g/g-os.h
gapDB.o: $(SRCROOT)/g/g-struct.h
gapDB.o: $(SRCROOT)/gap4/IO.h
gapDB.o: $(SRCROOT)/gap4/gap-dbstruct.h
gapDB.o: $(SRCROOT)/gap4/gap-if.h
gapDB.o: $(SRCROOT)/gap4/io_handle.h
gapDB.o: $(SRCROOT)/gap4/io_utils.h
list.o: $(SRCROOT)/convert/list.h
main.o: $(SRCROOT)/Misc/array.h
main.o: $(SRCROOT)/Misc/misc.h
main.o: $(SRCROOT)/Misc/os.h
main.o: $(SRCROOT)/Misc/xalloc.h
main.o: $(SRCROOT)/Misc/xerror.h
main.o: $(PWD)/staden_config.h
main.o: $(SRCROOT)/convert/list.h
main.o: $(SRCROOT)/convert/process.h
main.o: $(SRCROOT)/g/freetree.h
main.o: $(SRCROOT)/g/g-error.h
main.o: $(SRCROOT)/g/g-filedefs.h
main.o: $(SRCROOT)/g/g-os.h
main.o: $(SRCROOT)/g/g-struct.h
main.o: $(SRCROOT)/gap4/gap-if.h
process.o: $(SRCROOT)/Misc/misc.h
process.o: $(SRCROOT)/Misc/os.h
process.o: $(SRCROOT)/Misc/xalloc.h
process.o: $(PWD)/staden_config.h
process.o: $(SRCROOT)/convert/bapDB.h
process.o: $(SRCROOT)/convert/dapDB.h
process.o: $(SRCROOT)/convert/flat_sd.h
process.o: $(SRCROOT)/convert/gapDB.h
process.o: $(SRCROOT)/convert/list.h
process.o: $(SRCROOT)/convert/process.h
