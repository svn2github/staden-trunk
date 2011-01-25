#
# Makefile for Miscellaneous routines
#
LIBS = misc
PROGS= $(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

CFLAGS += $(SHLIB_CFLAGS) $(MISC_DLL)

INCLUDES_S += $(CORBA_INC)

# pread/pwrite implementations for MinGW/Windows
ifeq ($(findstring MINGW,$(OS)),MINGW)
WIN_FUNCS=win_funcs.o
endif

#
# Objects
#
OBJS = \
	crash.o \
	date.o \
	filenames.o \
	files.o \
	find.o \
	shell.o \
	strings.o \
	usleep.o \
	error.o \
	xerror.o \
	xalloc.o \
	vlen.o \
	FtoC.o \
	getfile.o \
	bitmap.o \
	array.o \
	locks.o \
	array_arith.o \
	parse_db.o\
	dstring.o \
	string_alloc.o\
	$(WIN_FUNCS)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(CORBA_OBJS) $(DEF_FILE) $(L)/.dir
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(CORBA_OBJS) $(CORBA_DEP) $(BIOLIMS_LIB) $(MATH_LIB)

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

#find.c:
#	ln -s ../io_lib/utils/find.c

DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)


# DO NOT DELETE THIS LINE -- make depend depends on it.

FtoC.o: $(SRCROOT)/Misc/os.h
FtoC.o: $(PWD)/staden_config.h
array.o: $(SRCROOT)/Misc/array.h
array.o: $(SRCROOT)/Misc/xalloc.h
array.o: $(SRCROOT)/Misc/xerror.h
array_arith.o: $(SRCROOT)/Misc/misc.h
array_arith.o: $(SRCROOT)/Misc/os.h
array_arith.o: $(SRCROOT)/Misc/xalloc.h
array_arith.o: $(PWD)/staden_config.h
bitmap.o: $(SRCROOT)/Misc/bitmap.h
bitmap.o: $(SRCROOT)/Misc/os.h
bitmap.o: $(SRCROOT)/Misc/xalloc.h
bitmap.o: $(SRCROOT)/Misc/xerror.h
bitmap.o: $(PWD)/staden_config.h
crash.o: $(SRCROOT)/Misc/misc.h
crash.o: $(SRCROOT)/Misc/os.h
crash.o: $(SRCROOT)/Misc/xalloc.h
crash.o: $(PWD)/staden_config.h
date.o: $(SRCROOT)/Misc/misc.h
date.o: $(SRCROOT)/Misc/os.h
date.o: $(SRCROOT)/Misc/xalloc.h
date.o: $(PWD)/staden_config.h
dstring.o: $(SRCROOT)/Misc/dstring.h
dstring.o: $(SRCROOT)/Misc/vlen.h
filenames.o: $(SRCROOT)/Misc/misc.h
filenames.o: $(SRCROOT)/Misc/os.h
filenames.o: $(SRCROOT)/Misc/xalloc.h
filenames.o: $(PWD)/staden_config.h
files.o: $(SRCROOT)/Misc/misc.h
files.o: $(SRCROOT)/Misc/os.h
files.o: $(SRCROOT)/Misc/xalloc.h
files.o: $(PWD)/staden_config.h
find.o: $(SRCROOT)/Misc/misc.h
find.o: $(SRCROOT)/Misc/os.h
find.o: $(SRCROOT)/Misc/xalloc.h
find.o: $(PWD)/staden_config.h
getfile.o: $(SRCROOT)/Misc/getfile.h
getfile.o: $(SRCROOT)/Misc/misc.h
getfile.o: $(SRCROOT)/Misc/os.h
getfile.o: $(SRCROOT)/Misc/xalloc.h
getfile.o: $(PWD)/staden_config.h
locks.o: $(SRCROOT)/Misc/locks.h
parse_db.o: $(SRCROOT)/Misc/misc.h
parse_db.o: $(SRCROOT)/Misc/os.h
parse_db.o: $(SRCROOT)/Misc/parse_db.h
parse_db.o: $(SRCROOT)/Misc/xalloc.h
parse_db.o: $(PWD)/staden_config.h
shell.o: $(SRCROOT)/Misc/misc.h
shell.o: $(SRCROOT)/Misc/os.h
shell.o: $(SRCROOT)/Misc/xalloc.h
shell.o: $(PWD)/staden_config.h
string_alloc.o: $(SRCROOT)/Misc/string_alloc.h
strings.o: $(SRCROOT)/Misc/misc.h
strings.o: $(SRCROOT)/Misc/os.h
strings.o: $(SRCROOT)/Misc/xalloc.h
strings.o: $(PWD)/staden_config.h
vlen.o: $(SRCROOT)/Misc/vlen.h
vlen.o: $(PWD)/staden_config.h
xalloc.o: $(SRCROOT)/Misc/error.h
xerror.o: $(SRCROOT)/Misc/misc.h
xerror.o: $(SRCROOT)/Misc/os.h
xerror.o: $(SRCROOT)/Misc/xalloc.h
xerror.o: $(SRCROOT)/Misc/xerror.h
xerror.o: $(PWD)/staden_config.h
