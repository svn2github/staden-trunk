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
	dstring.o

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

FtoC.o: $(PWD)/staden_config.h
FtoC.o: os.h
array.o: array.h xerror.h
array.o: xalloc.h
array_arith.o: $(PWD)/staden_config.h
array_arith.o: misc.h
array_arith.o: os.h
array_arith.o: xalloc.h
bitmap.o: $(PWD)/staden_config.h
bitmap.o: bitmap.h os.h
bitmap.o: xerror.h xalloc.h
crash.o: $(PWD)/staden_config.h
crash.o: misc.h
crash.o: os.h
crash.o: xalloc.h
date.o: $(PWD)/staden_config.h
date.o: misc.h
date.o: os.h
date.o: xalloc.h
dstring.o: $(SRCROOT)/Misc/dstring.h
dstring.o: vlen.h
filenames.o: $(PWD)/staden_config.h
filenames.o: misc.h
filenames.o: os.h
filenames.o: xalloc.h
files.o: $(PWD)/staden_config.h
files.o: misc.h
files.o: os.h
files.o: xalloc.h
find.o: $(PWD)/staden_config.h
find.o: misc.h
find.o: os.h
find.o: xalloc.h
getfile.o: $(PWD)/staden_config.h
getfile.o: misc.h os.h
getfile.o: xalloc.h getfile.h
locks.o: locks.h
parse_db.o: $(PWD)/staden_config.h
parse_db.o: misc.h os.h
parse_db.o: xalloc.h parse_db.h
shell.o: $(PWD)/staden_config.h
shell.o: misc.h
shell.o: os.h
shell.o: xalloc.h
strings.o: $(PWD)/staden_config.h
strings.o: misc.h
strings.o: os.h
strings.o: xalloc.h
vlen.o: $(PWD)/staden_config.h
vlen.o: vlen.h
xalloc.o: error.h
xerror.o: $(PWD)/staden_config.h
xerror.o: os.h
xerror.o: xalloc.h
xerror.o: xerror.h misc.h
