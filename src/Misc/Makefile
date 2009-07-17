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

install:
	cp $(PROGLIBS) $(INSTALLLIB)


include dependencies
