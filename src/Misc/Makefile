#
# Makefile for Miscellaneous routines
#
LIBS = misc
PROGS= $(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

CFLAGS += $(SHLIB_CFLAGS) $(MISC_DLL)

INCLUDES_S += $(CORBA_INC)

MISCBIN=$(O)

#
# Objects
#
OBJS = \
	$(MISCBIN)/crash.o \
	$(MISCBIN)/date.o \
	$(MISCBIN)/filenames.o \
	$(MISCBIN)/files.o \
	$(MISCBIN)/find.o \
	$(MISCBIN)/shell.o \
	$(MISCBIN)/strings.o \
	$(MISCBIN)/usleep.o \
	$(MISCBIN)/error.o \
	$(MISCBIN)/xerror.o \
	$(MISCBIN)/xalloc.o \
	$(MISCBIN)/vlen.o \
	$(MISCBIN)/FtoC.o \
	$(MISCBIN)/getfile.o \
	$(MISCBIN)/bitmap.o \
	$(MISCBIN)/array.o \
	$(MISCBIN)/locks.o \
	$(MISCBIN)/array_arith.o \
	$(MISCBIN)/parse_db.o\
	$(MISCBIN)/dstring.o

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(CORBA_OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(CORBA_OBJS) $(CORBA_DEP) $(BIOLIMS_LIB) $(MATH_LIB)

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

find.c:
	ln -s ../io_lib/utils/find.c

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

include dependencies
