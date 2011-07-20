LIBS = primer3
PROGS = $(LIBS) primer3_core
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/../..
include $(SRCROOT)/global.mk
include ../../system.mk

CFLAGS += $(SHLIB_CFLAGS)

MAX_PRIMER_LENGTH = 36
P_DEFINES = -DDPAL_MAX_ALIGN=$(MAX_PRIMER_LENGTH) -DMAX_PRIMER_LENGTH=$(MAX_PRIMER_LENGTH)
DEFINES += $(P_DEFINES)

LIBOBJS=\
	primer3_lib.o\
	primer3_release.o\
	oligotm.o\
	dpal_primer.o\
	format_output.o\
	boulder_input.o

P3OBJS=primer3_main.o
P3LIBS=$(P3_LIB) $(MATH_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(LIBOBJS)
	$(MKDEFL) $@ $(LIBOBJS)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(LIBOBJS) $(DEF_FILE) $(L)/.dir
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LIBOBJS) $(SHLIB_DEP) $(MATH_LIB)

primer3_core: $(P3OBJS) $(PROGLIBS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(P3OBJS) $(P3LIBS) $(LIBSC)

dpal_primer.o: dpal.c dpal.h primer3_release.h
	$(CC) $(CFLAGS) $(P_DEFINES) $(COBJFLAG)$@ -c $<
	$(MKDEFC) $(MKFLAGS) $@

install:
	cp $(PROGLIBS) $(INSTALLLIB)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

# DO NOT DELETE THIS LINE -- make depend depends on it
