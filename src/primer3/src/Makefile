LIBS = primer3
PROGS = $(LIBS) $(O)/primer3_core
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=../..

include	$(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

P3BIN=$(O)
#L=$(O)

MAX_PRIMER_LENGTH = 36
P_DEFINES = -DDPAL_MAX_ALIGN=$(MAX_PRIMER_LENGTH) -DMAX_PRIMER_LENGTH=$(MAX_PRIMER_LENGTH)
DEFINES += $(P_DEFINES)

LIBOBJS=\
	$(P3BIN)/primer3_lib.o\
	$(P3BIN)/primer3_release.o\
	$(P3BIN)/oligotm.o\
	$(P3BIN)/dpal_primer.o\
	$(P3BIN)/format_output.o\
	$(P3BIN)/boulder_input.o

P3OBJS=$(P3BIN)/primer3_main.o
P3LIBS=$(P3_LIB) $(MATH_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(LIBOBJS)
	$(MKDEFL) $@ $(LIBOBJS)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(LIBOBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LIBOBJS) $(SHLIB_DEP)

$(O)/primer3_core: $(P3OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(P3OBJS) $(P3LIBS) $(LIBSC)

$(O)/dpal_primer.o: dpal.c dpal.h primer3_release.h
	$(CC) $(CFLAGS) $(P_DEFINES) $(COBJFLAG)$@ -c $<
	$(MKDEFC) $(MKFLAGS) $@

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

distsrc: distsrc_dirs
	-cp -R *.[ch] Makefile dependencies \
		$(DIRNAME)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
