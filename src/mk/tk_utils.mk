#
# Makefile for tk_utils routines
#

LIBS = tk_utils
PROGS= $(LIBS) $(O)/stash $(WINSTASH)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(TKUTILS_INC) $(IOLIB_INC) $(TK_INC) $(SEQUTILS_INC) \
	$(BIOLIMS_INC)

CFLAGS += $(SHLIB_CFLAGS) $(TK_UTILS_DLL)

TKUTILSBIN=$(O)

DEFINES += -DCHECK_LICENCE

#
# Objects
#
OBJS = \
	$(TKUTILSBIN)/cli_arg.o \
	$(TKUTILSBIN)/tclXkeylist.o \
	$(TKUTILSBIN)/tclXutil.o \
	$(TKUTILSBIN)/tcl_utils.o \
	$(TKUTILSBIN)/tcl_debug.o \
	$(TKUTILSBIN)/misc.o \
	$(TKUTILSBIN)/init.o \
	$(TKUTILSBIN)/text_output.o \
	$(TKUTILSBIN)/tkRaster.o \
	$(TKUTILSBIN)/tkRasterBuiltIn.o \
	$(TKUTILSBIN)/sheet.o \
	$(TKUTILSBIN)/tkSheet.o \
	$(TKUTILSBIN)/tkSheet_common.o \
	$(TKUTILSBIN)/trace_print.o \
	$(TKUTILSBIN)/postscript.o \
	$(TKUTILSBIN)/split.o \
	$(TKUTILSBIN)/tkTrace.o \
	$(TKUTILSBIN)/tkTraceComp.o \
	$(TKUTILSBIN)/tkTraceIO.o \
	$(TKUTILSBIN)/tkTraceDisp.o \
	$(TKUTILSBIN)/capture.o \
	$(TKUTILSBIN)/canvas_box.o \
	$(TKUTILSBIN)/ruler_tick.o \
	$(TKUTILSBIN)/restriction_enzyme_map.o \
	$(TKUTILSBIN)/check_licence.o\
	$(TKUTILSBIN)/licence_utils.o\
	$(TKUTILSBIN)/valid_seq.o\
	$(TKUTILSBIN)/md52.o\
	$(TKUTILS_EXTRAS)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

DEPS = \
	$(IOLIB_DEP) \
	$(TK_DEP) \
	$(TCL_DEP) \
	$(MISC_DEP) \
	$(SOCKET)

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(DEPS)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for WINNT
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


#
# The Staden program shell - stash
#
STASH_OBJS = \
	$(TKUTILSBIN)/tkMain.o\
	$(TKUTILSBIN)/stash.o\
	$(TKUTILSBIN)/check_licence.o\
	$(TKUTILSBIN)/licence_utils.o\
	$(TKUTILSBIN)/valid_seq.o\
	$(TKUTILSBIN)/md52.o

STASHDEP = \
	$(IOLIB_DEP) \
	$(TKUTILS_DEP) \
	$(TK_DEP) \
	$(TCL_DEP) \
	$(MISC_DEP)

STASH_LIBS = \
	$(IOLIB_LIB) \
	$(TKUTILS_LIB) \
	$(TK_LIB) \
	$(TCL_LIB) \
	$(MISC_LIB)

# 7/1/99 johnt - added abstractions to support Visual C++ linker
$(O)/stash: $(STASH_OBJS) $(WINSTASH)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(STASH_OBJS) $(STASH_LIBS) $(LIBSC) $(STASHDEP)

# 7/1/99 johnt - added additional target winstash to allow stash to be run without windows console displayed
#              - need to maintain console version (stash) to allow non tk scripts out output to stdout
$(WINSTASH): $(STASH_OBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMWIN) $(STASH_OBJS) $(STASH_LIBS) $(LIBSC) $(STASHDEP)

DEPEND_OBJ = $(OBJS) $(STASH_OBJS)


# Copy the tk_utils source, but also copy some dummy licence files.
distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl tclIndex Makefile dependencies \
		$(DIRNAME)
	-rm $(DIRNAME)/resource.h
	rm $(DIRNAME)/check_licence.c \
	   $(DIRNAME)/licence_utils.c \
	   $(DIRNAME)/licence_utils.h \
	   $(DIRNAME)/licence.h \
	   $(DIRNAME)/valid_seq.c \
	   $(DIRNAME)/md52.c \
	   $(DIRNAME)/md52.h \
	   $(DIRNAME)/boxes.h
	-cp distrib/check_licence.c \
	    distrib/licence_utils.c \
	    distrib/licence_utils.h \
	    distrib/licence.h \
	    distrib/valid_seq.c \
	    distrib/md52.c \
	    distrib/md52.h \
	    distrib/boxes.h \
	    $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)
	-mkdir $(INSTALLLIB)/tk_utils
	cp *.tcl tclIndex $(INSTALLLIB)/tk_utils
	cp $(O)/stash $(INSTALLBIN)

include dependencies
# DO NOT DELETE
