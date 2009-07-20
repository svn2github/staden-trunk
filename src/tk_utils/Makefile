#
# Makefile for tk_utils routines
#

LIBS = tk_utils
ifneq ($(MACHINE),windows)
PROGS= $(LIBS) stash
else
PROGS= $(LIBS)
endif
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(TKUTILS_INC) $(IOLIB_INC) $(TK_INC) \
	      $(SEQUTILS_INC) $(BIOLIMS_INC) $(PNG_INC)

CFLAGS += $(SHLIB_CFLAGS) $(TK_UTILS_DLL)

#DEFINES += -DCHECK_LICENCE -DTCL_LOCAL_APPINIT=Stash_AppInit
DEFINES += -DTCL_LOCAL_APPINIT=Stash_AppInit

#
# Objects
#
OBJS = \
	cli_arg.o \
	tclXkeylist.o \
	tclXutil.o \
	tcl_utils.o \
	tcl_debug.o \
	misc.o \
	init.o \
	text_output.o \
	tkRaster.o \
	tkRasterBuiltIn.o \
	sheet.o \
	tkSheet.o \
	tkSheet_common.o \
	trace_print.o \
	postscript.o \
	split.o \
	tkTrace.o \
	tkTraceComp.o \
	tkTraceIO.o \
	tkTraceDisp.o \
	capture.o \
	canvas_box.o \
	ruler_tick.o \
	restriction_enzyme_map.o \
	container_cmds.o\
	container.o\
	container_ruler.o\
	element_canvas.o\
	tclCanvGraph.o\
	tkCanvGraph.o\
	matrix.o\
	$(TKUTILS_EXTRAS)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

TKU_LIBS = \
	$(IOLIB_LIB) \
	$(TK_LIB) \
	$(TCL_LIB) \
	$(MISC_LIB) \
	$(PNG_LIB) \
	$(ZLIB_LIB) \
	$(SOCKET)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(TKU_LIBS)


$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

# Stash interpreter executable
# Uses the main tcl interpreter, but overrides the appinit to use a local
# function.
ifneq ($(MACHINE),windows)
STASH_OBJS = \
	stash.o \
	tclAppInit.o

STASH_LIBS = \
	$(TKUTILS_LIB) \
	$(TK_LIB)

stash: $(STASH_OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMWIN) $(STASH_OBJS) $(STASH_LIBS) $(LIBSC)
endif

# Copy the tk_utils source
distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/tclIndex $(S)/Makefile $(DIRNAME)
	cp $(S)/*.bmp $(S)/help_config $(S)/tk_utilsrc $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)
	-mkdir $(INSTALLTCL)/tk_utils
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/tk_utils
ifneq ($(MACHINE),windows)
	cp stash $(INSTALLBIN)
endif
	cp $(S)/2bars_*.bmp $(S)/help_config $(S)/tk_utilsrc $(INSTALLETC)

include dependencies
# DO NOT DELETE
