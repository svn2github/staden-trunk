#
# Makefile for tk_utils routines
#

LIBS = tk_utils
ifneq ($(MACHINE),windows)
PROGS= $(LIBS)
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
#DEFINES += -DTCL_LOCAL_APPINIT=Stash_AppInit

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
	element_canvas.o\
	container.o\
	container_ruler.o\
	$(TKUTILS_EXTRAS)

# Spin2 only objects, which currently isn't compiled and has compilation
# issues on some systems.
#
#	tkCanvGraph.o\
#	container_cmds.o\
#	tclCanvGraph.o\
#	matrix.o

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

TKU_LIBS = \
	$(IOLIB_LIB) \
	$(TK_LIB) \
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
	cp $(S)/stash $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)
	-mkdir $(INSTALLTCL)/tk_utils
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/tk_utils
	cp $(S)/stash $(INSTALLBIN)
ifneq ($(MACHINE),windows)
endif
	cp $(S)/2bars_*.bmp $(S)/help_config $(S)/tk_utilsrc $(INSTALLETC)

DEPEND_OBJ = $(OBJS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

canvas_box.o: $(SRCROOT)/Misc/misc.h
canvas_box.o: $(SRCROOT)/Misc/os.h
canvas_box.o: $(SRCROOT)/Misc/xalloc.h
canvas_box.o: $(PWD)/staden_config.h
canvas_box.o: $(SRCROOT)/tk_utils/canvas_box.h
canvas_box.o: $(SRCROOT)/tk_utils/container.h
canvas_box.o: $(SRCROOT)/tk_utils/element_canvas.h
canvas_box.o: $(SRCROOT)/tk_utils/ruler_tick.h
canvas_box.o: $(SRCROOT)/tk_utils/tclCanvGraph.h
canvas_box.o: $(SRCROOT)/tk_utils/tclXkeylist.h
canvas_box.o: $(SRCROOT)/tk_utils/tcl_utils.h
capture.o: $(SRCROOT)/Misc/misc.h
capture.o: $(SRCROOT)/Misc/os.h
capture.o: $(SRCROOT)/Misc/xalloc.h
capture.o: $(PWD)/staden_config.h
capture.o: $(SRCROOT)/tk_utils/capture.h
capture.o: $(SRCROOT)/tk_utils/tcl_utils.h
cli_arg.o: $(SRCROOT)/Misc/misc.h
cli_arg.o: $(SRCROOT)/Misc/os.h
cli_arg.o: $(SRCROOT)/Misc/xalloc.h
cli_arg.o: $(PWD)/staden_config.h
cli_arg.o: $(SRCROOT)/tk_utils/cli_arg.h
container.o: $(SRCROOT)/Misc/xalloc.h
container.o: $(SRCROOT)/tk_utils/canvas_box.h
container.o: $(SRCROOT)/tk_utils/container.h
container.o: $(SRCROOT)/tk_utils/container_ruler.h
container.o: $(SRCROOT)/tk_utils/element_canvas.h
container.o: $(SRCROOT)/tk_utils/tclCanvGraph.h
container.o: $(SRCROOT)/tk_utils/tcl_utils.h
container.o: $(SRCROOT)/tk_utils/text_output.h
container_ruler.o: $(SRCROOT)/Misc/xalloc.h
container_ruler.o: $(SRCROOT)/tk_utils/canvas_box.h
container_ruler.o: $(SRCROOT)/tk_utils/container.h
container_ruler.o: $(SRCROOT)/tk_utils/tcl_utils.h
container_ruler.o: $(SRCROOT)/tk_utils/text_output.h
element_canvas.o: $(SRCROOT)/tk_utils/canvas_box.h
element_canvas.o: $(SRCROOT)/tk_utils/container.h
element_canvas.o: $(SRCROOT)/tk_utils/element_canvas.h
element_canvas.o: $(SRCROOT)/tk_utils/tclCanvGraph.h
element_canvas.o: $(SRCROOT)/tk_utils/text_output.h
init.o: $(SRCROOT)/tk_utils/tclCanvGraph.h
init.o: $(SRCROOT)/tk_utils/tclXkeylist.h
init.o: $(SRCROOT)/tk_utils/user_defaults.h
misc.o: $(SRCROOT)/Misc/getfile.h
misc.o: $(SRCROOT)/Misc/misc.h
misc.o: $(SRCROOT)/Misc/os.h
misc.o: $(SRCROOT)/Misc/xalloc.h
misc.o: $(PWD)/staden_config.h
misc.o: $(SRCROOT)/seq_utils/renz_utils.h
misc.o: $(SRCROOT)/tk_utils/capture.h
misc.o: $(SRCROOT)/tk_utils/cli_arg.h
misc.o: $(SRCROOT)/tk_utils/tcl_utils.h
postscript.o: $(SRCROOT)/Misc/misc.h
postscript.o: $(SRCROOT)/Misc/os.h
postscript.o: $(SRCROOT)/Misc/xalloc.h
postscript.o: $(PWD)/staden_config.h
postscript.o: $(SRCROOT)/tk_utils/cli_arg.h
postscript.o: $(SRCROOT)/tk_utils/postscript.h
postscript.o: $(SRCROOT)/tk_utils/split.h
restriction_enzyme_map.o: $(SRCROOT)/Misc/misc.h
restriction_enzyme_map.o: $(SRCROOT)/Misc/os.h
restriction_enzyme_map.o: $(SRCROOT)/Misc/xalloc.h
restriction_enzyme_map.o: $(PWD)/staden_config.h
restriction_enzyme_map.o: $(SRCROOT)/seq_utils/renz_utils.h
restriction_enzyme_map.o: $(SRCROOT)/tk_utils/canvas_box.h
restriction_enzyme_map.o: $(SRCROOT)/tk_utils/tcl_utils.h
restriction_enzyme_map.o: $(SRCROOT)/tk_utils/text_output.h
ruler_tick.o: $(SRCROOT)/Misc/misc.h
ruler_tick.o: $(SRCROOT)/Misc/os.h
ruler_tick.o: $(SRCROOT)/Misc/xalloc.h
ruler_tick.o: $(PWD)/staden_config.h
ruler_tick.o: $(SRCROOT)/tk_utils/canvas_box.h
ruler_tick.o: $(SRCROOT)/tk_utils/container.h
ruler_tick.o: $(SRCROOT)/tk_utils/ruler_tick.h
sheet.o: $(SRCROOT)/Misc/misc.h
sheet.o: $(SRCROOT)/Misc/os.h
sheet.o: $(SRCROOT)/Misc/xalloc.h
sheet.o: $(PWD)/staden_config.h
sheet.o: $(SRCROOT)/tk_utils/sheet.h
split.o: $(SRCROOT)/Misc/misc.h
split.o: $(SRCROOT)/Misc/os.h
split.o: $(SRCROOT)/Misc/xalloc.h
split.o: $(PWD)/staden_config.h
split.o: $(SRCROOT)/tk_utils/split.h
tclXkeylist.o: $(SRCROOT)/tk_utils/tclExtdInt.h
tclXutil.o: $(SRCROOT)/tk_utils/tclExtdInt.h
tcl_debug.o: $(SRCROOT)/tk_utils/tcl_utils.h
tcl_utils.o: $(SRCROOT)/Misc/misc.h
tcl_utils.o: $(SRCROOT)/Misc/os.h
tcl_utils.o: $(SRCROOT)/Misc/vlen.h
tcl_utils.o: $(SRCROOT)/Misc/xalloc.h
tcl_utils.o: $(PWD)/staden_config.h
tcl_utils.o: $(SRCROOT)/tk_utils/tclXkeylist.h
tcl_utils.o: $(SRCROOT)/tk_utils/tcl_utils.h
tcl_utils.o: $(SRCROOT)/tk_utils/text_output.h
text_output.o: $(SRCROOT)/Misc/FtoC.h
text_output.o: $(SRCROOT)/Misc/getfile.h
text_output.o: $(SRCROOT)/Misc/misc.h
text_output.o: $(SRCROOT)/Misc/os.h
text_output.o: $(SRCROOT)/Misc/vlen.h
text_output.o: $(SRCROOT)/Misc/xalloc.h
text_output.o: $(PWD)/staden_config.h
text_output.o: $(SRCROOT)/tk_utils/tcl_utils.h
text_output.o: $(SRCROOT)/tk_utils/text_output.h
tkRaster.o: $(SRCROOT)/Misc/xalloc.h
tkRaster.o: $(SRCROOT)/tk_utils/tcl_utils.h
tkRaster.o: $(SRCROOT)/tk_utils/tkRaster.h
tkRaster.o: $(SRCROOT)/tk_utils/tkRasterBuiltIn.h
tkRasterBuiltIn.o: $(SRCROOT)/tk_utils/tkRaster.h
tkRasterBuiltIn.o: $(SRCROOT)/tk_utils/tkRasterBuiltIn.h
tkSheet.o: $(SRCROOT)/Misc/xalloc.h
tkSheet.o: $(SRCROOT)/tk_utils/sheet.h
tkSheet.o: $(SRCROOT)/tk_utils/tcl_utils.h
tkSheet.o: $(SRCROOT)/tk_utils/tkSheet.h
tkSheet.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tkSheet.o: $(SRCROOT)/tk_utils/tkSheet_config.h
tkSheet.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tkSheet.o: $(SRCROOT)/tk_utils/tk_defs.h
tkSheet_common.o: $(SRCROOT)/Misc/xalloc.h
tkSheet_common.o: $(SRCROOT)/tk_utils/sheet.h
tkSheet_common.o: $(SRCROOT)/tk_utils/tkSheet.h
tkSheet_common.o: $(SRCROOT)/tk_utils/tkSheet_common.h
tkSheet_common.o: $(SRCROOT)/tk_utils/tkSheet_struct.h
tkSheet_common.o: $(SRCROOT)/tk_utils/tk_defs.h
tkTrace.o: $(SRCROOT)/Misc/misc.h
tkTrace.o: $(SRCROOT)/Misc/os.h
tkTrace.o: $(SRCROOT)/Misc/xalloc.h
tkTrace.o: $(PWD)/staden_config.h
tkTrace.o: $(SRCROOT)/tk_utils/cli_arg.h
tkTrace.o: $(SRCROOT)/tk_utils/postscript.h
tkTrace.o: $(SRCROOT)/tk_utils/split.h
tkTrace.o: $(SRCROOT)/tk_utils/tcl_utils.h
tkTrace.o: $(SRCROOT)/tk_utils/tkTrace.h
tkTrace.o: $(SRCROOT)/tk_utils/tkTraceIO.h
tkTrace.o: $(SRCROOT)/tk_utils/tk_defs.h
tkTrace.o: $(SRCROOT)/tk_utils/trace_print.h
tkTraceComp.o: $(SRCROOT)/tk_utils/cli_arg.h
tkTraceComp.o: $(SRCROOT)/tk_utils/postscript.h
tkTraceComp.o: $(SRCROOT)/tk_utils/tkTrace.h
tkTraceComp.o: $(SRCROOT)/tk_utils/tkTraceIO.h
tkTraceDisp.o: $(SRCROOT)/Misc/misc.h
tkTraceDisp.o: $(SRCROOT)/Misc/os.h
tkTraceDisp.o: $(SRCROOT)/Misc/xalloc.h
tkTraceDisp.o: $(PWD)/staden_config.h
tkTraceDisp.o: $(SRCROOT)/tk_utils/cli_arg.h
tkTraceDisp.o: $(SRCROOT)/tk_utils/postscript.h
tkTraceDisp.o: $(SRCROOT)/tk_utils/tkTrace.h
tkTraceDisp.o: $(SRCROOT)/tk_utils/tkTraceIO.h
tkTraceIO.o: $(SRCROOT)/Misc/misc.h
tkTraceIO.o: $(SRCROOT)/Misc/os.h
tkTraceIO.o: $(SRCROOT)/Misc/xalloc.h
tkTraceIO.o: $(PWD)/staden_config.h
tkTraceIO.o: $(SRCROOT)/tk_utils/cli_arg.h
tkTraceIO.o: $(SRCROOT)/tk_utils/postscript.h
tkTraceIO.o: $(SRCROOT)/tk_utils/tkTrace.h
tkTraceIO.o: $(SRCROOT)/tk_utils/trace_print.h
trace_print.o: $(SRCROOT)/Misc/misc.h
trace_print.o: $(SRCROOT)/Misc/os.h
trace_print.o: $(SRCROOT)/Misc/xalloc.h
trace_print.o: $(PWD)/staden_config.h
trace_print.o: $(SRCROOT)/tk_utils/cli_arg.h
trace_print.o: $(SRCROOT)/tk_utils/postscript.h
trace_print.o: $(SRCROOT)/tk_utils/split.h
trace_print.o: $(SRCROOT)/tk_utils/tkTrace.h
trace_print.o: $(SRCROOT)/tk_utils/trace_print.h
