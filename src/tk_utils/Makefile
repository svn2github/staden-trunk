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
endif
	cp $(S)/2bars_*.bmp $(S)/help_config $(S)/tk_utilsrc $(INSTALLETC)

DEPEND_OBJ = $(OBJS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

canvas_box.o: $(SRCROOT)/Misc/misc.h
canvas_box.o: $(SRCROOT)/Misc/os.h
canvas_box.o: $(SRCROOT)/Misc/xalloc.h
canvas_box.o: $(PWD)/staden_config.h
canvas_box.o: canvas_box.h
canvas_box.o: container.h
canvas_box.o: element_canvas.h
canvas_box.o: ruler_tick.h
canvas_box.o: tclCanvGraph.h
canvas_box.o: tclXkeylist.h
canvas_box.o: tcl_utils.h
capture.o: $(SRCROOT)/Misc/misc.h
capture.o: $(SRCROOT)/Misc/os.h
capture.o: $(SRCROOT)/Misc/xalloc.h
capture.o: $(PWD)/staden_config.h
capture.o: capture.h
capture.o: tcl_utils.h
cli_arg.o: $(SRCROOT)/Misc/misc.h
cli_arg.o: $(SRCROOT)/Misc/os.h
cli_arg.o: $(SRCROOT)/Misc/xalloc.h
cli_arg.o: $(PWD)/staden_config.h
cli_arg.o: cli_arg.h
container.o: $(SRCROOT)/Misc/xalloc.h
container.o: canvas_box.h
container.o: container.h
container.o: container_ruler.h
container.o: element_canvas.h
container.o: tclCanvGraph.h
container.o: tcl_utils.h
container.o: text_output.h
container_cmds.o: $(SRCROOT)/Misc/xalloc.h
container_cmds.o: canvas_box.h
container_cmds.o: cli_arg.h
container_cmds.o: container.h
container_cmds.o: tclXkeylist.h
container_cmds.o: tcl_utils.h
container_cmds.o: text_output.h
container_ruler.o: $(SRCROOT)/Misc/xalloc.h
container_ruler.o: canvas_box.h
container_ruler.o: container.h
container_ruler.o: tcl_utils.h
container_ruler.o: text_output.h
element_canvas.o: canvas_box.h
element_canvas.o: container.h
element_canvas.o: element_canvas.h
element_canvas.o: tclCanvGraph.h
element_canvas.o: text_output.h
init.o: tclCanvGraph.h
init.o: tclXkeylist.h
init.o: tkCanvGraph.h
init.o: user_defaults.h
matrix.o: matrix.h
misc.o: $(SRCROOT)/Misc/getfile.h
misc.o: $(SRCROOT)/Misc/misc.h
misc.o: $(SRCROOT)/Misc/os.h
misc.o: $(SRCROOT)/Misc/xalloc.h
misc.o: $(SRCROOT)/seq_utils/renz_utils.h
misc.o: $(PWD)/staden_config.h
misc.o: cli_arg.h
misc.o: tcl_utils.h capture.h
postscript.o: $(SRCROOT)/Misc/misc.h
postscript.o: $(SRCROOT)/Misc/os.h
postscript.o: $(SRCROOT)/Misc/xalloc.h
postscript.o: $(PWD)/staden_config.h
postscript.o: cli_arg.h
postscript.o: postscript.h
postscript.o: split.h
restriction_enzyme_map.o: $(SRCROOT)/Misc/misc.h
restriction_enzyme_map.o: $(SRCROOT)/Misc/os.h
restriction_enzyme_map.o: $(SRCROOT)/Misc/xalloc.h
restriction_enzyme_map.o: $(SRCROOT)/seq_utils/renz_utils.h
restriction_enzyme_map.o: $(PWD)/staden_config.h
restriction_enzyme_map.o: canvas_box.h
restriction_enzyme_map.o: tcl_utils.h
restriction_enzyme_map.o: text_output.h
ruler_tick.o: $(SRCROOT)/Misc/misc.h
ruler_tick.o: $(SRCROOT)/Misc/os.h
ruler_tick.o: $(SRCROOT)/Misc/xalloc.h
ruler_tick.o: $(PWD)/staden_config.h
ruler_tick.o: canvas_box.h
ruler_tick.o: container.h
ruler_tick.o: ruler_tick.h
sheet.o: $(SRCROOT)/Misc/misc.h
sheet.o: $(SRCROOT)/Misc/os.h
sheet.o: $(SRCROOT)/Misc/xalloc.h
sheet.o: $(PWD)/staden_config.h
sheet.o: sheet.h
split.o: $(SRCROOT)/Misc/misc.h
split.o: $(SRCROOT)/Misc/os.h
split.o: $(SRCROOT)/Misc/xalloc.h
split.o: $(PWD)/staden_config.h
split.o: split.h
stash.o: tclCanvGraph.h
stash.o: tkCanvGraph.h
tclCanvGraph.o: tclCanvGraph.h
tclXkeylist.o: tclExtdInt.h
tclXutil.o: tclExtdInt.h
tcl_debug.o: tcl_utils.h
tcl_utils.o: $(SRCROOT)/Misc/misc.h
tcl_utils.o: $(SRCROOT)/Misc/os.h
tcl_utils.o: $(SRCROOT)/Misc/vlen.h
tcl_utils.o: $(SRCROOT)/Misc/xalloc.h
tcl_utils.o: $(PWD)/staden_config.h
tcl_utils.o: tclXkeylist.h
tcl_utils.o: tcl_utils.h
tcl_utils.o: text_output.h
text_output.o: $(SRCROOT)/Misc/FtoC.h
text_output.o: $(SRCROOT)/Misc/getfile.h
text_output.o: $(SRCROOT)/Misc/misc.h
text_output.o: $(SRCROOT)/Misc/os.h
text_output.o: $(SRCROOT)/Misc/vlen.h
text_output.o: $(SRCROOT)/Misc/xalloc.h
text_output.o: $(PWD)/staden_config.h
text_output.o: tcl_utils.h
text_output.o: text_output.h
tkCanvGraph.o: $(SRCROOT)/Misc/os.h
tkCanvGraph.o: $(PWD)/staden_config.h
tkCanvGraph.o: matrix.h
tkCanvGraph.o: tclCanvGraph.h
tkCanvGraph.o: tkCanvGraph.h
tkRaster.o: $(SRCROOT)/Misc/xalloc.h
tkRaster.o: tcl_utils.h
tkRaster.o: tkRaster.h
tkRaster.o: tkRasterBuiltIn.h
tkRasterBuiltIn.o: tkRaster.h
tkRasterBuiltIn.o: tkRasterBuiltIn.h
tkSheet.o: $(SRCROOT)/Misc/xalloc.h
tkSheet.o: tcl_utils.h
tkSheet.o: tkSheet.h sheet.h
tkSheet.o: tkSheet_common.h
tkSheet.o: tkSheet_config.h
tkSheet.o: tkSheet_struct.h
tkSheet.o: tk_defs.h
tkSheet_common.o: $(SRCROOT)/Misc/xalloc.h
tkSheet_common.o: sheet.h
tkSheet_common.o: tkSheet.h
tkSheet_common.o: tkSheet_common.h
tkSheet_common.o: tkSheet_struct.h
tkSheet_common.o: tk_defs.h
tkTrace.o: $(SRCROOT)/Misc/misc.h
tkTrace.o: $(SRCROOT)/Misc/os.h
tkTrace.o: $(SRCROOT)/Misc/xalloc.h
tkTrace.o: $(PWD)/staden_config.h
tkTrace.o: cli_arg.h
tkTrace.o: postscript.h
tkTrace.o: split.h
tkTrace.o: tcl_utils.h
tkTrace.o: tkTrace.h
tkTrace.o: tkTraceIO.h
tkTrace.o: tk_defs.h
tkTrace.o: trace_print.h
tkTraceComp.o: cli_arg.h
tkTraceComp.o: postscript.h
tkTraceComp.o: tkTrace.h
tkTraceComp.o: tkTraceIO.h
tkTraceDisp.o: $(SRCROOT)/Misc/misc.h
tkTraceDisp.o: $(SRCROOT)/Misc/os.h
tkTraceDisp.o: $(SRCROOT)/Misc/xalloc.h
tkTraceDisp.o: $(PWD)/staden_config.h
tkTraceDisp.o: cli_arg.h
tkTraceDisp.o: postscript.h
tkTraceDisp.o: tkTrace.h
tkTraceDisp.o: tkTraceIO.h
tkTraceIO.o: $(SRCROOT)/Misc/misc.h
tkTraceIO.o: $(SRCROOT)/Misc/os.h
tkTraceIO.o: $(SRCROOT)/Misc/xalloc.h
tkTraceIO.o: $(PWD)/staden_config.h
tkTraceIO.o: cli_arg.h
tkTraceIO.o: postscript.h
tkTraceIO.o: tkTrace.h
tkTraceIO.o: trace_print.h
trace_print.o: $(SRCROOT)/Misc/misc.h
trace_print.o: $(SRCROOT)/Misc/os.h
trace_print.o: $(SRCROOT)/Misc/xalloc.h
trace_print.o: $(PWD)/staden_config.h
trace_print.o: cli_arg.h
trace_print.o: postscript.h
trace_print.o: split.h
trace_print.o: tkTrace.h
trace_print.o: trace_print.h
