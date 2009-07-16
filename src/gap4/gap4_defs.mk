#
# Here we define the necessary object files needed for an external program
# to link with the gap4 I/O routines
#
#
# The object file list are in  $(GAPDB_EXT_OBJS).
# The linker libraries are in  $(GAPDB_EXT_LIBS).
# The link dependencies are in $(GAPDB_EXT_DEPS).
# To use these, simply add "include $(SRCROOT)/mk/gap4_defs.mk" to your
# Makefile after the global and machine includes.
#
# For an example of usage, see the "convert" program Makefile.

GAPDB_LOW=\
	$(GAP4BIN)/gap-local.o\
	$(GAP4BIN)/gap-remote.o\
	$(GAP4BIN)/gap-if.o\
	$(GAP4BIN)/gap-init.o\
	$(GAP4BIN)/gap-dbstruct.o\
	$(GAP4BIN)/gap-create.o\
	$(GAP4BIN)/gap-error.o\
	$(GAP4BIN)/stack_dump.o\
	$(GAP4BIN)/gap-io.o

GAPDB_MID=\
        $(GAP4BIN)/IO.o \
        $(GAP4BIN)/io-reg.o \
        $(GAP4BIN)/actf.o

GAPDB_UTILS=\
	$(GAP4BIN)/io_handle.o \
	$(GAP4BIN)/io_utils.o

# GAPDB_EXT_OBJS is basically the low and mid level files plus the
# text-io-reg.o object. This is simply stub routines to enable safe linking
# and isn't needed by gap4 itself.

GAPDB_EXT_OBJS=\
	$(GAPDB_LOW) \
	$(GAPDB_MID) \
	$(GAPDB_UTILS) \
	$(GAP4BIN)/text-io-reg.o

GAPDB_EXT_INC=$(G_INC) -I$(GAP4SRC)

GAPDB_EXT_DEPS=\
	$(G_DEP) \
	$(TEXTUTILS_DEP) \
	$(MISC_DEP) \
	$(TCL_DEP)

GAPDB_EXT_LIBS=\
	$(G_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB) \
	$(TCL_LIB)
