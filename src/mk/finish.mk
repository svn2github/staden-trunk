LIBS = finish
PROGS = $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..

include	 $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

# Local comment: Comment out next line for remote compilation
FINBIN=$(O)

#L=$(O)

INCLUDES_E := $(GAP4_INC) $(G_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(INCLUDES_E)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	$(FINBIN)/finish.o \
	$(FINBIN)/finish_hash.o \
	$(FINBIN)/finish_long.o \
	$(FINBIN)/finish_main.o \
	$(FINBIN)/finish_utils.o \
	$(FINBIN)/finish_walk.o

FIN_LIBS=\
	$(TK_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(MISC_LIB) \
	$(GAP_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	-$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(FIN_LIBS) $(SHLIB_DEP)

install:
	-mkdir $(INSTALLLIB)/finish
	cp finish.tcl sanger_names.tcl $(INSTALLLIB)/finish
	cp METHODS CHEMISTRY $(INSTALLLIB)/finish
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

distsrc: distsrc_dirs
	-cp -R *.[ch] finish.tcl sanger_names.tcl \
		Makefile dependencies METHODS CHEMISTRY \
		$(DIRNAME)

DEPEND_OBJ = $(OBJS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
