LIBS = haplo
PROGS = x$(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

# Local comment: Comment out next line for remote compilation
HAPLOBIN=$(O)

#L=$(O)

INCLUDES_E := $(MISC_INC) $(GAP4_INC) $(G_INC) $(TK_INC) $(ITK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(INCLUDES_E)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	$(HAPLOBIN)/haplo_tcl.o \
	$(HAPLOBIN)/haplo_snps.o \
	$(HAPLOBIN)/haplo_split.o \
	$(HAPLOBIN)/haplo_cons.o

HAPLO_DEP=\
	$(TK_LIB) \
	$(TKUTILS_LIB) \
	$(MISC_LIB) \
	$(GAP_LIB)

x$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(HAPLO_DEP) $(SHLIB_DEP)

include $(SRCROOT)/gap4/gap4_defs.mk

install:
	-mkdir $(INSTALLLIB)/haplo
	cp *.tcl tclIndex $(INSTALLLIB)/haplo
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl Makefile dependencies tclIndex \
	    $(DIRNAME)

DEPEND_OBJ = $(OBJS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
