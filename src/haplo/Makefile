LIBS = haplo
PROGS = x$(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E := $(MISC_INC) $(GAP4_INC) $(G_INC) $(TK_INC) $(ITK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(INCLUDES_E)
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

OBJS=\
	haplo_tcl.o \
	haplo_snps.o \
	haplo_split.o \
	haplo_cons.o

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

install:
	-mkdir $(INSTALLTCL)/haplo
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/haplo
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/haplorc $(INSTALLETC)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/Makefile $(S)/tclIndex $(S)/haplorc \
	   $(DIRNAME)

DEPEND_OBJ = $(OBJS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
