LIBS = tgap
PROGS = $(LIBS) $(PROGBINS)
PROGBINS= $(O)/tg_index $(O)/tg_view
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
#SRCROOT=/home5/pubseq/share/src
#REMOTESRC=$(GAP4SRC)

include	$(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

# Local comment: Comment out next line for remote compilation
TGAPBIN=$(O)

INCLUDES_E += $(MISC_INC) $(TCL_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(IOLIB_INC)\
	-I$(SRCROOT)/samtools-0.1.2-$(MACHINE) -I$(SRCROOT)/tgap
CFLAGS += $(SHLIB_CFLAGS)

# For gprof profiling
#CFLAGS   += -pg
#CLDFLAGS += -pg

# For debugging g-alloc
#CFLAGS += -DVALGRIND
#CFLAGS += -D_FORTIFY_SOURCE=2 -fstack-protector

TGAP_LIB = -L$(L) -ltgap

#-----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# "g" library components (from Gap4's src/g)

GSIO = \
	$(TGAPBIN)/g-files.o \
	$(TGAPBIN)/g-db.o \
	$(TGAPBIN)/g-struct.o

GSMISC = \
	$(TGAPBIN)/g-error.o \
	$(TGAPBIN)/g-alloc.o \
	$(TGAPBIN)/g-io.o

GSREQ = \
	$(TGAPBIN)/g-connect.o \
	$(TGAPBIN)/g-request.o

GOBJS = $(GSIO) $(GSMISC) $(GSREQ)

TG_IO = \
	$(TGAPBIN)/tg_iface_g.o \
	$(TGAPBIN)/b+tree2.o \
	$(TGAPBIN)/hache_table.o \
	$(TGAPBIN)/tg_gio.o \
	$(TGAPBIN)/tg_cache.o \
	$(TGAPBIN)/tg_contig.o \
	$(TGAPBIN)/tg_sequence.o \
	$(TGAPBIN)/tg_track.o \
	$(TGAPBIN)/tg_bin.o \
	$(TGAPBIN)/tg_utils.o \
	$(TGAPBIN)/tg_register.o \
	$(TGAPBIN)/tg_tracks.o \
	$(TGAPBIN)/tg_anno.o \
	$(TGAPBIN)/tg_library.o

TGSH_OBJ = \
	$(TG_IO) \
	$(TGAPBIN)/gap_cli_arg.o \
	$(TGAPBIN)/tg_tcl.o

TGAPOBJS = $(TGSH_OBJ) $(GOBJS) $(L)/libptmalloc3.a


$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(TGAPOBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LINK_PATHFLAG)$(L) $(TGAPOBJS) $(TGAP_LIBS) $(SHLIB_DEP)

TG_IND_OBJ = \
	$(TG_IO) \
	$(TGAPBIN)/tg_index.o \
	$(TGAPBIN)/baf.o \
	$(TGAPBIN)/maq.o \
	$(TGAPBIN)/maqmap.o \
	$(TGAPBIN)/sam.o \
	$(TGAPBIN)/ace.o
#	mapass2.o

TGILIBS = \
	$(TGAP_LIB) \
	$(MISC_LIB) \
	$(ZLIB_LIB) \
	$(TKUTILS_LIB) \
	$(SEQUTILS_LIB) \
	$(TCL_LIB) \
	-L../samtools-0.1.2-$(MACHINE) -lbam

$(O)/tg_index: $(TG_IND_OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_IND_OBJ) $(TGILIBS) $(LIBSC)

TG_VIEW_OBJ = \
	$(TG_IO) \
	$(TGAPBIN)/tg_view.o

TGVLIBS = \
	$(TGILIBS) \
	-lcurses

$(O)/tg_view: $(TG_VIEW_OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_VIEW_OBJ) $(TGVLIBS) $(LIBSC)

TG_DEPTH_OBJ = \
	$(TG_IO) \
	$(TGAPBIN)/tg_depth.o

TGDLIBS = \
	$(TGILIBS) \

$(O)/tg_depth: $(TG_DEPTH_OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_DEPTH_OBJ) $(TGDLIBS) $(LIBSC)

DEPEND_OBJ = $(TG_IND_OBJ) $(TG_VIEW_OBJ) $(TGAPOBJS)

install:
	cp $(PROGBINS) $(INSTALLBIN)
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)
