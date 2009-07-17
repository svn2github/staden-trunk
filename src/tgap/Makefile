LIBS = tgap
PROGS = $(LIBS) $(PROGBINS)
PROGBINS= tg_index tg_view
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

# Local comment: Comment out next line for remote compilation

INCLUDES_E += $(MISC_INC) $(TCL_INC) $(TKUTILS_INC) $(SEQUTILS_INC) \
	$(IOLIB_INC) $(SAMTOOLS_INC)
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
	g-files.o \
	g-db.o \
	g-struct.o

GSMISC = \
	g-error.o \
	g-alloc.o \
	g-io.o

GSREQ = \
	g-connect.o \
	g-request.o

GOBJS = $(GSIO) $(GSMISC) $(GSREQ)

TG_IO = \
	tg_iface_g.o \
	b+tree2.o \
	hache_table.o \
	tg_gio.o \
	tg_cache.o \
	tg_contig.o \
	tg_sequence.o \
	tg_track.o \
	tg_bin.o \
	tg_utils.o \
	tg_register.o \
	tg_tracks.o \
	tg_anno.o \
	tg_library.o

TGSH_OBJ = \
	$(TG_IO) \
	gap_cli_arg.o \
	tg_tcl.o

#TGAPOBJS = $(TGSH_OBJ) $(GOBJS) $(L)/libptmalloc3.a
TGAPOBJS = $(TGSH_OBJ) $(GOBJS)


$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(TGAPOBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LINK_PATHFLAG)$(L) $(TGAPOBJS) $(TGAP_LIBS) $(SHLIB_DEP)

TG_IND_OBJ = \
	$(TG_IO) \
	tg_index.o \
	baf.o \
	maq.o \
	maqmap.o \
	sam.o \
	ace.o
#	mapass2.o

TGILIBS = \
	$(TGAP_LIB) \
	$(MISC_LIB) \
	$(ZLIB_LIB) \
	$(TKUTILS_LIB) \
	$(SEQUTILS_LIB) \
	$(TCL_LIB) \
	$(SAMTOOLS_LIB)

tg_index: $(TG_IND_OBJ) $(LIBS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_IND_OBJ) $(TGILIBS) $(LIBSC)

TG_VIEW_OBJ = \
	$(TG_IO) \
	tg_view.o

TGVLIBS = \
	$(TGILIBS) \
	-lcurses

tg_view: $(TG_VIEW_OBJ) $(LIBS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_VIEW_OBJ) $(TGVLIBS) $(LIBSC)

DEPEND_OBJ = $(TG_IND_OBJ) $(TG_VIEW_OBJ) $(TGAPOBJS)

install:
	cp $(PROGBINS) $(INSTALLBIN)
	cp $(PROGLIBS) $(INSTALLLIB)
