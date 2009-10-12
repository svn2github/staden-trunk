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
	tg_library.o \
	tg_index_common.o\
	zfio.o

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
	sam_index.o \
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

tg_index: $(TG_IND_OBJ) $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_IND_OBJ) $(TGILIBS) $(LIBSC)

TG_VIEW_OBJ = \
	$(TG_IO) \
	tg_view.o

TGVLIBS = \
	$(TGILIBS) \
	-lcurses

tg_view: $(TG_VIEW_OBJ) $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_VIEW_OBJ) $(TGVLIBS) $(LIBSC)

DEPEND_OBJ = $(TG_IND_OBJ) $(TG_VIEW_OBJ) $(TGAPOBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.pl $(S)/BAF $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGBINS) $(INSTALLBIN)
	cp $(PROGLIBS) $(INSTALLLIB)

# DO NOT DELETE THIS LINE -- make depend depends on it.

ace.o: $(SRCROOT)/Misc/array.h
ace.o: $(SRCROOT)/Misc/os.h
ace.o: $(SRCROOT)/Misc/tree.h
ace.o: $(SRCROOT)/seq_utils/dna_utils.h
ace.o: $(SRCROOT)/tgap/g-connect.h
ace.o: $(SRCROOT)/tgap/g-db.h
ace.o: $(SRCROOT)/tgap/g-defs.h
ace.o: $(SRCROOT)/tgap/g-error.h
ace.o: $(SRCROOT)/tgap/g-io.h
ace.o: $(SRCROOT)/tgap/g-misc.h
ace.o: $(SRCROOT)/tgap/g-request.h
ace.o: $(SRCROOT)/tgap/g-struct.h
ace.o: $(SRCROOT)/tgap/tg_index.h
ace.o: $(PWD)/staden_config.h
ace.o: ace.h
ace.o: g-alloc.h
ace.o: g-filedefs.h g-os.h
ace.o: tg_gio.h tg_struct.h g.h
ace.o: tg_iface.h tg_cache_item.h
ace.o: tg_library.h misc.h
ace.o: tg_register.h tg_contig.h
ace.o: tg_sequence.h tg_bin.h
ace.o: tg_tcl.h
ace.o: tg_track.h tg_anno.h
ace.o: xalloc.h tg_utils.h
ace.o: xerror.h hache_table.h
b+tree2.o: b+tree2.h
baf.o: $(SRCROOT)/Misc/array.h
baf.o: $(SRCROOT)/Misc/os.h
baf.o: $(SRCROOT)/Misc/tree.h
baf.o: $(SRCROOT)/tgap/g-connect.h
baf.o: $(SRCROOT)/tgap/g-db.h
baf.o: $(SRCROOT)/tgap/g-defs.h
baf.o: $(SRCROOT)/tgap/g-error.h
baf.o: $(SRCROOT)/tgap/g-io.h
baf.o: $(SRCROOT)/tgap/g-misc.h
baf.o: $(SRCROOT)/tgap/g-request.h
baf.o: $(SRCROOT)/tgap/g-struct.h
baf.o: $(SRCROOT)/tgap/tg_index.h
baf.o: $(PWD)/staden_config.h
baf.o: baf.h
baf.o: g-alloc.h
baf.o: g-filedefs.h g-os.h
baf.o: tg_gio.h tg_struct.h g.h
baf.o: tg_iface.h tg_cache_item.h
baf.o: tg_library.h misc.h
baf.o: tg_register.h tg_contig.h
baf.o: tg_sequence.h tg_bin.h
baf.o: tg_tcl.h
baf.o: tg_track.h tg_anno.h
baf.o: xalloc.h tg_utils.h
baf.o: xerror.h hache_table.h
g-alloc.o: $(SRCROOT)/Misc/os.h
g-alloc.o: $(PWD)/staden_config.h
g-alloc.o: g-alloc.h
g-connect.o: $(SRCROOT)/Misc/array.h
g-connect.o: $(SRCROOT)/Misc/os.h
g-connect.o: $(SRCROOT)/tgap/g-connect.h
g-connect.o: $(SRCROOT)/tgap/g-db.h
g-connect.o: $(SRCROOT)/tgap/g-defs.h
g-connect.o: $(SRCROOT)/tgap/g-error.h
g-connect.o: $(SRCROOT)/tgap/g-struct.h
g-connect.o: $(PWD)/staden_config.h
g-connect.o: g-alloc.h
g-connect.o: g-filedefs.h g-os.h
g-connect.o: hache_table.h
g-connect.o: xerror.h
g-db.o: $(SRCROOT)/Misc/array.h
g-db.o: $(SRCROOT)/Misc/os.h
g-db.o: $(SRCROOT)/tgap/g-db.h
g-db.o: $(SRCROOT)/tgap/g-defs.h
g-db.o: $(SRCROOT)/tgap/g-error.h
g-db.o: $(SRCROOT)/tgap/g-io.h
g-db.o: $(SRCROOT)/tgap/g-request.h
g-db.o: $(SRCROOT)/tgap/g-struct.h
g-db.o: $(PWD)/staden_config.h
g-db.o: g-alloc.h
g-db.o: g-filedefs.h g-os.h
g-db.o: g-files.h
g-db.o: hache_table.h
g-db.o: xerror.h
g-error.o: $(SRCROOT)/tgap/g-error.h
g-error.o: $(SRCROOT)/tgap/g-misc.h
g-error.o: xerror.h
g-files.o: $(SRCROOT)/Misc/array.h
g-files.o: $(SRCROOT)/Misc/os.h
g-files.o: $(SRCROOT)/tgap/g-db.h
g-files.o: $(SRCROOT)/tgap/g-defs.h
g-files.o: $(SRCROOT)/tgap/g-error.h
g-files.o: $(SRCROOT)/tgap/g-io.h
g-files.o: $(SRCROOT)/tgap/g-struct.h
g-files.o: $(PWD)/staden_config.h
g-files.o: g-filedefs.h g-os.h
g-files.o: g-files.h
g-files.o: hache_table.h
g-files.o: xalloc.h
g-files.o: xerror.h g-alloc.h
g-io.o: $(SRCROOT)/Misc/array.h
g-io.o: $(SRCROOT)/Misc/os.h
g-io.o: $(SRCROOT)/tgap/g-io.h
g-io.o: $(SRCROOT)/tgap/g-struct.h
g-io.o: $(PWD)/staden_config.h
g-io.o: g-alloc.h
g-io.o: g-filedefs.h g-os.h
g-io.o: g-files.h
g-io.o: xerror.h hache_table.h
g-request.o: $(SRCROOT)/Misc/array.h
g-request.o: $(SRCROOT)/Misc/os.h
g-request.o: $(SRCROOT)/tgap/g-db.h
g-request.o: $(SRCROOT)/tgap/g-defs.h
g-request.o: $(SRCROOT)/tgap/g-error.h
g-request.o: $(SRCROOT)/tgap/g-misc.h
g-request.o: $(SRCROOT)/tgap/g-request.h
g-request.o: $(SRCROOT)/tgap/g-struct.h
g-request.o: $(PWD)/staden_config.h
g-request.o: g-filedefs.h g-os.h
g-request.o: g-files.h
g-request.o: hache_table.h
g-request.o: xalloc.h
g-request.o: xerror.h g-alloc.h
g-struct.o: $(SRCROOT)/Misc/array.h
g-struct.o: $(SRCROOT)/Misc/os.h
g-struct.o: $(SRCROOT)/tgap/g-io.h
g-struct.o: $(SRCROOT)/tgap/g-struct.h
g-struct.o: $(PWD)/staden_config.h
g-struct.o: g-filedefs.h g-os.h
g-struct.o: g-files.h
g-struct.o: hache_table.h
g-struct.o: xalloc.h
g-struct.o: xerror.h g-alloc.h
gap_cli_arg.o: $(SRCROOT)/Misc/array.h
gap_cli_arg.o: $(SRCROOT)/Misc/os.h
gap_cli_arg.o: $(SRCROOT)/Misc/tree.h
gap_cli_arg.o: $(SRCROOT)/tgap/g-connect.h
gap_cli_arg.o: $(SRCROOT)/tgap/g-db.h
gap_cli_arg.o: $(SRCROOT)/tgap/g-defs.h
gap_cli_arg.o: $(SRCROOT)/tgap/g-error.h
gap_cli_arg.o: $(SRCROOT)/tgap/g-io.h
gap_cli_arg.o: $(SRCROOT)/tgap/g-misc.h
gap_cli_arg.o: $(SRCROOT)/tgap/g-request.h
gap_cli_arg.o: $(SRCROOT)/tgap/g-struct.h
gap_cli_arg.o: $(PWD)/staden_config.h
gap_cli_arg.o: g-alloc.h
gap_cli_arg.o: g-filedefs.h
gap_cli_arg.o: g-os.h
gap_cli_arg.o: gap_cli_arg.h
gap_cli_arg.o: hache_table.h
gap_cli_arg.o: misc.h
gap_cli_arg.o: tg_anno.h
gap_cli_arg.o: tg_bin.h
gap_cli_arg.o: tg_cache_item.h
gap_cli_arg.o: tg_contig.h
gap_cli_arg.o: tg_iface.h
gap_cli_arg.o: tg_library.h
gap_cli_arg.o: tg_register.h
gap_cli_arg.o: tg_sequence.h
gap_cli_arg.o: tg_struct.h g.h
gap_cli_arg.o: tg_tcl.h
gap_cli_arg.o: tg_track.h
gap_cli_arg.o: tg_utils.h
gap_cli_arg.o: xalloc.h tg_gio.h
gap_cli_arg.o: xerror.h
hache_table.o: $(SRCROOT)/Misc/os.h
hache_table.o: $(PWD)/staden_config.h
hache_table.o: hache_table.h
maq.o: $(SRCROOT)/Misc/array.h
maq.o: $(SRCROOT)/Misc/os.h
maq.o: $(SRCROOT)/Misc/tree.h
maq.o: $(SRCROOT)/tgap/g-connect.h
maq.o: $(SRCROOT)/tgap/g-db.h
maq.o: $(SRCROOT)/tgap/g-defs.h
maq.o: $(SRCROOT)/tgap/g-error.h
maq.o: $(SRCROOT)/tgap/g-io.h
maq.o: $(SRCROOT)/tgap/g-misc.h
maq.o: $(SRCROOT)/tgap/g-request.h
maq.o: $(SRCROOT)/tgap/g-struct.h
maq.o: $(SRCROOT)/tgap/tg_index.h
maq.o: $(PWD)/staden_config.h
maq.o: g-alloc.h
maq.o: g-filedefs.h g-os.h
maq.o: hache_table.h
maq.o: maqmap.h maq.h
maq.o: tg_iface.h tg_cache_item.h
maq.o: tg_library.h misc.h
maq.o: tg_register.h tg_contig.h
maq.o: tg_sequence.h tg_bin.h
maq.o: tg_struct.h g.h
maq.o: tg_tcl.h
maq.o: tg_track.h tg_anno.h
maq.o: xalloc.h tg_utils.h
maq.o: xerror.h tg_gio.h
maqmap.o: maqmap.h
sam_index.o: $(SRCROOT)/Misc/array.h
sam_index.o: $(SRCROOT)/Misc/os.h
sam_index.o: $(SRCROOT)/Misc/tree.h
sam_index.o: $(SRCROOT)/tgap/g-connect.h
sam_index.o: $(SRCROOT)/tgap/g-db.h
sam_index.o: $(SRCROOT)/tgap/g-defs.h
sam_index.o: $(SRCROOT)/tgap/g-error.h
sam_index.o: $(SRCROOT)/tgap/g-io.h
sam_index.o: $(SRCROOT)/tgap/g-misc.h
sam_index.o: $(SRCROOT)/tgap/g-request.h
sam_index.o: $(SRCROOT)/tgap/g-struct.h
sam_index.o: $(SRCROOT)/tgap/tg_index.h
sam_index.o: $(PWD)/staden_config.h
sam_index.o: g-alloc.h
sam_index.o: g-filedefs.h g-os.h
sam_index.o: sam_index.h
sam_index.o: tg_gio.h tg_struct.h g.h
sam_index.o: tg_iface.h tg_cache_item.h
sam_index.o: tg_library.h misc.h
sam_index.o: tg_register.h tg_contig.h
sam_index.o: tg_sequence.h tg_bin.h
sam_index.o: tg_tcl.h
sam_index.o: tg_track.h tg_anno.h
sam_index.o: xalloc.h tg_utils.h
sam_index.o: xerror.h hache_table.h
tg_anno.o: $(SRCROOT)/Misc/array.h
tg_anno.o: $(SRCROOT)/Misc/os.h
tg_anno.o: $(SRCROOT)/Misc/tree.h
tg_anno.o: $(SRCROOT)/tgap/g-connect.h
tg_anno.o: $(SRCROOT)/tgap/g-db.h
tg_anno.o: $(SRCROOT)/tgap/g-defs.h
tg_anno.o: $(SRCROOT)/tgap/g-error.h
tg_anno.o: $(SRCROOT)/tgap/g-io.h
tg_anno.o: $(SRCROOT)/tgap/g-misc.h
tg_anno.o: $(SRCROOT)/tgap/g-request.h
tg_anno.o: $(SRCROOT)/tgap/g-struct.h
tg_anno.o: $(PWD)/staden_config.h
tg_anno.o: g-alloc.h
tg_anno.o: g-filedefs.h g-os.h
tg_anno.o: g.h
tg_anno.o: tg_cache_item.h
tg_anno.o: tg_contig.h
tg_anno.o: tg_gio.h tg_struct.h
tg_anno.o: tg_iface.h
tg_anno.o: tg_library.h misc.h
tg_anno.o: tg_register.h
tg_anno.o: tg_sequence.h tg_bin.h
tg_anno.o: tg_tcl.h
tg_anno.o: tg_track.h tg_anno.h
tg_anno.o: xalloc.h tg_utils.h
tg_anno.o: xerror.h hache_table.h
tg_bin.o: $(SRCROOT)/Misc/array.h
tg_bin.o: $(SRCROOT)/Misc/os.h
tg_bin.o: $(SRCROOT)/Misc/tree.h
tg_bin.o: $(SRCROOT)/tgap/g-connect.h
tg_bin.o: $(SRCROOT)/tgap/g-db.h
tg_bin.o: $(SRCROOT)/tgap/g-defs.h
tg_bin.o: $(SRCROOT)/tgap/g-error.h
tg_bin.o: $(SRCROOT)/tgap/g-io.h
tg_bin.o: $(SRCROOT)/tgap/g-misc.h
tg_bin.o: $(SRCROOT)/tgap/g-request.h
tg_bin.o: $(SRCROOT)/tgap/g-struct.h
tg_bin.o: $(PWD)/staden_config.h
tg_bin.o: g-alloc.h
tg_bin.o: g-filedefs.h g-os.h
tg_bin.o: g.h
tg_bin.o: tg_cache_item.h
tg_bin.o: tg_contig.h
tg_bin.o: tg_gio.h tg_struct.h
tg_bin.o: tg_iface.h
tg_bin.o: tg_library.h misc.h
tg_bin.o: tg_register.h
tg_bin.o: tg_sequence.h tg_bin.h
tg_bin.o: tg_track.h tg_anno.h
tg_bin.o: tg_tracks.h
tg_bin.o: tg_utils.h tg_tcl.h
tg_bin.o: xalloc.h
tg_bin.o: xerror.h hache_table.h
tg_cache.o: $(SRCROOT)/Misc/array.h
tg_cache.o: $(SRCROOT)/Misc/os.h
tg_cache.o: $(SRCROOT)/Misc/tree.h
tg_cache.o: $(SRCROOT)/tgap/g-connect.h
tg_cache.o: $(SRCROOT)/tgap/g-db.h
tg_cache.o: $(SRCROOT)/tgap/g-defs.h
tg_cache.o: $(SRCROOT)/tgap/g-error.h
tg_cache.o: $(SRCROOT)/tgap/g-io.h
tg_cache.o: $(SRCROOT)/tgap/g-misc.h
tg_cache.o: $(SRCROOT)/tgap/g-request.h
tg_cache.o: $(SRCROOT)/tgap/g-struct.h
tg_cache.o: $(PWD)/staden_config.h
tg_cache.o: g-alloc.h
tg_cache.o: g-filedefs.h g-os.h
tg_cache.o: g.h
tg_cache.o: hache_table.h
tg_cache.o: tg_anno.h
tg_cache.o: tg_bin.h tg_track.h
tg_cache.o: tg_cache_item.h
tg_cache.o: tg_contig.h
tg_cache.o: tg_gio.h tg_struct.h
tg_cache.o: tg_iface.h
tg_cache.o: tg_library.h misc.h
tg_cache.o: tg_register.h
tg_cache.o: tg_sequence.h
tg_cache.o: tg_tcl.h
tg_cache.o: xalloc.h tg_utils.h
tg_cache.o: xerror.h
tg_contig.o: $(SRCROOT)/Misc/array.h
tg_contig.o: $(SRCROOT)/Misc/os.h
tg_contig.o: $(SRCROOT)/Misc/tree.h
tg_contig.o: $(SRCROOT)/tgap/g-connect.h
tg_contig.o: $(SRCROOT)/tgap/g-db.h
tg_contig.o: $(SRCROOT)/tgap/g-defs.h
tg_contig.o: $(SRCROOT)/tgap/g-error.h
tg_contig.o: $(SRCROOT)/tgap/g-io.h
tg_contig.o: $(SRCROOT)/tgap/g-misc.h
tg_contig.o: $(SRCROOT)/tgap/g-request.h
tg_contig.o: $(SRCROOT)/tgap/g-struct.h
tg_contig.o: $(PWD)/staden_config.h
tg_contig.o: g-alloc.h
tg_contig.o: g-filedefs.h g-os.h
tg_contig.o: g.h
tg_contig.o: hache_table.h
tg_contig.o: tg_anno.h
tg_contig.o: tg_bin.h tg_track.h
tg_contig.o: tg_cache_item.h
tg_contig.o: tg_contig.h
tg_contig.o: tg_gio.h tg_struct.h
tg_contig.o: tg_iface.h
tg_contig.o: tg_library.h misc.h
tg_contig.o: tg_register.h
tg_contig.o: tg_sequence.h
tg_contig.o: tg_tcl.h
tg_contig.o: xalloc.h tg_utils.h
tg_contig.o: xerror.h
tg_gio.o: $(SRCROOT)/Misc/array.h
tg_gio.o: $(SRCROOT)/Misc/os.h
tg_gio.o: $(SRCROOT)/Misc/tree.h
tg_gio.o: $(SRCROOT)/tgap/g-connect.h
tg_gio.o: $(SRCROOT)/tgap/g-db.h
tg_gio.o: $(SRCROOT)/tgap/g-defs.h
tg_gio.o: $(SRCROOT)/tgap/g-error.h
tg_gio.o: $(SRCROOT)/tgap/g-io.h
tg_gio.o: $(SRCROOT)/tgap/g-misc.h
tg_gio.o: $(SRCROOT)/tgap/g-request.h
tg_gio.o: $(SRCROOT)/tgap/g-struct.h
tg_gio.o: $(PWD)/staden_config.h
tg_gio.o: g-alloc.h
tg_gio.o: g-filedefs.h g-os.h
tg_gio.o: g.h
tg_gio.o: tg_cache_item.h
tg_gio.o: tg_contig.h
tg_gio.o: tg_gio.h tg_struct.h
tg_gio.o: tg_iface.h
tg_gio.o: tg_iface_g.h
tg_gio.o: tg_library.h misc.h
tg_gio.o: tg_register.h
tg_gio.o: tg_sequence.h tg_bin.h
tg_gio.o: tg_tcl.h
tg_gio.o: tg_track.h tg_anno.h
tg_gio.o: xalloc.h tg_utils.h
tg_gio.o: xerror.h hache_table.h
tg_iface_g.o: $(SRCROOT)/Misc/array.h
tg_iface_g.o: $(SRCROOT)/Misc/os.h
tg_iface_g.o: $(SRCROOT)/Misc/tree.h
tg_iface_g.o: $(SRCROOT)/seq_utils/dna_utils.h
tg_iface_g.o: $(SRCROOT)/tgap/g-connect.h
tg_iface_g.o: $(SRCROOT)/tgap/g-db.h
tg_iface_g.o: $(SRCROOT)/tgap/g-defs.h
tg_iface_g.o: $(SRCROOT)/tgap/g-error.h
tg_iface_g.o: $(SRCROOT)/tgap/g-io.h
tg_iface_g.o: $(SRCROOT)/tgap/g-misc.h
tg_iface_g.o: $(SRCROOT)/tgap/g-request.h
tg_iface_g.o: $(SRCROOT)/tgap/g-struct.h
tg_iface_g.o: $(PWD)/staden_config.h
tg_iface_g.o: b+tree2.h
tg_iface_g.o: g-alloc.h
tg_iface_g.o: g-filedefs.h g-os.h
tg_iface_g.o: g.h
tg_iface_g.o: hache_table.h
tg_iface_g.o: misc.h xalloc.h
tg_iface_g.o: tg_anno.h
tg_iface_g.o: tg_bin.h tg_track.h
tg_iface_g.o: tg_cache_item.h
tg_iface_g.o: tg_contig.h
tg_iface_g.o: tg_gio.h
tg_iface_g.o: tg_iface.h
tg_iface_g.o: tg_iface_g.h
tg_iface_g.o: tg_library.h
tg_iface_g.o: tg_register.h
tg_iface_g.o: tg_sequence.h
tg_iface_g.o: tg_struct.h
tg_iface_g.o: tg_tcl.h
tg_iface_g.o: tg_utils.h
tg_iface_g.o: xerror.h
tg_index.o: $(SRCROOT)/Misc/array.h
tg_index.o: $(SRCROOT)/Misc/os.h
tg_index.o: $(SRCROOT)/Misc/tree.h
tg_index.o: $(SRCROOT)/tgap/g-connect.h
tg_index.o: $(SRCROOT)/tgap/g-db.h
tg_index.o: $(SRCROOT)/tgap/g-defs.h
tg_index.o: $(SRCROOT)/tgap/g-error.h
tg_index.o: $(SRCROOT)/tgap/g-io.h
tg_index.o: $(SRCROOT)/tgap/g-misc.h
tg_index.o: $(SRCROOT)/tgap/g-request.h
tg_index.o: $(SRCROOT)/tgap/g-struct.h
tg_index.o: $(SRCROOT)/tgap/tg_index.h
tg_index.o: $(PWD)/staden_config.h
tg_index.o: ace.h baf.h sam_index.h
tg_index.o: g-alloc.h
tg_index.o: g-filedefs.h g-os.h
tg_index.o: hache_table.h
tg_index.o: maq.h
tg_index.o: maqmap.h
tg_index.o: tg_anno.h
tg_index.o: tg_bin.h tg_track.h
tg_index.o: tg_cache_item.h
tg_index.o: tg_contig.h
tg_index.o: tg_iface.h
tg_index.o: tg_library.h misc.h
tg_index.o: tg_register.h
tg_index.o: tg_sequence.h
tg_index.o: tg_struct.h g.h
tg_index.o: tg_tcl.h
tg_index.o: xalloc.h tg_utils.h
tg_index.o: xerror.h tg_gio.h
tg_library.o: $(SRCROOT)/Misc/array.h
tg_library.o: $(SRCROOT)/Misc/os.h
tg_library.o: $(SRCROOT)/Misc/tree.h
tg_library.o: $(SRCROOT)/tgap/g-connect.h
tg_library.o: $(SRCROOT)/tgap/g-db.h
tg_library.o: $(SRCROOT)/tgap/g-defs.h
tg_library.o: $(SRCROOT)/tgap/g-error.h
tg_library.o: $(SRCROOT)/tgap/g-io.h
tg_library.o: $(SRCROOT)/tgap/g-misc.h
tg_library.o: $(SRCROOT)/tgap/g-request.h
tg_library.o: $(SRCROOT)/tgap/g-struct.h
tg_library.o: $(PWD)/staden_config.h
tg_library.o: g-alloc.h
tg_library.o: g-filedefs.h g-os.h
tg_library.o: hache_table.h
tg_library.o: tg_anno.h
tg_library.o: tg_bin.h tg_track.h
tg_library.o: tg_cache_item.h
tg_library.o: tg_contig.h
tg_library.o: tg_gio.h
tg_library.o: tg_iface.h
tg_library.o: tg_library.h misc.h
tg_library.o: tg_register.h
tg_library.o: tg_sequence.h
tg_library.o: tg_struct.h g.h
tg_library.o: tg_utils.h tg_tcl.h
tg_library.o: xalloc.h
tg_library.o: xerror.h
tg_register.o: $(SRCROOT)/Misc/array.h
tg_register.o: $(SRCROOT)/Misc/os.h
tg_register.o: $(SRCROOT)/Misc/tree.h
tg_register.o: $(SRCROOT)/tgap/g-connect.h
tg_register.o: $(SRCROOT)/tgap/g-db.h
tg_register.o: $(SRCROOT)/tgap/g-defs.h
tg_register.o: $(SRCROOT)/tgap/g-error.h
tg_register.o: $(SRCROOT)/tgap/g-io.h
tg_register.o: $(SRCROOT)/tgap/g-misc.h
tg_register.o: $(SRCROOT)/tgap/g-request.h
tg_register.o: $(SRCROOT)/tgap/g-struct.h
tg_register.o: $(PWD)/staden_config.h
tg_register.o: g-alloc.h
tg_register.o: g-filedefs.h
tg_register.o: g-os.h
tg_register.o: hache_table.h
tg_register.o: misc.h xalloc.h
tg_register.o: tg_anno.h
tg_register.o: tg_bin.h
tg_register.o: tg_cache_item.h
tg_register.o: tg_contig.h
tg_register.o: tg_gio.h
tg_register.o: tg_iface.h
tg_register.o: tg_library.h
tg_register.o: tg_register.h
tg_register.o: tg_sequence.h
tg_register.o: tg_struct.h g.h
tg_register.o: tg_tcl.h
tg_register.o: tg_track.h
tg_register.o: tg_utils.h
tg_register.o: xerror.h
tg_sequence.o: $(SRCROOT)/Misc/array.h
tg_sequence.o: $(SRCROOT)/Misc/os.h
tg_sequence.o: $(SRCROOT)/Misc/tree.h
tg_sequence.o: $(SRCROOT)/tgap/g-connect.h
tg_sequence.o: $(SRCROOT)/tgap/g-db.h
tg_sequence.o: $(SRCROOT)/tgap/g-defs.h
tg_sequence.o: $(SRCROOT)/tgap/g-error.h
tg_sequence.o: $(SRCROOT)/tgap/g-io.h
tg_sequence.o: $(SRCROOT)/tgap/g-misc.h
tg_sequence.o: $(SRCROOT)/tgap/g-request.h
tg_sequence.o: $(SRCROOT)/tgap/g-struct.h
tg_sequence.o: $(PWD)/staden_config.h
tg_sequence.o: g-alloc.h
tg_sequence.o: g-filedefs.h
tg_sequence.o: g-os.h
tg_sequence.o: hache_table.h
tg_sequence.o: misc.h xalloc.h
tg_sequence.o: tg_anno.h
tg_sequence.o: tg_bin.h
tg_sequence.o: tg_cache_item.h
tg_sequence.o: tg_contig.h
tg_sequence.o: tg_gio.h
tg_sequence.o: tg_iface.h
tg_sequence.o: tg_library.h
tg_sequence.o: tg_register.h
tg_sequence.o: tg_sequence.h
tg_sequence.o: tg_struct.h g.h
tg_sequence.o: tg_tcl.h
tg_sequence.o: tg_track.h
tg_sequence.o: tg_utils.h
tg_sequence.o: xerror.h
tg_tcl.o: $(SRCROOT)/Misc/array.h
tg_tcl.o: $(SRCROOT)/Misc/os.h
tg_tcl.o: $(SRCROOT)/Misc/tree.h
tg_tcl.o: $(SRCROOT)/tgap/g-connect.h
tg_tcl.o: $(SRCROOT)/tgap/g-db.h
tg_tcl.o: $(SRCROOT)/tgap/g-defs.h
tg_tcl.o: $(SRCROOT)/tgap/g-error.h
tg_tcl.o: $(SRCROOT)/tgap/g-io.h
tg_tcl.o: $(SRCROOT)/tgap/g-misc.h
tg_tcl.o: $(SRCROOT)/tgap/g-request.h
tg_tcl.o: $(SRCROOT)/tgap/g-struct.h
tg_tcl.o: $(SRCROOT)/tk_utils/tcl_utils.h
tg_tcl.o: $(PWD)/staden_config.h
tg_tcl.o: g-alloc.h g-filedefs.h
tg_tcl.o: g-os.h hache_table.h
tg_tcl.o: tg_cache_item.h
tg_tcl.o: tg_contig.h
tg_tcl.o: tg_iface.h
tg_tcl.o: tg_library.h tg_utils.h
tg_tcl.o: tg_register.h
tg_tcl.o: tg_sequence.h tg_bin.h
tg_tcl.o: tg_struct.h g.h
tg_tcl.o: tg_tcl.h gap_cli_arg.h
tg_tcl.o: tg_track.h tg_anno.h
tg_tcl.o: xalloc.h tg_gio.h
tg_tcl.o: xerror.h misc.h
tg_track.o: $(SRCROOT)/Misc/array.h
tg_track.o: $(SRCROOT)/Misc/os.h
tg_track.o: $(SRCROOT)/Misc/tree.h
tg_track.o: $(SRCROOT)/tgap/g-connect.h
tg_track.o: $(SRCROOT)/tgap/g-db.h
tg_track.o: $(SRCROOT)/tgap/g-defs.h
tg_track.o: $(SRCROOT)/tgap/g-error.h
tg_track.o: $(SRCROOT)/tgap/g-io.h
tg_track.o: $(SRCROOT)/tgap/g-misc.h
tg_track.o: $(SRCROOT)/tgap/g-request.h
tg_track.o: $(SRCROOT)/tgap/g-struct.h
tg_track.o: $(PWD)/staden_config.h
tg_track.o: g-alloc.h
tg_track.o: g-filedefs.h g-os.h
tg_track.o: g.h
tg_track.o: hache_table.h
tg_track.o: tg_anno.h
tg_track.o: tg_bin.h tg_track.h
tg_track.o: tg_cache_item.h
tg_track.o: tg_contig.h
tg_track.o: tg_gio.h tg_struct.h
tg_track.o: tg_iface.h
tg_track.o: tg_library.h misc.h
tg_track.o: tg_register.h
tg_track.o: tg_sequence.h
tg_track.o: tg_tcl.h
tg_track.o: xalloc.h tg_utils.h
tg_track.o: xerror.h
tg_tracks.o: $(SRCROOT)/Misc/array.h
tg_tracks.o: $(SRCROOT)/Misc/os.h
tg_tracks.o: $(SRCROOT)/Misc/tree.h
tg_tracks.o: $(SRCROOT)/tgap/g-connect.h
tg_tracks.o: $(SRCROOT)/tgap/g-db.h
tg_tracks.o: $(SRCROOT)/tgap/g-defs.h
tg_tracks.o: $(SRCROOT)/tgap/g-error.h
tg_tracks.o: $(SRCROOT)/tgap/g-io.h
tg_tracks.o: $(SRCROOT)/tgap/g-misc.h
tg_tracks.o: $(SRCROOT)/tgap/g-request.h
tg_tracks.o: $(SRCROOT)/tgap/g-struct.h
tg_tracks.o: $(PWD)/staden_config.h
tg_tracks.o: g-alloc.h
tg_tracks.o: g-filedefs.h g-os.h
tg_tracks.o: g.h
tg_tracks.o: hache_table.h
tg_tracks.o: tg_anno.h
tg_tracks.o: tg_bin.h tg_track.h
tg_tracks.o: tg_cache_item.h
tg_tracks.o: tg_contig.h
tg_tracks.o: tg_gio.h tg_struct.h
tg_tracks.o: tg_iface.h
tg_tracks.o: tg_library.h misc.h
tg_tracks.o: tg_register.h
tg_tracks.o: tg_sequence.h
tg_tracks.o: tg_tcl.h
tg_tracks.o: tg_tracks.h
tg_tracks.o: xalloc.h tg_utils.h
tg_tracks.o: xerror.h
tg_utils.o: $(SRCROOT)/Misc/array.h
tg_utils.o: $(SRCROOT)/Misc/os.h
tg_utils.o: $(SRCROOT)/Misc/tree.h
tg_utils.o: $(SRCROOT)/tgap/g-connect.h
tg_utils.o: $(SRCROOT)/tgap/g-db.h
tg_utils.o: $(SRCROOT)/tgap/g-defs.h
tg_utils.o: $(SRCROOT)/tgap/g-error.h
tg_utils.o: $(SRCROOT)/tgap/g-io.h
tg_utils.o: $(SRCROOT)/tgap/g-misc.h
tg_utils.o: $(SRCROOT)/tgap/g-request.h
tg_utils.o: $(SRCROOT)/tgap/g-struct.h
tg_utils.o: $(PWD)/staden_config.h
tg_utils.o: g-alloc.h
tg_utils.o: g-filedefs.h g-os.h
tg_utils.o: g.h
tg_utils.o: hache_table.h
tg_utils.o: misc.h
tg_utils.o: tg_anno.h
tg_utils.o: tg_bin.h tg_track.h
tg_utils.o: tg_cache_item.h
tg_utils.o: tg_contig.h
tg_utils.o: tg_gio.h tg_struct.h
tg_utils.o: tg_iface.h
tg_utils.o: tg_library.h
tg_utils.o: tg_register.h
tg_utils.o: tg_sequence.h
tg_utils.o: tg_utils.h tg_tcl.h
tg_utils.o: xalloc.h
tg_utils.o: xerror.h
tg_view.o: $(SRCROOT)/Misc/array.h
tg_view.o: $(SRCROOT)/Misc/os.h
tg_view.o: $(SRCROOT)/Misc/tree.h
tg_view.o: $(SRCROOT)/tgap/g-connect.h
tg_view.o: $(SRCROOT)/tgap/g-db.h
tg_view.o: $(SRCROOT)/tgap/g-defs.h
tg_view.o: $(SRCROOT)/tgap/g-error.h
tg_view.o: $(SRCROOT)/tgap/g-io.h
tg_view.o: $(SRCROOT)/tgap/g-misc.h
tg_view.o: $(SRCROOT)/tgap/g-request.h
tg_view.o: $(SRCROOT)/tgap/g-struct.h
tg_view.o: $(PWD)/staden_config.h
tg_view.o: g-alloc.h g-filedefs.h
tg_view.o: g-os.h hache_table.h
tg_view.o: g.h
tg_view.o: tg_cache_item.h
tg_view.o: tg_contig.h
tg_view.o: tg_gio.h tg_struct.h
tg_view.o: tg_iface.h
tg_view.o: tg_library.h
tg_view.o: tg_register.h
tg_view.o: tg_sequence.h tg_bin.h
tg_view.o: tg_track.h tg_anno.h
tg_view.o: tg_utils.h tg_tcl.h
tg_view.o: xalloc.h
tg_view.o: xerror.h misc.h
