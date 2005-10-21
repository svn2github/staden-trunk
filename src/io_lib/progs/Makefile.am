bin_PROGRAMS = convert_trace makeSCF extract_seq index_tar scf_dump scf_info scf_update get_comment hash_tar hash_extract trace_dump hash_sff

convert_trace_SOURCES = convert_trace.c
convert_trace_LDADD = $(top_builddir)/read/libread.a @LIBZ@

get_comment_SOURCES = get_comment.c
get_comment_LDADD = $(top_builddir)/read/libread.a @LIBZ@

makeSCF_SOURCES = makeSCF.c
makeSCF_LDADD = $(top_builddir)/read/libread.a @LIBZ@

extract_seq_SOURCES = extract_seq.c
extract_seq_LDADD = $(top_builddir)/read/libread.a @LIBZ@

index_tar_SOURCES = index_tar.c

scf_dump_SOURCES = scf_dump.c
scf_dump_LDADD = $(top_builddir)/read/libread.a @LIBZ@

scf_info_SOURCES = scf_info.c
scf_info_LDADD = $(top_builddir)/read/libread.a @LIBZ@

scf_update_SOURCES = scf_update.c
scf_update_LDADD = $(top_builddir)/read/libread.a @LIBZ@

hash_tar_SOURCES = hash_tar.c
hash_tar_LDADD = $(top_builddir)/read/libread.a @LIBZ@

hash_sff_SOURCES = hash_sff.c
hash_sff_LDADD = $(top_builddir)/read/libread.a @LIBZ@

hash_extract_SOURCES = hash_extract.c
hash_extract_LDADD = $(top_builddir)/read/libread.a @LIBZ@

trace_dump_SOURCES = trace_dump.c
trace_dump_LDADD = $(top_builddir)/read/libread.a @LIBZ@

INCLUDES= \
	-I${top_srcdir} \
	-I${top_srcdir}/read \
	-I${top_srcdir}/alf \
	-I${top_srcdir}/abi \
	-I${top_srcdir}/ctf \
	-I${top_srcdir}/ztr \
	-I${top_srcdir}/plain \
	-I${top_srcdir}/scf \
	-I${top_srcdir}/sff \
	-I${top_srcdir}/exp_file \
	-I${top_srcdir}/utils
