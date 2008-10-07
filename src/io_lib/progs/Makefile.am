bin_PROGRAMS = convert_trace makeSCF extract_seq extract_qual extract_fastq index_tar scf_dump scf_info scf_update get_comment hash_tar hash_extract hash_list trace_dump hash_sff append_sff ztr_dump illumina2srf srf2illumina srf_dump_all srf_index_hash srf_extract_linear srf_extract_hash srf2fastq srf2fasta srf_filter srf_info srf_list

convert_trace_SOURCES = convert_trace.c
convert_trace_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

get_comment_SOURCES = get_comment.c
get_comment_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

makeSCF_SOURCES = makeSCF.c
makeSCF_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

extract_seq_SOURCES = extract_seq.c
extract_seq_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

extract_qual_SOURCES = extract_qual.c
extract_qual_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

extract_fastq_SOURCES = extract_fastq.c
extract_fastq_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

index_tar_SOURCES = index_tar.c

scf_dump_SOURCES = scf_dump.c
scf_dump_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

scf_info_SOURCES = scf_info.c
scf_info_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

scf_update_SOURCES = scf_update.c
scf_update_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

hash_tar_SOURCES = hash_tar.c
hash_tar_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

hash_list_SOURCES = hash_list.c
hash_list_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

hash_sff_SOURCES = hash_sff.c
hash_sff_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

append_sff_SOURCES = append_sff.c
append_sff_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

hash_extract_SOURCES = hash_extract.c
hash_extract_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

trace_dump_SOURCES = trace_dump.c
trace_dump_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

ztr_dump_SOURCES = ztr_dump.c
ztr_dump_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

illumina2srf_SOURCES = solexa2srf.c
illumina2srf_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf2illumina_SOURCES = srf2solexa.c
srf2illumina_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf2fastq_SOURCES = srf2fastq.c
srf2fastq_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf2fasta_SOURCES = srf2fasta.c
srf2fasta_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf_dump_all_SOURCES = srf_dump_all.c
srf_dump_all_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf_extract_linear_SOURCES = srf_extract_linear.c
srf_extract_linear_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf_extract_hash_SOURCES = srf_extract_hash.c
srf_extract_hash_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf_index_hash_SOURCES = srf_index_hash.c
srf_index_hash_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf_filter_SOURCES = srf_filter.c
srf_filter_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf_info_SOURCES = srf_info.c
srf_info_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

srf_list_SOURCES = srf_list.c
srf_list_LDADD = $(top_builddir)/io_lib/libread.a @LIBZ@ @LIBCURL@

INCLUDES= -I${top_srcdir}
