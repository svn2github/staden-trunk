lib_LIBRARIES = libread.a

libread_a_SOURCES = \
	Read.c \
	Read.h \
	scf_extras.c \
	scf_extras.h \
	translate.c \
	translate.h \
	../ztr/compression.c \
	../ztr/compression.h \
	../ztr/ztr.c \
	../ztr/ztr.h \
	../ztr/ztr_translate.c \
	../ztr/deflate_simple.c \
	../ztr/deflate_simple.h \
	../abi/abi.h \
	../abi/fpoint.c \
	../abi/fpoint.h \
	../abi/seqIOABI.c \
	../abi/seqIOABI.h \
	../alf/seqIOALF.c \
	../alf/alf.h \
	../ctf/ctfCompress.c \
	../ctf/seqIOCTF.c \
	../ctf/seqIOCTF.h \
	../exp_file/expFileIO.c \
	../exp_file/expFileIO.h \
	../plain/plain.h \
	../plain/seqIOPlain.c \
	../scf/misc_scf.c \
	../scf/read_scf.c \
	../scf/scf.h \
	../scf/write_scf.c \
	../sff/sff.c \
	../sff/sff.h \
	../utils/array.c \
	../utils/array.h \
	../utils/compress.c \
	../utils/compress.h \
	../utils/error.c \
	../utils/error.h \
	../utils/files.c \
	../utils/find.c \
	../utils/mach-io.c \
	../utils/mach-io.h \
	../utils/misc.h \
	../utils/open_trace_file.c \
	../utils/open_trace_file.h \
	../utils/read_alloc.c \
	../utils/strings.c \
	../utils/tar_format.h \
	../utils/traceType.c \
	../utils/traceType.h \
	../utils/xalloc.c \
	../utils/xalloc.h \
	../utils/vlen.c \
	../utils/vlen.h \
	../utils/hash_table.c \
	../utils/hash_table.h \
	../utils/jenkins_lookup3.c \
	../utils/jenkins_lookup3.h \
	../utils/mFILE.c \
	../utils/mFILE.h \
	../utils/os.h \
	../utils/stdio_hack.h

libread_a_CPPFLAGS = @LIBCURL_CPPFLAGS@

#libread_a_LDFLAGS = -avoid-version

INCLUDES = \
	-I${top_srcdir} \
	-I${top_srcdir}/include \
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
