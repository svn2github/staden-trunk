noinst_LTLIBRARIES = libztr.la

sources = \
	compression.c \
	compression.h \
	ztr.c \
	ztr.h \
	ztr_translate.c

if HAVE_ZLIB
libztr_la_SOURCES = ${sources}
else
libztr_la_SOURCES =
endif

libztr_la_LIBADD = @LIBZ@ @LIBM@

INCLUDES = \
	-I${top_srcdir} \
	-I${top_srcdir}/read \
	-I${top_srcdir}/alf \
	-I${top_srcdir}/abi \
	-I${top_srcdir}/ctf \
	-I${top_srcdir}/ztr \
	-I${top_srcdir}/plain \
	-I${top_srcdir}/scf \
	-I${top_srcdir}/exp_file \
	-I${top_srcdir}/utils

EXTRA_DIST = FORMAT

pkginclude_HEADERS = \
	compression.h \
	ztr.h