lib_LTLIBRARIES = libread.la

libread_la_SOURCES = \
	Read.c \
	Read.h \
	scf_extras.c \
	scf_extras.h \
	translate.c \
	translate.h

libread_la_LIBADD = \
	${top_builddir}/utils/libutils.la \
	${top_builddir}/abi/libabi.la \
	${top_builddir}/alf/libalf.la \
	${top_builddir}/ctf/libctf.la \
	${top_builddir}/exp_file/libexp_file.la \
	${top_builddir}/plain/libplain.la \
	${top_builddir}/scf/libscf.la \
	${top_builddir}/ztr/libztr.la

libread_la_LDFLAGS = -avoid-version

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
	-I${top_srcdir}/exp_file \
	-I${top_srcdir}/utils

pkginclude_HEADERS = \
	Read.h \
	scf_extras.h \
	translate.h

EXTRA_DIST = README
