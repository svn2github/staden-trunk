noinst_LTLIBRARIES = libutils.la

libutils_la_SOURCES = \
	array.c \
	array.h \
	compress.c \
	compress.h \
	error.c \
	error.h \
	files.c \
	find.c \
	mach-io.c \
	mach-io.h \
	misc.h \
	open_trace_file.c \
	open_trace_file.h \
	read_alloc.c \
	strings.c \
	tar_format.h \
	traceType.c \
	traceType.h \
	xalloc.c \
	xalloc.h

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

EXTRA_DIST = README

pkginclude_HEADERS = \
	array.h \
	compress.h \
	error.h \
	mach-io.h \
	misc.h \
	open_trace_file.h \
	tar_format.h \
	traceType.h \
	xalloc.h
