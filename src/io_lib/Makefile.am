AUTOMAKE_OPTIONS = foreign no-dependencies

SUBDIRS = read progs

man_MANS = \
	man/man3/exp2read.3 \
	man/man3/ExperimentFile.3 \
	man/man3/fread_reading.3 \
	man/man3/fread_scf.3 \
	man/man3/fwrite_reading.3 \
	man/man3/fwrite_scf.3 \
	man/man3/read2exp.3 \
	man/man3/read2scf.3 \
	man/man4/Read.4 \
	man/man3/read_allocate.3 \
	man/man3/read_deallocate.3 \
	man/man3/read_reading.3 \
	man/man3/read_scf.3 \
	man/man3/read_scf_header.3 \
	man/man3/scf2read.3 \
	man/man3/write_reading.3 \
	man/man3/write_scf.3 \
	man/man3/write_scf_header.3

pkginclude_HEADERS = \
	read/Read.h \
	read/scf_extras.h \
	read/translate.h \
	abi/abi.h \
	abi/fpoint.h \
	abi/seqIOABI.h \
	alf/alf.h \
	ctf/seqIOCTF.h \
	exp_file/expFileIO.h \
	plain/plain.h \
	scf/scf.h \
	sff/sff.h \
	utils/array.h \
	utils/compress.h \
	utils/error.h \
	utils/mach-io.h \
	utils/misc.h \
	utils/open_trace_file.h \
	utils/tar_format.h \
	utils/traceType.h \
	utils/xalloc.h \
	utils/mFILE.h \
	utils/stdio_hack.h \
	utils/vlen.h \
	utils/hash_table.h \
	utils/jenkins_lookup3.h \
	utils/os.h \
	ztr/compression.h \
	ztr/ztr.h \
	ztr/deflate_simple.h \
	io_lib_config.h

bin_SCRIPTS = io_lib-config

EXTRA_DIST = README COPYRIGHT CHANGES man options.mk bootstrap \
	ztr/FORMAT io_lib-config.in

