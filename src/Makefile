#
# Root level Makefile for Staden package.
# MACHINE define must be set before using this makefile. Set either using
# "gmake MACHINE=alpha" or as an environment variable.
#
# Valid values are "alpha", "al", "dec", "sgi", "solaris", and "sun".
#
# Invoking this makefile with JOB set will pass that as the item to make onto
# all sublevel makefiles. eg "gmake JOB=clean".
#

include	 mk/global.mk

DIRS =	Misc\
	io_lib\
	tk_utils\
	text_utils\
	seq_utils\
	mutlib\
	abi\
	alf\
	g\
	trev\
	primer3\
	primer3/src\
	gap4\
	prefinish\
	eba\
	init_exp\
	get_scf_field\
	spin\
	vector_clip\
	tracediff\
	screen_seq\
	qclip\
	polyA_clip\
	stops\
	pregap4\
	make_weights\
	find_renz\
	scripts\
	spin_emboss\
	spin2_emboss\
	copy_reads\
	hetins\
	mutscan

#	seqed\
#	spin2

DIRS_OLD = \
	convert\
	ted\
	gap/osp-bits\
	gap\
	frog\
	getMCH\
	indexseqlibs\
	vepe\
	trace_diff\
	trace_clip\
	clip\
	repe\
	lookup

JOB = all

all:	$(DIRS)

old:	$(DIRS_OLD)

clean:
	$(MAKE) $(FLAGS) JOB=clean

spotless:
	$(MAKE) $(FLAGS) JOB=spotless
	-rm lib/$(O)/so_locations

relink:
	-rm -f lib/$(O)/so_locations
	$(MAKE) $(FLAGS) JOB=relink

depend:
	$(MAKE) $(FLAGS) JOB=depend

install:
	$(MAKE) $(FLAGS) JOB=install


.PHONY: distcommon distwindows distunix


# Distribute just enough source to allow third party extensions to build.
# In other words, the header files and perhaps the free bits (like io_lib).
distcommon:
	-mkdir $(DIST)/src

	-cp README.subset $(DIST)/src/README
	-cp -R Makefile $(DIST)/src

	-mkdir $(DIST)/src/Misc
	-cp -R Misc/*.h $(DIST)/src/Misc/

#	-mkdir $(DIST)/src/composition
#	cd composition;\
#	cp -R Makefile *.[ch] *.tcl dependencies $(DIST)/src/composition;\
#	cp -R *.html *.topic *.index $(DIST)/src/composition;\
#	cp -R compositionrc tclIndex $(DIST)/src/composition;
#	cd ..

	-mkdir $(DIST)/src/g
	-cp -R g/*.h $(DIST)/src/g

	-mkdir $(DIST)/src/gap4
	-cp -R gap4/*.h $(DIST)/src/gap4
	-cp -R gap4/*.h $(DIST)/src/gap4

	# FIXME: It's downloadable from our ftpsite anyway.
	#
	# -mkdir $(DIST)/src/io_lib
	# -cp -R io_lib/*.h $(DIST)/src/io_lib

	-mkdir $(DIST)/src/lib
	-cp -R lib/README $(DIST)/src/lib

	-mkdir $(DIST)/src/mk
	-cp -R mk/*.mk mk/README $(DIST)/src/mk

	-mkdir $(DIST)/src/scripts
	cd scripts;\
	cp FE-stealdata FE-stealdata.csh Makefile       $(DIST)/src/scripts;\
	cp assemble assemble4 assemble4.doc bzip2_trace $(DIST)/src/scripts;\
	cp bzip_trace clip-seq.awk clip-seqs            $(DIST)/src/scripts;\
	cp clipping.doc dependencies dub exp2fasta      $(DIST)/src/scripts;\
	cp fasta-split fasta2exp make_fasta gzip_trace  $(DIST)/src/scripts;\
	cp ncbi_gc_parse.tcl splitseq_da staden2exp     $(DIST)/src/scripts;\
	cp stealdata stealdata.doc stealdata.stl        $(DIST)/src/scripts;\
	cd ..

	-mkdir $(DIST)/src/seq_utils
	-cp -R seq_utils/*.h $(DIST)/src/seq_utils

	-mkdir $(DIST)/src/text_utils
	-cp -R text_utils/*.h $(DIST)/src/text_utils

	-mkdir $(DIST)/src/tk_utils
	-cp -R tk_utils/*.h $(DIST)/src/tk_utils

	# Bits of TCL that are needed for compilation of libraries (eg:PRIDE)
	-mkdir $(DIST)/src/tcl8.4.0
	-mkdir $(DIST)/src/tcl8.4.0/generic
	-cp -R tcl8.4.0/generic/tcl*.h $(DIST)/src/tcl8.4.0/generic
	-mkdir $(DIST)/src/tk8.4.0
	-mkdir $(DIST)/src/tk8.4.0/generic
	-mkdir $(DIST)/src/tk8.4.0/compat
	-cp -R tk8.4.0/generic/tk*.h $(DIST)/src/tk8.4.0/generic
	-cp -R tk8.4.0/compat/*.h $(DIST)/src/tk8.4.0/compat



distunix: distcommon
	-mkdir $(DIST)/src/composition/alpha-binaries
	-mkdir $(DIST)/src/composition/linux-binaries
	-mkdir $(DIST)/src/composition/sgi-binaries
	-mkdir $(DIST)/src/composition/solaris-binaries
	-mkdir $(DIST)/src/composition/install
	-mkdir $(DIST)/src/lib/alpha-binaries
	-mkdir $(DIST)/src/lib/linux-binaries
	-mkdir $(DIST)/src/lib/sgi-binaries
	-mkdir $(DIST)/src/lib/solaris-binaries

	# Biolims - for Solaris only
	-cp -R lib/solaris-binaries/libbiolimsio.so $(DIST)/lib/solaris-binaries
	-cp -R lib/solaris-binaries/libbiolimsgap.so $(DIST)/lib/solaris-binaries

	# Bits of TCL that are needed for compilation of libraries (eg:PRIDE)
	-mkdir $(DIST)/src/tcl8.4.0/unix
	-cp -R tcl8.4.0/unix/*.h $(DIST)/src/tcl8.4.0/unix
	-mkdir $(DIST)/src/tk8.4.0/unix
	-cp -R tk8.4.0/unix/*.h $(DIST)/src/tk8.4.0/unix



distwindows: distcommon
	-mkdir $(DIST)/src/composition/windows-binaries
	-cp windows/run/Release/sprun.exe $(DIST)
	-cp trev/trev.bat $(DIST)/windows-bin
	-cp gap4/gap4.bat $(DIST)/windows-bin
	-cp spin/spin.bat $(DIST)/windows-bin
	-cp pregap4/pregap4.bat $(DIST)/windows-bin
	-cp ../windows-bin/*.exe $(DIST)/windows-bin
	-cp scripts/splitseq_da $(DIST)/windows-bin
	-mkdir $(DIST)/src/tcl8.4.0/win
	-cp -R tcl8.4.0/win/*.h $(DIST)/src/tcl8.4.0/win
	-mkdir $(DIST)/src/tk8.4.0/win
	-cp -R tk8.4.0/win/*.h $(DIST)/src/tk8.4.0/win




distsrc:
	-mkdir $(DISTSRC)/src
	-mkdir $(DISTSRC)/src/mk
	-mkdir $(DISTSRC)/src/lib
	-mkdir $(DISTSRC)/src/lib/alpha-binaries
	-mkdir $(DISTSRC)/src/lib/linux-binaries
	-mkdir $(DISTSRC)/src/lib/sgi-binaries
	-mkdir $(DISTSRC)/src/lib/solaris-binaries
	cp -R mk/*.mk mk/README $(DISTSRC)/src/mk
	cp -R Makefile $(DISTSRC)/src
	mkdir $(DISTSRC)/src/ftp
	cp -R ftp/tcl8.4.0sp.tar.gz \
	      ftp/tk8.4.0sp.tar.gz \
	      ftp/cvs-itcl3.1.0.tar.gz \
	      ftp/zlib.tar.gz \
	      ftp/libpng-1.2.5.tar.bz2 \
	      $(DISTSRC)/src/ftp
	$(MAKE) $(FLAGS) JOB=distsrc
	# A few bits of Tcl that are needed for compilation
	-mkdir $(DISTSRC)/src/tcl8.4.0
	-mkdir $(DISTSRC)/src/tcl8.4.0/generic
	-cp -R tcl8.4.0/generic/tcl*.h $(DISTSRC)/src/tcl8.4.0/generic/
	-mkdir $(DISTSRC)/src/tcl8.4.0/unix
	-cp -R tcl8.4.0/unix/tclAppInit.c $(DISTSRC)/src/tcl8.4.0/unix
	# A few bits of Tk that are needed for compilation
	-mkdir $(DISTSRC)/src/tk8.4.0
	-mkdir $(DISTSRC)/src/tk8.4.0/generic
	-mkdir $(DISTSRC)/src/tk8.4.0/unix
	-mkdir $(DISTSRC)/src/tk8.4.0/compat
	-cp -R tk8.4.0/generic/tk*.h $(DISTSRC)/src/tk8.4.0/generic/
	-cp -R tk8.4.0/unix/tkUnixPort.h $(DISTSRC)/src/tk8.4.0/unix
	-cp -R tk8.4.0/compat/*.h $(DISTSRC)/src/tk8.4.0/compat
	# Some itcl stuff
	-mkdir -p $(DISTSRC)/src/incrTcl-3.3cvs/itcl/generic
	-cp -R incrTcl-3.3cvs/itcl/generic/*.h \
	       $(DISTSRC)/src/incrTcl-3.3cvs/itcl/generic
	-mkdir -p $(DISTSRC)/src/incrTcl-3.3cvs/itk/generic
	-cp -R incrTcl-3.3cvs/itk/generic/*.h \
	       $(DISTSRC)/src/incrTcl-3.3cvs/itk/generic
	# Dummy licence code
	-cp -R licence/distrib $(DISTSRC)/src/licence
	## STLPort
	#-mkdir $(DISTSRC)/src/stlport
	#-cp -R stlport/alpha $(DISTSRC)/src/stlport
	#-cp -R stlport/linux $(DISTSRC)/src/stlport
	#-cp -R stlport/macosx $(DISTSRC)/src/stlport

	# Windows compilation bits
	mkdir $(DISTSRC)/src/windows
	cp -R windows/include $(DISTSRC)/src/windows
	mkdir $(DISTSRC)/src/windows/mkdef
	cp -R windows/mkdef/mkdef.* \
	      $(DISTSRC)/src/windows/mkdef
	mkdir $(DISTSRC)/src/windows/run
	cp -R windows/run/run.* \
	      windows/run/main.cpp\
	      windows/run/notes.txt\
	      windows/run/sprun.rc\
	      windows/run/staden.ico\
	      $(DISTSRC)/src/windows/run
	mkdir $(DISTSRC)/src/windows/dumpenv
	cp -R windows/dumpenv/dumpenv.*\
	      windows/dumpenv/main.cpp\
	      $(DISTSRC)/src/windows/dumpenv
	mkdir $(DISTSRC)/src/windows/stashmsg
	cp -R windows/stashmsg/stashmsg.*\
	      windows/stashmsg/icon1.ico\
	      windows/stashmsg/messages.*\
	      windows/stashmsg/msg00001.bin\
	      windows/stashmsg/resource.h\
	      windows/stashmsg/stadenmsg.rc\
	      $(DISTSRC)/src/windows/stashmsg
	cp -R windows/install-2003 $(DISTSRC)/src/windows

$(DIRS): FORCE
	cd $@ && $(MAKE) $(FLAGS) $(MFLAGS) $(JOB)

FORCE:
# DO NOT DELETE THIS LINE -- make depend depends on it.
