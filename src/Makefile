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

xrelink:
	-rm -f lib/$(O)/so_locations
	$(MAKE) $(FLAGS) JOB=relink

depend:
	$(MAKE) $(FLAGS) JOB=depend

# Installs newly built executables into $STADENROOT/$(MACHINE)-bin
install:
	$(MAKE) $(FLAGS) JOB=install


.PHONY: distcommon distwindows distunix


# Full source is available now, so we no longer need to copy any of this
# when building a binary distribution.

# Make install does most of the work, but this doesn't install the third
# party libraries. So we'll do that here. NOTE: this only installs the
# binary components. The tcl/itcl files are already in the ../lib directory
# and will be copied over separately. This is because in some cases they will
# have been manually hacked/bug fixed.
distunix:
	-cp -R lib/$(MACHINE)-binaries/libpng* $(DIST)/lib/$(MACHINE)-binaries
	-cp -R lib/$(MACHINE)-binaries/libtcl* $(DIST)/lib/$(MACHINE)-binaries
	-cp -R lib/$(MACHINE)-binaries/libtk* $(DIST)/lib/$(MACHINE)-binaries
	-cp -R lib/$(MACHINE)-binaries/libitcl* $(DIST)/lib/$(MACHINE)-binaries
	-cp -R lib/$(MACHINE)-binaries/libitk* $(DIST)/lib/$(MACHINE)-binaries


distwindows:
	-mkdir $(DIST)/src/composition/windows-binaries
	-cp windows/run/Release/sprun.exe $(DIST)
	-cp trev/trev.bat $(DIST)/windows-bin
	-cp gap4/gap4.bat $(DIST)/windows-bin
	-cp spin/spin.bat $(DIST)/windows-bin
	-cp pregap4/pregap4.bat $(DIST)/windows-bin
	-cp ../windows-bin/*.exe $(DIST)/windows-bin
	-cp scripts/splitseq_da $(DIST)/windows-bin




distsrc:
	-mkdir $(DISTSRC)/src
	-mkdir $(DISTSRC)/src/mk
	-mkdir $(DISTSRC)/src/lib
	-mkdir $(DISTSRC)/src/lib/alpha-binaries
	-mkdir $(DISTSRC)/src/lib/linux-binaries
	-mkdir $(DISTSRC)/src/lib/sgi-binaries
	-mkdir $(DISTSRC)/src/lib/solaris-binaries
	cp -R mk/*.mk mk/README $(DISTSRC)/src/mk
	cp -R Makefile Makefile.thirdparty $(DISTSRC)/src
	$(MAKE) $(FLAGS) JOB=distsrc
#	 A few bits of Tcl that are needed for compilation
#	-mkdir $(DISTSRC)/src/tcl8.4.6
#	-mkdir $(DISTSRC)/src/tcl8.4.6/generic
#	-cp -R tcl8.4.6/generic/tcl*.h $(DISTSRC)/src/tcl8.4.6/generic/
#	-mkdir $(DISTSRC)/src/tcl8.4.6/unix
#	-cp -R tcl8.4.6/unix/tclAppInit.c $(DISTSRC)/src/tcl8.4.6/unix
#	 A few bits of Tk that are needed for compilation
#	-mkdir $(DISTSRC)/src/tk8.4.6
#	-mkdir $(DISTSRC)/src/tk8.4.6/generic
#	-mkdir $(DISTSRC)/src/tk8.4.6/unix
#	-mkdir $(DISTSRC)/src/tk8.4.6/compat
#	-cp -R tk8.4.6/generic/tk*.h $(DISTSRC)/src/tk8.4.6/generic/
#	-cp -R tk8.4.6/unix/tkUnixPort.h $(DISTSRC)/src/tk8.4.6/unix
#	-cp -R tk8.4.6/compat/*.h $(DISTSRC)/src/tk8.4.6/compat
#	 Some itcl stuff
#	-mkdir -p $(DISTSRC)/src/incrTcl-3.3cvs/itcl/generic
#	-cp -R incrTcl-3.3cvs/itcl/generic/*.h \
#	       $(DISTSRC)/src/incrTcl-3.3cvs/itcl/generic
#	-mkdir -p $(DISTSRC)/src/incrTcl-3.3cvs/itk/generic
#	-cp -R incrTcl-3.3cvs/itk/generic/*.h \
#	       $(DISTSRC)/src/incrTcl-3.3cvs/itk/generic
#	 Dummy licence code
	-cp -R licence $(DISTSRC)/src/licence
#	# STLPort
#	-mkdir $(DISTSRC)/src/stlport
#	-cp -R stlport/alpha $(DISTSRC)/src/stlport
#	-cp -R stlport/linux $(DISTSRC)/src/stlport
#	-cp -R stlport/macosx $(DISTSRC)/src/stlport

#	 Windows compilation bits
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
