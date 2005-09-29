.PHONY: distcommon distunix distwindows install

distcommon:
	-mkdir $(DIST)/lib
	# Third party libraries
	cp -R tcl $(DIST)/lib
	cp -R tk $(DIST)/lib
	cp -R itcl3.3 $(DIST)/lib
	cp -R itk3.3 $(DIST)/lib
	cp -R iwidgets $(DIST)/lib
	cp -r tablelist3.4 $(DIST)/lib
	-rm -r $(DIST)/lib/iwidgets/demos
	-rm -r $(DIST)/lib/tk/demos
	-rm -r $(DIST)/lib/tk/orig
	-rm -r $(DIST)/lib/tablelist2.0/CHANGES.txt
	-rm -r $(DIST)/lib/tablelist2.0/README.txt
	-rm -r $(DIST)/lib/tablelist2.0/demos
	-rm -r $(DIST)/lib/tablelist2.0/doc

	@# Tcl interfaces to third party source - we don't want to use the
	@# Makefiles for these in the real source dir.
	-mkdir $(DIST)/lib/cap2
	cp cap2/*.tcl cap2/tclIndex $(DIST)/lib/cap2
	-mkdir $(DIST)/lib/cap3
	cp cap3/*.tcl cap3/tclIndex $(DIST)/lib/cap3
	-mkdir $(DIST)/lib/fak2
	cp fak2/*.tcl fak2/tclIndex $(DIST)/lib/fak2
	-mkdir $(DIST)/lib/phrap
	cp phrap/*.tcl phrap/tclIndex $(DIST)/lib/phrap


distunix: distcommon
	-mkdir $(DIST)/lib/alpha-binaries
	cp -R alpha-binaries/* $(DIST)/lib/alpha-binaries
	-mkdir $(DIST)/lib/linux-binaries
	cp -R linux-binaries/* $(DIST)/lib/linux-binaries
#	cp -R ../lib.third/sgi-binaries $(DIST)/lib
#	cp -R ../lib.third/solaris-binaries $(DIST)/lib
#	cp -R ../lib.third/macosx-binaries $(DIST)/lib

distmacosx: distcommon
	cp -R macosx-binaries $(DIST)/lib

distwindows: distcommon
	mv $(DIST)/lib/windows-binaries/itcl33.dll $(DIST)/lib/itcl3.3
	mv $(DIST)/lib/windows-binaries/itk33.dll $(DIST)/lib/itk3.3

distsrc:
	-mkdir $(DISTSRC)/lib

	cp Makefile $(DISTSRC)/lib

	@# Third party libraries of tcl code.
	cp -R tcl $(DISTSRC)/lib
	cp -R tk $(DISTSRC)/lib
	cp -R itcl3.3 $(DISTSRC)/lib
	cp -R itk3.3 $(DISTSRC)/lib
	cp -R iwidgets $(DISTSRC)/lib
	cp -r tablelist3.4 $(DISTSRC)/lib
	-rm -r $(DISTSRC)/lib/iwidgets/demos
	-rm -r $(DISTSRC)/lib/tk/demos
	-rm -r $(DISTSRC)/lib/tk/orig
	-rm -r $(DISTSRC)/lib/tablelist2.0/CHANGES.txt
	-rm -r $(DISTSRC)/lib/tablelist2.0/README.txt
	-rm -r $(DISTSRC)/lib/tablelist2.0/demos
	-rm -r $(DISTSRC)/lib/tablelist2.0/doc

	@# Tcl interfaces to third party source - we don't want to use the
	@# Makefiles for these in the real source dir.
	-mkdir $(DISTSRC)/lib/cap2
	cp cap2/*.tcl cap2/tclIndex $(DISTSRC)/lib/cap2
	-mkdir $(DISTSRC)/lib/cap3
	cp cap3/*.tcl cap3/tclIndex $(DISTSRC)/lib/cap3
	-mkdir $(DISTSRC)/lib/fak2
	cp fak2/*.tcl fak2/tclIndex $(DISTSRC)/lib/fak2
	-mkdir $(DISTSRC)/lib/phrap
	cp phrap/*.tcl phrap/tclIndex $(DISTSRC)/lib/phrap

	cd $(DISTSRC)/lib; ln -s ../src/lib/linux-binaries .
	cd $(DISTSRC)/lib; ln -s ../src/lib/alpha-binaries .
	cd $(DISTSRC)/lib; ln -s ../src/lib/solaris-binaries .
	cd $(DISTSRC)/lib; ln -s ../src/lib/sgi-binaries .
	cd $(DISTSRC)/lib; ln -s ../src/lib/macosx-binaries .

install:

