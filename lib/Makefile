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
	-mkdir $(DIST)/lib/gap
	cp ../lib/gap/*.tcl $(DIST)/lib/gap
	cp ../lib/gap/tclIndex $(DIST)/lib/gap
	-mkdir $(DIST)/lib/pregap4
	-mkdir $(DIST)/lib/pregap4/modules
	-mkdir $(DIST)/lib/pregap4/templates
	-mkdir $(DIST)/lib/pregap4/naming_schemes
	cp ../lib/pregap4/widget_licence $(DIST)/lib/pregap4
	cp ../lib/pregap4/*.tcl $(DIST)/lib/pregap4
	cp ../lib/pregap4/pregap4 $(DIST)/lib/pregap4
	cp ../lib/pregap4/tclIndex $(DIST)/lib/pregap4
	cp ../lib/pregap4/modules/*.p4m $(DIST)/lib/pregap4/modules
	cp ../lib/pregap4/naming_schemes/*.p4t $(DIST)/lib/pregap4/naming_schemes
	-mkdir $(DIST)/lib/spin
	cp ../lib/spin/*.tcl $(DIST)/lib/spin
	cp ../lib/spin/*.itk $(DIST)/lib/spin
	cp ../lib/spin/tclIndex $(DIST)/lib/spin
	-mkdir $(DIST)/lib/spin2
	cp ../lib/spin2/*.tcl $(DIST)/lib/spin2
	cp ../lib/spin2/*.itk $(DIST)/lib/spin2
	cp ../lib/spin2/tclIndex $(DIST)/lib/spin2
	-mkdir $(DIST)/lib/seqed
	cp ../lib/seqed/*.tcl $(DIST)/lib/seqed
	cp ../lib/seqed/*.itcl $(DIST)/lib/seqed
	cp ../lib/seqed/*.itk $(DIST)/lib/seqed
	cp ../lib/seqed/*.gif $(DIST)/lib/seqed
	cp ../lib/seqed/tclIndex $(DIST)/lib/seqed
	-mkdir $(DIST)/lib/spin_emboss
	cp -r ../lib/spin_emboss/acdtcl $(DIST)/lib/spin_emboss
	cp ../lib/spin_emboss/*.tcl $(DIST)/lib/spin_emboss
	cp ../lib/spin_emboss/tclIndex $(DIST)/lib/spin_emboss
	-mkdir $(DIST)/lib/spin2_emboss
	cp -r ../lib/spin2_emboss/acdtcl $(DIST)/lib/spin2_emboss
	cp ../lib/spin2_emboss/*.tcl $(DIST)/lib/spin2_emboss
	cp ../lib/spin2_emboss/tclIndex $(DIST)/lib/spin2_emboss
	-mkdir $(DIST)/lib/tk_utils
	cp ../lib/tk_utils/*.tcl $(DIST)/lib/tk_utils
	cp ../lib/tk_utils/tclIndex $(DIST)/lib/tk_utils
	-mkdir $(DIST)/lib/trev
	cp ../lib/trev/*.tcl $(DIST)/lib/trev
	-mkdir $(DIST)/lib/prefinish
	-mkdir $(DIST)/lib/prefinish/rules
	cp ../lib/prefinish/*.tcl $(DIST)/lib/prefinish
	cp ../lib/prefinish/*.itcl $(DIST)/lib/prefinish
	cp ../lib/prefinish/*.itk $(DIST)/lib/prefinish
	cp ../lib/prefinish/tclIndex $(DIST)/lib/prefinish
	cp ../lib/prefinish/config.prefinish $(DIST)/lib/prefinish
	cp ../lib/prefinish/args.template $(DIST)/lib/prefinish
	cp ../lib/prefinish/help $(DIST)/lib/prefinish
	cp ../lib/prefinish/README $(DIST)/lib/prefinish
	cp ../lib/prefinish/CHEMISTRY $(DIST)/lib/prefinish
	cp ../lib/prefinish/METHODS $(DIST)/lib/prefinish
	cp ../lib/prefinish/rules/*.rule $(DIST)/lib/prefinish/rules
	-mkdir $(DIST)/lib/copy_reads
	cp ../lib/copy_reads/*.tcl $(DIST)/lib/copy_reads
	cp ../lib/copy_reads/tclIndex $(DIST)/lib/copy_reads
	-mkdir $(DIST)/lib/windows-binaries
	cp ../lib/windows-binaries/*.dll $(DIST)/lib/windows-binaries

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

