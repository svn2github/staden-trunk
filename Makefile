#-------------
# Directories
#-------------

VERSION=$(shell cat Version)
#VERSION=2002.0a3
DISTROOT=$(shell pwd)/distrib

ifeq ($(MACHINE),windows)
DIST=$(DISTROOT)/windows-$(VERSION)
else
ifeq ($(MACHINE),macosx)
DIST=$(DISTROOT)/macosx-$(VERSION)
else
DIST=$(DISTROOT)/unix-$(VERSION)
endif
endif
DISTV=$(DIST)v

DISTSRC=$(DIST)-src
SUBFLAGS=DIST=$(DIST)



#---------
# Targets
#---------

.PHONY: all install dist distsrc depend distwindows distunix distmacosx distcommon

all:
	cd src; $(MAKE) all
	cd doc; $(MAKE) all
	cd course; $(MAKE) all

install:
	cd src; $(MAKE) $@

clean:
	cd src; $(MAKE) clean

depend:
	cd src; $(MAKE) depend




#----------------------
# Distribution Targets
#----------------------

distmsg:
	@echo Ready to copy distribution to $(DISTROOT)/$(VERSION)
	@echo Hit control-C now to abort.
	@sleep 3

dist: distmsg distunix distwindows distmacosx
	echo Distribution built.


distcommon:
	-mkdir $(DIST)
	cp Version-* $(DIST)
	cd tables; $(MAKE) $(SUBFLAGS) dist
	cd userdata; $(MAKE) $(SUBFLAGS) dist
	cd demo; $(MAKE) $(SUBFLAGS) dist


distunix: distcommon
	cp staden.login staden.profile $(DIST)

	# cp install $(DIST)

	cd course; $(MAKE) $(SUBFLAGS) distunix
	cd doc; $(MAKE) $(SUBFLAGS) distunix
	cd lib; $(MAKE) $(SUBFLAGS) distunix
	#cd man; $(MAKE) $(SUBFLAGS) distunix
	cd src; $(MAKE) $(SUBFLAGS) distunix

	#-----------------------------------------------------
	# This is inefficient as it install the lib Tcl files multiple times,
	# but we can live with that as it provides simplicity.
	#-----------------------------------------------------
	-mkdir $(DIST)/alpha-bin
	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=alpha install
	#-----------------------------------------------------
#	-mkdir $(DIST)/linux-bin
#	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=linux install
	#-----------------------------------------------------
#	-mkdir $(DIST)/solaris-bin
#	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=solaris install
	#-----------------------------------------------------
#	-mkdir $(DIST)/sgi-bin
#	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=sgi install
	#-----------------------------------------------------
#	-mkdir $(DIST)/macosx-bin
#	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=macosx install


distmacosx: distcommon
	cp staden.login staden.profile $(DIST)

	# cp install $(DIST)
	cp README $(DIST)

	cd course; $(MAKE) $(SUBFLAGS) distunix
	cd doc; $(MAKE) $(SUBFLAGS) distunix
	cd lib; $(MAKE) $(SUBFLAGS) distmacosx
	cd man; $(MAKE) $(SUBFLAGS) distunix
	cd src; $(MAKE) $(SUBFLAGS) distunix

	-mkdir $(DIST)/macosx-bin
	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=macosx install

	cd src/macosx; $(MAKE) $(SUBFLAGS) MACHINE=macosx install


distwindows: distcommon
	-mkdir $(DIST)/windows-bin
	cd course; $(MAKE) $(SUBFLAGS) distwindows
	cd doc; $(MAKE) $(SUBFLAGS) distwindows
	cd lib; $(MAKE) $(SUBFLAGS) distwindows
	cd src; $(MAKE) $(SUBFLAGS) distwindows

#
# Produce a cutdown viewer version from an existing full distribution
#
distviewercommon:
	-mkdir $(DISTV)
	-cp $(DIST)/Version-* $(DISTV)

	-mkdir $(DISTV)/lib
	-cp -R $(DIST)/lib/gap       $(DISTV)/lib
	-cp -R $(DIST)/lib/tcl       $(DISTV)/lib
	-cp -R $(DIST)/lib/tk        $(DISTV)/lib
	-cp -R $(DIST)/lib/tk_utils  $(DISTV)/lib
	-cp -R $(DIST)/lib/trev      $(DISTV)/lib

	-mkdir $(DISTV)/doc
	-cp -R $(DIST)/doc/Acknowledgements \
	       $(DIST)/doc/Copyright \
	       $(DIST)/doc/LICENCE \
	       $(DIST)/doc/manual \
	       $(DISTV)/doc
	-rm $(DISTV)/doc/manual/*.ps
	-cp doc/ReleaseNotes.viewer $(DISTV)/doc/ReleaseNotes

	-mkdir $(DISTV)/tables
	-cp -R $(DIST)/tables/GTAGDB \
	       $(DIST)/tables/NOTEDB \
	       $(DIST)/tables/RENZYM.* \
	       $(DIST)/tables/follow.bitmap \
	       $(DIST)/tables/gaprc \
	       $(DIST)/tables/gaprc_menu_viewer \
	       $(DIST)/tables/gcodes \
	       $(DIST)/tables/help_config \
	       $(DIST)/tables/nuc_matrix \
	       $(DIST)/tables/shlib.conf \
	       $(DIST)/tables/tk_utilsrc \
	       $(DIST)/tables/trevrc \
	       $(DISTV)/tables

	cp -R viewer_example $(DISTV)/example

distviewerunix: distviewercommon
	-cp doc/INSTALL.viewer $(DISTV)/INSTALL

	cp tables/licence.viewer.unix $(DISTV)/tables/licence.txt
	cp $(DIST)/staden.profile $(DIST)/staden.login $(DISTV)

	-mkdir $(DISTV)/lib/alpha-binaries
	-mkdir $(DISTV)/alpha-bin
	-cp $(DIST)/alpha-bin/stash \
	    $(DIST)/alpha-bin/gap4 \
	    $(DIST)/alpha-bin/trev \
	    $(DISTV)/alpha-bin
	-cp $(DIST)/lib/alpha-binaries/libg.so* \
	    $(DIST)/lib/alpha-binaries/libgap.so* \
	    $(DIST)/lib/alpha-binaries/libio-utils.so* \
	    $(DIST)/lib/alpha-binaries/libmisc.so* \
	    $(DIST)/lib/alpha-binaries/libread.so* \
	    $(DIST)/lib/alpha-binaries/libseq_utils.so* \
	    $(DIST)/lib/alpha-binaries/libtcl8.3.so* \
	    $(DIST)/lib/alpha-binaries/libtk8.3.so* \
	    $(DIST)/lib/alpha-binaries/libtk_utils.so* \
	    $(DISTV)/lib/alpha-binaries

	-mkdir $(DISTV)/lib/linux-binaries
	-mkdir $(DISTV)/linux-bin
	-cp $(DIST)/linux-bin/stash \
	    $(DIST)/linux-bin/gap4 \
	    $(DIST)/linux-bin/trev \
	    $(DISTV)/linux-bin
	-cp $(DIST)/lib/linux-binaries/libg.so* \
	    $(DIST)/lib/linux-binaries/libgap.so* \
	    $(DIST)/lib/linux-binaries/libio-utils.so* \
	    $(DIST)/lib/linux-binaries/libmisc.so* \
	    $(DIST)/lib/linux-binaries/libread.so* \
	    $(DIST)/lib/linux-binaries/libseq_utils.so* \
	    $(DIST)/lib/linux-binaries/libtcl8.3.so* \
	    $(DIST)/lib/linux-binaries/libtk8.3.so* \
	    $(DIST)/lib/linux-binaries/libtk_utils.so* \
	    $(DISTV)/lib/linux-binaries

	-mkdir $(DISTV)/lib/sgi-binaries
	-mkdir $(DISTV)/sgi-bin
	-cp $(DIST)/sgi-bin/stash \
	    $(DIST)/sgi-bin/gap4 \
	    $(DIST)/sgi-bin/trev \
	    $(DISTV)/sgi-bin
	-cp $(DIST)/lib/sgi-binaries/libg.so* \
	    $(DIST)/lib/sgi-binaries/libgap.so* \
	    $(DIST)/lib/sgi-binaries/libio-utils.so* \
	    $(DIST)/lib/sgi-binaries/libmisc.so* \
	    $(DIST)/lib/sgi-binaries/libread.so* \
	    $(DIST)/lib/sgi-binaries/libseq_utils.so* \
	    $(DIST)/lib/sgi-binaries/libtcl8.3.so* \
	    $(DIST)/lib/sgi-binaries/libtk8.3.so* \
	    $(DIST)/lib/sgi-binaries/libtk_utils.so* \
	    $(DISTV)/lib/sgi-binaries

	-mkdir $(DISTV)/lib/solaris-binaries
	-mkdir $(DISTV)/solaris-bin
	-cp $(DIST)/solaris-bin/stash \
	    $(DIST)/solaris-bin/gap4 \
	    $(DIST)/solaris-bin/trev \
	    $(DISTV)/solaris-bin
	-cp $(DIST)/lib/solaris-binaries/libg.so* \
	    $(DIST)/lib/solaris-binaries/libgap.so* \
	    $(DIST)/lib/solaris-binaries/libio-utils.so* \
	    $(DIST)/lib/solaris-binaries/libmisc.so* \
	    $(DIST)/lib/solaris-binaries/libread.so* \
	    $(DIST)/lib/solaris-binaries/libseq_utils.so* \
	    $(DIST)/lib/solaris-binaries/libtcl8.3.so* \
	    $(DIST)/lib/solaris-binaries/libtk8.3.so* \
	    $(DIST)/lib/solaris-binaries/libtk_utils.so* \
	    $(DISTV)/lib/solaris-binaries

distviewerwindows: distviewercommon
	cp tables/licence.viewer.windows $(DISTV)/tables/licence.txt

	-mkdir $(DISTV)/lib/windows-binaries
	-mkdir $(DISTV)/windows-bin

	-cp $(DIST)/windows-bin/winstash.exe $(DISTV)/windows-bin
	-cp src/windows/run/run.exe $(DISTV)

	-cp $(DIST)/lib/windows-binaries/g.dll \
	    $(DIST)/lib/windows-binaries/gap.dll \
	    $(DIST)/lib/windows-binaries/io-utils.dll \
	    $(DIST)/lib/windows-binaries/misc.dll \
	    $(DIST)/lib/windows-binaries/msvcrt.dll \
	    $(DIST)/lib/windows-binaries/read.dll \
	    $(DIST)/lib/windows-binaries/seq_utils.dll \
	    $(DIST)/lib/windows-binaries/stashmsg.dll \
	    $(DIST)/lib/windows-binaries/tcl83.dll \
	    $(DIST)/lib/windows-binaries/tk83.dll \
	    $(DIST)/lib/windows-binaries/tk_utils.dll \
	    $(DIST)/lib/windows-binaries/z.dll \
	    $(DISTV)/lib/windows-binaries

#
# Creates a distribution of the source code suitable for sending to third
# parties.
#
distsrc:
	# Strictly source stuff
	-mkdir $(DISTSRC)
	-cp Version* Makefile $(DISTSRC)
	cd src; $(MAKE) DISTSRC=$(DISTSRC) $@
	cp -R README.build $(DISTSRC)

	# Other things needed to run the package - tables, configs, docs, etc
	cp staden.profile staden.login $(DISTSRC)
	cd tables; $(MAKE) DISTSRC=$(DISTSRC) $@
	cd userdata; $(MAKE) DISTSRC=$(DISTSRC) $@
	cd demo; $(MAKE) DISTSRC=$(DISTSRC) $@
	cd course; $(MAKE) DISTSRC=$(DISTSRC) $@
	cd doc; $(MAKE) DISTSRC=$(DISTSRC) $@
	cd lib; $(MAKE) DISTSRC=$(DISTSRC) $@
