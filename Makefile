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
	-mkdir -p $(DIST)
	cp Version $(DIST)
	cd tables; $(MAKE) $(SUBFLAGS) dist
	cd userdata; $(MAKE) $(SUBFLAGS) dist
	cd demo; $(MAKE) $(SUBFLAGS) dist
	cd lib; $(MAKE) $(SUBFLAGS) distcommon


distunix: distcommon
	cp staden.login staden.profile $(DIST)

	# cp install $(DIST)

	cd course; $(MAKE) $(SUBFLAGS) distunix
	cd doc; $(MAKE) $(SUBFLAGS) distunix
	#cd man; $(MAKE) $(SUBFLAGS) distunix
	-mkdir -p $(DIST)/lib/linux-binaries
	cd src; $(MAKE) $(SUBFLAGS) distunix

	#-----------------------------------------------------
	# This is inefficient as it install the lib Tcl files multiple times,
	# but we can live with that as it provides simplicity.
	#-----------------------------------------------------
#	-mkdir $(DIST)/alpha-bin
#	-mkdir -p $(DIST)/lib/alpha-bin
#	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=alpha install
	#-----------------------------------------------------
	-mkdir $(DIST)/linux-bin
	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=linux install
	#-----------------------------------------------------
#	-mkdir $(DIST)/solaris-bin
#	-mkdir -p $(DIST)/lib/solaris-bin
#	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=solaris install
	#-----------------------------------------------------
#	-mkdir $(DIST)/sgi-bin
#	-mkdir -p $(DIST)/lib/sgi-bin
#	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=sgi install
	#-----------------------------------------------------
#	-mkdir $(DIST)/macosx-bin
#	-mkdir -p $(DIST)/lib/macosx-bin
#	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=macosx install


distmacosx: distcommon
	cp staden.login staden.profile $(DIST)

	# cp install $(DIST)
	cp README $(DIST)

	cd course; $(MAKE) $(SUBFLAGS) distunix
	cd doc; $(MAKE) $(SUBFLAGS) distunix
	cd man; $(MAKE) $(SUBFLAGS) distunix
	cd src; $(MAKE) $(SUBFLAGS) distunix

	-mkdir $(DIST)/macosx-bin
	cd src; $(MAKE) STADENROOT=$(DIST) MACHINE=macosx install

	cd src/macosx; $(MAKE) $(SUBFLAGS) MACHINE=macosx install


distwindows: distcommon
	-mkdir $(DIST)/windows-bin
	cd course; $(MAKE) $(SUBFLAGS) distwindows
	cd doc; $(MAKE) $(SUBFLAGS) distwindows
	cd src; $(MAKE) $(SUBFLAGS) distwindows

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
