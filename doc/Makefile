VERSION=2.0.0b8
prefix=/usr

# GNU standard dir names
datarootdir = ${prefix}/share
docdir      = ${datarootdir}/doc/staden
mandir      = ${datarootdir}/man
man1dir     = ${mandir}/man1
man4dir     = ${mandir}/man4
htmldir     = ${docdir}

# Dir names used below, incorporating DESTDIR
DOCDIR     = ${DESTDIR}${docdir}
MANDIR     = ${DESTDIR}${mandir}
MAN1DIR    = ${DESTDIR}${man1dir}
MAN4DIR    = ${DESTDIR}${man4dir}
HTMLDIR    = ${DESTDIR}${htmldir}

all:
	@echo
	@echo Please rerun make specifying either target \"unix\" or \"windows\".
	@echo

unix: SYSTEM=unix
unix: common

windows: SYSTEM=windows
windows: common

common:
	cd manual; $(MAKE) $(SUBFLAGS) $(SYSTEM)
	cd scripting_manual; $(MAKE) $(SUBFLAGS)
	./parse_template $(SYSTEM) < overview.html.template > overview.html

install:
	# Man pages
	-mkdir -p            $(MAN1DIR)
	cp manual/man/man1/* $(MAN1DIR)

	-mkdir -p            $(MAN4DIR)
	cp manual/man/man4/* $(MAN4DIR)

	# Main PDF docs
	-mkdir -p                            $(DOCDIR)
	cp manual/manual.pdf manual/mini.pdf $(DOCDIR)

	# HTML pages
	cp *.html $(HTMLDIR)

	-mkdir            $(HTMLDIR)/manual
	cp manual/*.html  $(HTMLDIR)/manual
	cp manual/*.png   $(HTMLDIR)/manual
	cp manual/*.index $(HTMLDIR)/manual
	-mkdir            $(HTMLDIR)/manual/i
	cp i/*            $(HTMLDIR)/manual/i

	-mkdir                      $(HTMLDIR)/scripting_manual
	-cp scripting_manual/*.html $(HTMLDIR)/scripting_manual
	-mkdir                      $(HTMLDIR)/scripting_manual/i
	cp i/*                      $(HTMLDIR)/scripting_manual/i

	# Other bits and pieces
	cp Acknowledgements *.txt *.pdf $(DOCDIR)
	-mkdir -p $(DOCDIR)/i
	cp i/*    $(DOCDIR)/i


DISTSRC=staden_doc-$(VERSION)-src
distsrc:
	-mkdir -p $(DISTSRC)/manual
	-mkdir -p $(DISTSRC)/scripting_manual
	-cp -R i $(DISTSRC)/i
	-cp -R manual/Makefile \
	       manual/*.texi \
	       manual/*.png \
	       manual/*.m4 \
	       manual/README \
	       manual/dependencies \
	       manual/tools \
	       manual/man \
	       $(DISTSRC)/manual
	-cp -R scripting_manual/*.texi \
	       scripting_manual/*.m4 \
	       scripting_manual/Makefile \
	       scripting_manual/dependencies \
	       scripting_manual/i \
	       scripting_manual/tools \
	       $(DISTSRC)/scripting_manual
	-cp Acknowledgements $(DISTSRC)
	-cp README $(DISTSRC)
	-cp *.pdf $(DISTSRC)
	-cp emboss.txt $(DISTSRC)
	-cp *.gif $(DISTSRC)
	-cp Makefile $(DISTSRC)
	-cp *template $(DISTSRC)
	-find $(DISTSRC) -name .svn -exec rm -rf {} \;
	tar cfz $(DISTSRC).tar.gz $(DISTSRC)

install:

