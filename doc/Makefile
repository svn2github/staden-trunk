VERSION=2.0.0b1
PREFIX=_inst

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
	-mkdir -p $(PREFIX)/share/man
	cp -R manual/man/* $(PREFIX)/share/man/
	#
	# Main PDF docs
	-mkdir -p $(PREFIX)/share/doc/staden
	cp manual/manual.pdf manual/mini.pdf $(PREFIX)/share/doc/staden
	#
	# HTML pages
	-mkdir $(PREFIX)/share/doc/staden/manual
	cp manual/*.html $(PREFIX)/share/doc/staden/manual
	cp manual/*.png $(PREFIX)/share/doc/staden/manual
	cp manual/*.index $(PREFIX)/share/doc/staden/manual
	cp -r i $(PREFIX)/share/doc/staden/manual
	-mkdir $(PREFIX)/share/doc/staden/scripting_manual
	-cp scripting_manual/*.html $(PREFIX)/share/doc/staden/scripting_manual
	cp -r i $(PREFIX)/share/doc/staden/manual/scripting_manual
	cp *.html $(PREFIX)/share/doc/staden
	#
	# Other bits and pieces
	cp -r Acknowledgements *.txt *.pdf i $(PREFIX)/share/doc/staden
	#
	# Tidy up
	-test "$(PREFIX)" != "." && find $(PREFIX) -name .svn -exec rm -rf {} \;


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

