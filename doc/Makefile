all:
	cd scripting_manual; $(MAKE) $(SUBFLAGS)
	cd manual; $(MAKE) $(SUBFLAGS)

.PHONY: distcommon distwindows distunix install

distcommon:
	-mkdir $(DIST)/doc
	-mkdir $(DIST)/doc/manual
	-mkdir $(DIST)/doc/scripting_manual
	-cp -R i $(DIST)/doc/i
	-cp -R manual/i $(DIST)/doc/manual/i
	-cp -R scripting_manual/i $(DIST)/doc/scripting_manual/i
	-cp scripting_manual/*.html $(DIST)/doc/scripting_manual
	-cp Acknowledgements $(DIST)/doc
	-cp Copyright $(DIST)/doc
	-cp LICENCE $(DIST)/doc
	-cp README $(DIST)/doc
	-cp gkb547_gml.pdf $(DIST)/doc
	-cp scripting_manual/*.pdf $(DIST)/doc
	-cp emboss.txt $(DIST)/doc
	-cp *.gif $(DIST)/doc



distwindows: distcommon
	-cp windows/*.html $(DIST)/doc
	-cp manual/*_windows.index $(DIST)/doc/manual
	-cp manual/*_windows_*.html $(DIST)/doc/manual
	-cp manual/manual_windows.pdf $(DIST)/doc
	-cp manual/mini_windows.pdf $(DIST)/doc
	-cp manual/*.unix.png $(DIST)/doc/manual
	-cp manual/*.unix.png.html $(DIST)/doc/manual



distunix: distcommon
	-cp unix/*.html $(DIST)/doc
	-cp manual/*_unix.index $(DIST)/doc/manual
	-cp manual/*_unix_*.html $(DIST)/doc/manual
	-cp manual/manual_unix.pdf $(DIST)/doc
	-cp manual/mini_unix.pdf $(DIST)/doc
	-cp manual/*.unix.png $(DIST)/doc/manual
	-cp manual/*.unix.png.html $(DIST)/doc/manual

distsrc:
	-mkdir -p $(DISTSRC)/doc/manual
	-mkdir -p $(DISTSRC)/doc/scripting_manual
	-cp -R i $(DISTSRC)/doc/i
	-cp -R manual/i \
	       manual/Makefile \
	       manual/*.texi \
	       manual/*.png \
	       manual/*.m4 \
	       manual/README \
	       manual/dependencies \
	       manual/data \
	       manual/tools \
	       manual/man \
	       $(DISTSRC)/doc/manual
	-cp -R scripting_manual/i \
	       scripting_manual/*.texi \
	       scripting_manual/*.m4 \
	       scripting_manual/Makefile \
	       scripting_manual/dependencies \
	       scripting_manual/i \
	       scripting_manual/tools \
	       $(DISTSRC)/doc/scripting_manual
	-cp Acknowledgements $(DISTSRC)/doc
	-cp Copyright $(DISTSRC)/doc
	-cp LICENCE $(DISTSRC)/doc
	-cp README $(DISTSRC)/doc
	-cp gkb547_gml.pdf $(DISTSRC)/doc
	-cp emboss.txt $(DISTSRC)/doc
	-cp *.gif $(DISTSRC)/doc
	-cp Makefile $(DISTSRC)/doc

install:

