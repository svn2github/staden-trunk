all:
	cd texi; $(MAKE) $(SUBFLAGS) spotless all
	cd mutation_texi; $(MAKE) $(SUBFLAGS) spotless all

.PHONY: distcommon distunix distwindows install

distcommon:
	-mkdir $(DIST)/course
	-cp -R README $(DIST)/course
	-cp -R data $(DIST)/course


distunix: distcommon
	mkdir $(DIST)/course/unix_docs
	-cp -R texi/course_unix.pdf $(DIST)/course/unix_docs
	-cp -R texi/course_unix.ps $(DIST)/course/unix_docs
	-cp -R mutation_texi/notes.ps $(DIST)/course/unix_docs/mutation_notes.ps
	-cp -R mutation_texi/notes.pdf $(DIST)/course/unix_docs/mutation_notes.pdf


distwindows: distcommon
	mkdir $(DIST)/course/windows_docs
	-cp -R texi/course_windows.pdf $(DIST)/course/windows_docs
	-cp -R texi/course_windows.ps $(DIST)/course/windows_docs
	-cp -R mutation_texi/notes.ps $(DIST)/course/windows_docs/mutation_notes.ps
	-cp -R mutation_texi/notes.pdf $(DIST)/course/windows_docs/mutation_notes.pdf

distsrc:
	-mkdir $(DISTSRC)/course
	-cp -R README Makefile $(DISTSRC)/course
	-cp -R data $(DISTSRC)/course
	-mkdir $(DISTSRC)/course/texi
	-cp -R texi/*.texi \
	       texi/Makefile \
	       texi/*.eps \
	       texi/*.gif \
	       texi/*.m4 \
	       texi/Makefile \
	       $(DISTSRC)/course/texi
	-mkdir $(DISTSRC)/course/mutation_texi
	-cp -R mutation_texi/*.texi \
	       mutation_texi/*.eps \
	       mutation_texi/*.gif \
	       mutation_texi/*.m4 \
	       mutation_texi/Makefile \
	       $(DISTSRC)/course/mutation_texi
	-mkdir $(DISTSRC)/course/OpenOffice
	-cp -R OpenOffice/*.sxi $(DISTSRC)/course/OpenOffice
	-cp -R Makefile $(DISTSRC)/course/OpenOffice
	-mkdir $(DISTSRC)/course/unix_docs
	-mkdir $(DISTSRC)/course/windows_docs

install:

