all:
	cd texi; $(MAKE) $(SUBFLAGS) spotless all
	cd mutation_texi; $(MAKE) $(SUBFLAGS) spotless all

.PHONY: distcommon distunix distwindows install

distcommon:
	-mkdir $(DIST)/course
	-cp -R README $(DIST)/course
	-cp -R data $(DIST)/course


distunix: distcommon
	-cp -R unix_docs $(DIST)/course


distwindows: distcommon
	-cp -R windows_docs $(DIST)/course

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

