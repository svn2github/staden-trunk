all: start

#
# Sorry if this Makefile doesn't work correctly regarding dependencies. GNU
# make causes all sorts of headaches with it's inbuilt rules (which I seem
# unable to remove, even when using -d) and it has some quirky ideas as to
# which files are created temporarily (and thus should be removed). It's
# best to usually do 'gmake spotless all' or some such.
#
# However you should try "gmake depend" to keep the dependencies file up to
# date as this does solve many (if not all) dependency problems.
#
# The input files are always .texi
#
# .texinfo files are expanded up .texi files. They have the macros replaced
# and have been processed by m4 to include Unix or Windows specific components.
#
# The _us.* files are US-letter format copies, otherwise we assume everything
# is in A4 format.
#

# M4 preprocessor. Various buggy versions of this have caused problems in the
# past, so you may need to redefine this. On Digital Unix 4.0E the system m4
# does not work with our files. Certain versions (which?) of GNU m4 also fail,
# but this has now been patched.
M4=m4

#-----------------------------------------------------------------------------
# General rules

#
# M4 processed texinfo
#
%.texinfo:	%.texi
	$(M4) -D_tex < $< > $@
	../manual/tools/update-nodes $@

%.htmlinfo:	%.texi
	$(M4) -D_html < $< > $@
	../manual/tools/update-nodes $@

# How to build .dvi files from our m4-expanded .texinfo
# files
%.dvi:	%.texinfo
	texi2dvi $<

# US-letter PostScript from DVI
%_us.ps:	%_us.dvi
	dvips -t letter -o $@ $<

# A4 PostScript from DVI
%.ps:	%.dvi
	dvips -t a4 -Ppdf -o $@ $<

# PDF generation. I've never really had a great deal of luck with getting
# good fonts out of ps2pdf, but it's worth a try.
%.pdf:	%.ps
	ps2pdf -dPDFSETTINGS=/prepress -dEncodeColorImages=false $<

# Converts an A4 formatted .texi file into a US-letter formatted file.
%_us.texinfo: %.texinfo
	egrep -v '^@afourpaper' < $< > $@

# HTML files - built from an expanded .texinfo file with the -D_html m4 macro
# defined. We need the *_toc.html and the index files.
# For ease of browsing we create a separate html document for each of the main
# programs. The htmlinfo version is identical to texinfo except with a few
# tweaks to the cross-references (to allow cross-referencing between top-level
# documents) and the addition of an _split() command to request splitting an
# html page at a specific point).
%_toc.html:	%.htmlinfo
	../manual/tools/texi2html -menu -verbose -split_chapter -index_chars $<

# Man pages - taken from the .texi files directly (NB: may need m4 expansion,
# but this is not yet applied.)
man/man1/%.1:	%.1.texi
	../manual/tools/texi2man.pl $< > $@

man/man4/ExperimentFile.4: exp-t.texi
	../manual/tools/texi2man.pl $< > $@

man/man4/scf.4: scf-t.texi
	../manual/tools/texi2man.pl $< > $@

man/man4/ztr.4: ztr-t.texi
	../manual/tools/texi2man.pl $< > $@

# For any large pictures (_lpicure() macro) we create an html page containing
# the full-size picture to link from the small picture embedded in the page.
gifs:
	../manual/tools/make_gif_html

# PostScript versions of the gif images.
%.ps: %.gif
	../manual/tools/make_ps $<

# Help within the programs is HTML based. Given a topic the appropriate URL
# is obtained by looking it up in the .index file.
%.index:	%_toc.html
	../manual/tools/html_index.pl $<

#-----------------------------------------------------------------------------
# The main make targets.

clean:
	-rm *.aux *.cp *.fn *.ky *.log *.pg *.toc *.tp *.vr *.cps *.fns *.pgs *.vrs
	-rm core _tmp.texi _tmp.texi~ *.texinfo *.texinfo.tmp *.texinfo~
	-rm *.htmlinfo *.htmlinfo~

spotless:	clean
	-rm *.dvi *.html *.info *.info-[0-9] *.index *.topic
	-rm *.ps

depend:
	../manual/tools/make_dependencies > dependencies

include dependencies

start: scripting.texinfo scripting.dvi scripting.ps scripting_toc.html
