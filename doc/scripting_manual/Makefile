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

# PDF generation - directly from the texinfo file.
%.pdf:	%.texinfo
	texi2pdf $<


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

start: scripting.texinfo scripting.pdf scripting_toc.html
