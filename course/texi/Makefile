#
# Sorry if this Makefile doesn't work correctly regarding dependencies. GNU
# make causes all sorts of headaches with it's inbuilt rules (which I seem
# unable to remove, even when using -d) and it has some quirky ideas as to
# which files are created temporarily (and thus should be removed). It's
# best to usually do 'gmake spotless all' or some such.
#
# The input files are always .texi
#
# .texinfo files are expanded up .texi files. They have the macros replaced
# and have been processed by m4 to include Unix or Windows specific components.
#
# The _us.* files are US-letter format copies, otherwise we assume everything
# is in A4 format.
#

#-----------------------------------------------------------------------------
# General rules

#
# Unix vs Windows rules.
# The Unix and Windows manuals stem from the same text, but using m4 as a
# preprocessor to generate different .texinfo files.
#
%_unix.texinfo:	%.texi
	m4 -Uunix -D_unix -D_tex -D_html < $< > $@

%_windows.texinfo:	%.texi
	m4 -Uunix -D_windows -D_tex -D_html < $< > $@

# How to build .dvi files from our macro-expanded .texinfo
# files
%.dvi:	%.texinfo
	texi2dvi $<

# US-letter PostScript from DVI
%_us.ps:	%_us.dvi
	dvips -t letter -o $@ $<

# A4 PostScript from DVI
#%.ps:	%.dvi
#	dvips -t a4 -o $@ $<

%.ps:	%.dvi
	dvips -t a4 -Ppdf -o $@ $<

%.pdf:	%.ps
#	PATH=$$PATH:/usr/local/bin ps2pdf -dPDFSETTINGS=/prepress -dEncodeColorImages=false $<
	PATH=$$PATH:/usr/local/bin ps2pdf $<

# Converts an A4 formatted .texi file into a US-letter formatted file.
%_us.texinfo: %.texinfo
	egrep -v '^@afourpaper' < $< > $@

#-----------------------------------------------------------------------------
# The main make targets.

all: unix windows

unix: course_unix.ps course_unix.pdf
	cp course_unix.ps course_unix.pdf ../unix_docs

windows: course_windows.ps course_windows.pdf
	cp course_windows.ps course_windows.pdf ../windows_docs

clean:
	-rm *.aux *.cp *.fn *.ky *.log *.pg *.toc *.tp *.vr *.cps
	-rm core _tmp.texi _tmp.texi~ *.texinfo *.texinfo.tmp

spotless:	clean
	-rm *.dvi *_toc.html *_[0-9]*.html *.info *.info-[0-9] *.index *.topic
	-rm master_index.html master_contents.html master_brief.html
	-rm *.pdf *.ps
