#all: unix windows unix_us windows_us gzip
all: unix windows

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
# Unix vs Windows rules.
# The Unix and Windows manuals stem from the same text, but using m4 as a
# preprocessor to generate different .texinfo files. We also generate different
# documents designed for passing into TeX or texi2html. The changes are
# minimal and relate to cross-referencing and page-splitting.
#
%_unix.texinfo:	%.texi header.m4
	$(M4) -Uunix -D_unix -D_tex < $< > $@
	./tools/update-nodes $@

%_unix.htmlinfo:	%.texi header.m4
	$(M4) -Uunix -D_unix -D_html < $< > $@
	./tools/update-nodes $@

%_windows.texinfo:	%.texi header.m4
	$(M4) -Uunix -D_windows -D_tex < $< > $@
	./tools/update-nodes $@

%_windows.htmlinfo:	%.texi header.m4
	$(M4) -Uunix -D_windows -D_html < $< > $@
	./tools/update-nodes $@

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

# PDF generation. Directly from the texinfo.
# Note that as the texinfo is reordered, this produces out of order
# data too (contents page at the end).
# One solution is (eg):
#   pdftk manual_unix.pdf cat 1-2 469-end 3-468 output new.pdf
# however this loses the bookmarks.
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
	./tools/texi2html -menu -verbose -split_chapter -index_chars $<

# Man pages - taken from the .texi files directly (NB: may need m4 expansion,
# but this is not yet applied.)
man/man1/%.1:	%.1.texi
	./tools/texi2man.pl $< > $@

man/man4/ExperimentFile.4: exp-t.texi
	./tools/texi2man.pl $< > $@

man/man4/scf.4: scf-t.texi
	./tools/texi2man.pl $< > $@

man/man4/ztr.4: ztr-t.texi
	./tools/texi2man.pl $< > $@

# For any large pictures (_lpicure() macro) we create an html page containing
# the full-size picture to link from the small picture embedded in the page.
gifs:
	./tools/make_gif_html

pngs:
	./tools/make_png_html

# PostScript versions of the gif images.
%.ps: %.gif
	./tools/make_ps $<

%.pdf: %.png
	./tools/make_pdf $<


# The mini manuals need to have the xrefs removed. Internal references are
# kept "as is".
mini_unix.texinfo: mini_manual_unix.texinfo
	./tools/remove_xrefs.pl < $< > $@
mini_unix.htmlinfo: mini_manual_unix.texinfo
	./tools/remove_xrefs.pl < $< > $@

mini_windows.texinfo: mini_manual_windows.texinfo
	./tools/remove_xrefs.pl < $< > $@
mini_windows.htmlinfo: mini_manual_windows.texinfo
	./tools/remove_xrefs.pl < $< > $@

# Help within the programs is HTML based. Given a topic the appropriate URL
# is obtained by looking it up in the .index file.
%.index:	%_toc.html
	./tools/html_index.pl $<

# Resolves cross references between separate html texinfo documents by
# searching for <!-- XREF:name --> comments in the html code and matching
# these up with the node names.
# We need to have a different index for unix and windows versions as some
# pages are missing in unix which upsets the indexing.
xref_unix: \
	formats_unix.index \
	pregap4_unix.index \
	read_clipping_unix.index \
	vector_clip_unix.index \
	interface_unix.index \
	trev_unix.index \
	gap4_unix.index \
	manpages_unix.index \
	spin_unix.index \
	filebrowser_unix.index 

	./tools/xref_update.pl *_unix_*.html

xref_windows: \
	formats_windows.index \
	pregap4_windows.index \
	read_clipping_windows.index \
	vector_clip_windows.index \
	interface_windows.index \
	trev_windows.index \
	gap4_windows.index \
	manpages_windows.index \
	spin_windows.index \
	filebrowser_windows.index 

	./tools/xref_update.pl *_windows_*.html

# Master full and brief contents page generation. Basically this combines
# the separate _toc.html pages into global ones.
contents_unix:
	./tools/merge_indexes.pl unix \
		gap4_unix_toc.html \
		mutations_unix_toc.html \
		pregap4_unix_toc.html \
		read_clipping_unix_toc.html \
		vector_clip_unix_toc.html \
		trev_unix_toc.html \
		spin_unix_toc.html \
		interface_unix_toc.html \
		formats_unix_toc.html \
		manpages_unix_toc.html \
		references_unix_toc.html

contents_windows:
	./tools/merge_indexes.pl windows \
		gap4_windows_toc.html \
		mutations_windows_toc.html \
		pregap4_windows_toc.html \
		read_clipping_windows_toc.html \
		vector_clip_windows_toc.html \
		trev_windows_toc.html \
		spin_windows_toc.html \
		interface_windows_toc.html \
		formats_windows_toc.html \
		manpages_windows_toc.html \
		references_windows_toc.html

# Backup only the 'source' files and generating scripts
backup::
	tar cvf - Makefile README README.system *.template *.tcl tclIndex \
	          docmake list_nodes lowersection html_index update-nodes \
		  update-nodes.el *.gif *.texi *.pl texinfo.tex \
	  | gzip > backup/`date +"%d_%m_%y"`.tar.gz

#-----------------------------------------------------------------------------
# The main make targets.

unix: manual_unix xref_unix contents_unix mini_unix man pngs
windows: manual_windows xref_windows contents_windows mini_windows man pngs

unix_us: manual_unix_us.pdf mini_unix_us.pdf
windows_us: manual_windows_us.pdf mini_windows_us.pdf

mini_unix: mini_unix.pdf mini_unix_toc.html
	./tools/edit_mini_contents.pl unix < mini_unix_toc.html > tmp.html
	mv tmp.html mini_unix_toc.html

mini_windows: mini_windows.pdf mini_windows_toc.html
	./tools/edit_mini_contents.pl windows < mini_windows_toc.html > tmp.html
	mv tmp.html mini_windows_toc.html

manual_unix:	 manual_unix.pdf manual_unix_html
manual_unix_html: \
	interface_unix_toc.html \
	gap4_unix_toc.html \
	formats_unix_toc.html \
	vector_clip_unix_toc.html \
	trev_unix_toc.html \
	manpages_unix_toc.html \
	read_clipping_unix_toc.html \
	references_unix_toc.html \
	spin_unix_toc.html \
	pregap4_unix_toc.html \
	mutations_unix_toc.html

manual_windows:	manual_windows.pdf manual_windows_html
manual_windows_html: \
	interface_windows_toc.html \
	gap4_windows_toc.html \
	formats_windows_toc.html \
	vector_clip_windows_toc.html \
	trev_windows_toc.html \
	manpages_windows_toc.html \
	read_clipping_windows_toc.html \
	references_windows_toc.html \
	spin_windows_toc.html \
	pregap4_windows_toc.html \
	mutations_windows_toc.html

man:	man/man1/makeSCF.1 man/man1/eba.1 man/man1/makeSCF.1 \
	man/man1/get_scf_field.1 man/man1/init_exp.1 man/man4/scf.4 \
	man/man4/ExperimentFile.4 \
	man/man1/extract_seq.1 man/man1/copy_db.1 \
	man/man4/ztr.4 man/man1/convert_trace.1 \
	man/man1/find_renz.1 man/man1/get_comment.1 \
	man/man1/getABIfield.1 man/man1/trace_dump.1

# For distributions
gzip:
	-gzip < mini_windows.ps > mini_windows.ps.gz
	-gzip < mini_windows_us.ps > mini_windows_us.ps.gz
	-gzip < mini_unix.ps > mini_unix.ps.gz
	-gzip < mini_unix_us.ps > mini_unix_us.ps.gz
	-gzip < manual_windows.ps > manual_windows.ps.gz
	-gzip < manual_windows_us.ps > manual_windows_us.ps.gz
	-gzip < manual_unix.ps > manual_unix.ps.gz
	-gzip < manual_unix_us.ps > manual_unix_us.ps.gz

clean:
	-rm *.aux *.cp *.fn *.ky *.log *.pg *.toc *.tp *.vr *.cps *.fns *.pgs *.vrs
	-rm core _tmp.texi _tmp.texi~ *.texinfo *.texinfo.tmp *.texinfo~
	-rm *.htmlinfo *.htmlinfo~

spotless:	clean
	-rm *.dvi *.html *.info *.info-[0-9] *.index *.topic
	-rm master_index.html master_contents.html master_brief.html
	-rm manual_unix.ps manual_windows.ps mini_unix.ps mini_windows.ps
	-rm man/man1/* man/man4/*
	-rm *.ps *.pdf

depend:
	./tools/make_dependencies > dependencies

include dependencies
