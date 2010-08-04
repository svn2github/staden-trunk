all: 
	@echo
	@echo Please rerun make specifying either target \"unix\" or \"windows\".
	@echo

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
%.texinfo:	%.texi header.m4
	$(M4) $(M4OPT) -D_tex < $< > $@
	./tools/update-nodes $@

%.htmlinfo:	%.texi header.m4
	$(M4) $(M4OPT) -D_html < $< > $@
	./tools/update-nodes $@

# Remove implcit rules
%.dvi:  %.texi
%.dvi:  %.texinfo
%.pdf:	%.texinfo

%.eps_done:

# How to build .dvi files from our m4-expanded .texinfo files
%.dvi:	%.texinfo %.eps_done
	texi2dvi $<

# A4 or US Letter PostScript from DVI
%.ps:	%.dvi
	dvips -t $(PAPER) -Ppdf -o $@ $<

# PDF generation. Directly from the texinfo.
# Note that as the texinfo is reordered, this produces out of order
# data too (contents page at the end).
# One solution is (eg):
#   pdftk manual.pdf cat 1-2 469-end 3-468 output new.pdf
# however this loses the bookmarks.
%.pdf:	%.texinfo %.pdf_done
	texi2pdf $<

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

%.eps: %.png
	./tools/make_eps $<

%.pdf: %.png
	./tools/make_pdf $<


# The mini manuals need to have the xrefs removed. Internal references are
# kept "as is".
mini.texinfo: mini_manual.texinfo
	./tools/remove_xrefs.pl < $< > $@
mini.htmlinfo: mini_manual.texinfo
	./tools/remove_xrefs.pl < $< > $@

# Help within the programs is HTML based. Given a topic the appropriate URL
# is obtained by looking it up in the .index file.
%.index:	%_toc.html
	./tools/html_index.pl $<

# Resolves cross references between separate html texinfo documents by
# searching for <!-- XREF:name --> comments in the html code and matching
# these up with the node names.

xref: \
	formats.index \
	pregap4.index \
	read_clipping.index \
	vector_clip.index \
	interface.index \
	trev.index \
	gap4.index \
	manpages.index \
	spin.index \
	filebrowser.index 

	./tools/xref_update.pl *_*.html

# Master full and brief contents page generation. Basically this combines
# the separate _toc.html pages into global ones.
contents:
	./tools/merge_indexes.pl \
		gap4_toc.html \
		mutations_toc.html \
		pregap4_toc.html \
		read_clipping_toc.html \
		vector_clip_toc.html \
		trev_toc.html \
		spin_toc.html \
		interface_toc.html \
		formats_toc.html \
		manpages_toc.html \
		references_toc.html

# Backup only the 'source' files and generating scripts
backup::
	tar cvf - Makefile README README.system *.template *.tcl tclIndex \
	          docmake list_nodes lowersection html_index update-nodes \
		  update-nodes.el *.gif *.texi *.pl texinfo.tex \
	  | gzip > backup/`date +"%d_%m_%y"`.tar.gz

#-----------------------------------------------------------------------------
# The main make targets.

unix: M4OPT+=-Uunix -D_unix
unix: common

windows: M4OPT+=-Uunix -D_windows
windows: common

ifeq ($(PAPER),us)
M4OPT += -Dafourpaper=c
PAPER=letter
else
M4OPT += -Dafourpaper=c
PAPER=m4
endif

common: manual xref contents mini man pngs


mini.pdf_done:

mini: mini.pdf mini_toc.html
	./tools/edit_mini_contents.pl < mini_toc.html > tmp.html
	mv tmp.html mini_toc.html

manual:	 manual.pdf manual_html
manual_html: \
	interface_toc.html \
	gap4_toc.html \
	formats_toc.html \
	vector_clip_toc.html \
	trev_toc.html \
	manpages_toc.html \
	read_clipping_toc.html \
	references_toc.html \
	spin_toc.html \
	pregap4_toc.html \
	mutations_toc.html

man:	man/man1/makeSCF.1 man/man1/eba.1 man/man1/makeSCF.1 \
	man/man1/get_scf_field.1 man/man1/init_exp.1 man/man4/scf.4 \
	man/man4/ExperimentFile.4 \
	man/man1/extract_seq.1 man/man1/copy_db.1 \
	man/man4/ztr.4 man/man1/convert_trace.1 \
	man/man1/find_renz.1 man/man1/get_comment.1 \
	man/man1/getABIfield.1 man/man1/trace_dump.1


clean:
	-rm *.aux *.cp *.fn *.ky *.log *.pg *.toc *.tp *.vr *.cps *.fns *.pgs *.vrs
	-rm core _tmp.texi _tmp.texi~ *.texinfo *.texinfo.tmp *.texinfo~
	-rm *.htmlinfo *.htmlinfo~

spotless:	clean
	-rm *.dvi *.html *.info *.info-[0-9] *.index *.topic
	-rm master_index.html master_contents.html master_brief.html
	-rm manual.ps mini.ps
	-rm man/man1/* man/man4/*
	-rm *.ps *.pdf *.eps

depend:
	./tools/make_dependencies > dependencies

include dependencies
