.PHONY: dist install

dist:
	-mkdir $(DIST)/tables
	-cp -R alus $(DIST)/tables
	-cp -R blast_data $(DIST)/tables
	-cp -R gcodes $(DIST)/tables
	-cp -R vectors $(DIST)/tables
	-cp *.wts $(DIST)/tables
	-cp *.bitmap $(DIST)/tables
	-cp README $(DIST)/tables
	-cp GTAGDB $(DIST)/tables
	-cp NOTEDB $(DIST)/tables
	-cp FEAT*DB $(DIST)/tables
	-cp RENZYM.4 $(DIST)/tables
	-cp RENZYM.6 $(DIST)/tables
	-cp RENZYM.ALL $(DIST)/tables
	-cp renzyme_bairoch $(DIST)/tables
	-cp cap2rc $(DIST)/tables
	-cp cap3rc $(DIST)/tables
	-cp fak2rc $(DIST)/tables
	-cp gaprc $(DIST)/tables
	-cp gaprc_menu_demo $(DIST)/tables
	-cp gaprc_menu_full $(DIST)/tables
	-cp gaprc_menu_mito $(DIST)/tables
	-cp gaprc_menu_viewer $(DIST)/tables
	-cp niprc $(DIST)/tables
	-cp phraprc $(DIST)/tables
	-cp pregap4rc $(DIST)/tables
	-cp siprc $(DIST)/tables
	-cp tk_utilsrc $(DIST)/tables
	-cp spinrc $(DIST)/tables
	-cp spin2rc $(DIST)/tables
	-cp seqedrc $(DIST)/tables
	-cp itcl3.3rc $(DIST)/tables
	-cp itk3.3rc $(DIST)/tables
	-cp iwidgetsrc $(DIST)/tables
	-cp spin_embossrc $(DIST)/tables
	-cp spin2_embossrc $(DIST)/tables
	-cp emboss_menu $(DIST)/tables
	-cp trevrc $(DIST)/tables
	-cp shlib.conf $(DIST)/tables
	-cp help_config $(DIST)/tables
	-cp pam250 $(DIST)/tables
	-cp m13mp18_primers $(DIST)/tables
	-cp vector_primer $(DIST)/tables
	-cp nuc_matrix $(DIST)/tables
	-cp align_lib_nuc_matrix $(DIST)/tables
	-cp copy_readsrc $(DIST)/tables
	-cp prefinishrc $(DIST)/tables
	-cp config.prefinish $(DIST)/tables
	-cp *.bmp $(DIST)/tables

distsrc:
	$(MAKE) $(SUBFLAGS) DIST=$(DISTSRC) dist
	cp Makefile $(DISTSRC)/tables

install:


