.PHONY: dist install

dist:
	-mkdir $(DIST)/demo
	-cp -R gap4 $(DIST)/demo

distsrc:
	$(MAKE) $(SUBFLAGS) DIST=$(DISTSRC) dist
	cp Makefile $(DISTSRC)/demo

install:
