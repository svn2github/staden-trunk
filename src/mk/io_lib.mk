SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

FLAGS	 =
DEFS	 = INCLUDES_S="-I../include" SRCROOT=../..

#
# Define the location where you wish libraries, include files and man
# pages to be placed.
#
INST_PREFIX=/usr/local
MAN_DIR	= $(INST_PREFIX)/man
LIB_DIR = $(INST_PREFIX)/lib
HEA_DIR = $(INST_PREFIX)/include

# For compilation external to the rest of the staden package the following
# need uncommenting and editing to suit the local system. (Set RANLIB to
# /bin/true on systems not requiring it.)
#
#DEFS	 = CC=cc\
#	   CFLAGS="-g -I../include"\
#	   AR=ar\
#	   ARFLAGS=rv\
#	   RANLIB=ranlib

# Libraries to compile.
# Define utils at the start, and read as the last.
#
DIRS	= utils scf exp_file plain abi alf read progs

JOB = all 
 
all: $(DIRS) 
 
$(DIRS): FORCE 
	cd $@ && $(MAKE) $(MFLAGS) $(FLAGS) $(DEFS) $(JOB) 
 
FORCE: 

SUBTARGETS = depend clean relink spotless install 
$(SUBTARGETS): FORCE 
	$(MAKE) $(FLAGS) JOB=$@ 

ext_install:
	if [ ! -d "$(INST_PREFIX)" ]; then mkdir -p $(INST_PREFIX); fi
	if [ ! -d "$(LIB_DIR)" ]; then mkdir $(LIB_DIR); fi
	-cp lib/$(O)/*.a $(LIB_DIR)
	if [ ! -d "$(HEA_DIR)" ]; then mkdir $(HEA_DIR); fi
	-cp include/*.h $(HEA_DIR)
	if [ ! -d "$(MAN_DIR)" ]; then mkdir $(MAN_DIR); fi
	-cp -r man/* $(MAN_DIR)

install:
	for i in $(DIRS);\
	do (cd $$i;\
	    $(MAKE) $(FLAGS) $(DEFS) install);\
	done
	cd progs; $(MAKE) $(FLAGS) $(DEFS) install

.PHONY: distsrc
distsrc:
	-cp -R Makefile COPYRIGHT README man include $(DIRNAME)
	-for i in abi alf exp_file plain progs read scf utils; \
	do (cd $$i; $(MAKE) $(FLAGS) $(DEFS) $@); done
	-rmdir $(DIRNAME)/*-binaries

dist:
	rm -rf /tmp/io_lib
	mkdir /tmp/io_lib
	mkdir /tmp/io_lib/scf /tmp/io_lib/plain /tmp/io_lib/exp_file \
		/tmp/io_lib/utils /tmp/io_lib/read /tmp/io_lib/progs \
		/tmp/io_lib/abi /tmp/io_lib/alf
	for i in /tmp/io_lib/*; \
	    do \
		mkdir $$i/alpha-binaries; \
		mkdir $$i/sun-binaries; \
		mkdir $$i/solaris-binaries; \
		mkdir $$i/sgi-binaries; \
		mkdir $$i/linux-binaries; \
		mkdir $$i/WINNT-binaries; \
	    done
	mkdir /tmp/io_lib/mk /tmp/io_lib/data
	cp README COPYRIGHT dependencies /tmp/io_lib
	cp Makefile.standalone /tmp/io_lib/Makefile
	cp -R man include /tmp/io_lib/
	rm /tmp/io_lib/include/os.h /tmp/io_lib/include/misc.h
	cp include/os.h include/misc.h /tmp/io_lib/include
	cp data/* /tmp/io_lib/data
	cp -R ../mk/alpha.mk ../mk/global.mk ../mk/solaris.mk \
		../mk/sun.mk ../mk/sgi.mk ../mk/linux.mk \
		../mk/windows.mk ../mk/windows /tmp/io_lib/mk
	cp -R scf/*.[ch] scf/README scf/scf1 scf/dependencies \
		scf/Makefile scf/Makefile.static /tmp/io_lib/scf
	cp plain/*.[ch] plain/Makefile plain/dependencies \
		plain/Makefile.static /tmp/io_lib/plain
	cp exp_file/*.[ch] exp_file/Makefile exp_file/dependencies \
		exp_file/Makefile.static /tmp/io_lib/exp_file
	cp utils/*.[ch] utils/README utils/dependencies utils/Makefile \
		utils/Makefile.static /tmp/io_lib/utils
	cp -R read/*.[ch] read/README read/dependencies \
		read/Makefile.static read/scf1 /tmp/io_lib/read
	cp read/Makefile.standalone /tmp/io_lib/read/Makefile
	cp -R progs/*.[ch] progs/Makefile progs/dependencies \
		/tmp/io_lib/progs
	cp abi/*.[ch] abi/Makefile abi/dependencies \
		abi/Makefile.static /tmp/io_lib/abi
	cp alf/*.[ch] alf/Makefile alf/dependencies \
		alf/Makefile.static /tmp/io_lib/alf
	-rm /tmp/io_lib/*/*.old.c /tmp/io_lib/*/*.old.h /tmp/io_lib/*/*~
	mkdir /tmp/io_lib/lib
	mkdir /tmp/io_lib/lib/alpha-binaries
	mkdir /tmp/io_lib/lib/solaris-binaries
	mkdir /tmp/io_lib/lib/sun-binaries
	mkdir /tmp/io_lib/lib/sgi-binaries
	mkdir /tmp/io_lib/lib/linux-binaries
	mkdir /tmp/io_lib/lib/WINNT-binaries

	@echo
	@echo
	@echo '******'
	@echo 'Don't forget to update the version number in the README file'
	@echo '******'
	@echo
