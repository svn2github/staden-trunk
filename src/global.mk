# Global makefile.
# We must have $(PROGS) defined before this is included otherwise the target
# dependancies for all are not set up correctly.

# 7/1/99 johnt - defining this causes problems with windows, should be in $(MACHINE).mk
# Just incase (as in the IRIX situation)...
#SHELL	= /bin/sh

# 26/3/99 johnt - added Corba support when CORBA defined
#         hopefully all corba dependant stuff is in this file
#         although other makefiles contain corba references
#         everything is defined here

#------------------------------------------------------------------------------
# Installation directories - used by make install.
#
# INSTALLBIN		Location to install binary executable files
#
# INSTALLSCRIPT		Location to install scripts. These can be executed
#			by any system, so you may wish to have
#			scripts in .../bin and binaries in .../bin/`arch`.
#
# INSTALLLIB		Location to install dynamic libraries.
#
# INSTALLSEQBIN		As for INSTALLBIN, except holding the sequence library
#			indexing software
#
# INSTALLSEQSCRIPT	As for INSTALLSCRIPT, except holding the sequence
#			library indexing software
#
INSTALLBIN	= $(STADENROOT)/$(MACHINE)-bin
INSTALLSCRIPT	= $(STADENROOT)/$(MACHINE)-bin
INSTALLLIB	= $(STADENROOT)/lib
INSTALLSEQBIN	= $(STADENROOT)/seqlibs-$(MACHINE)-bin
INSTALLSEQSCRIPT= $(STADENROOT)/seqlibs-$(MACHINE)-bin

#------------------------------------------------------------------------------
# Options listed here ought be available to all architectures. If a system 
# does not have a particular option then leave it out here and add it to all
# the relevant system makefiles. In cases where one system is the exception
# it may be best to simply treat these as defaults and then override the
# parameters in the system makefile.

# Debugging and optimising switches; define [CF]OPTDEBUG to be one or the other
CDEBUG		= -g
FDEBUG		= -g -C
COPT		= -O2 -g3
FOPT		= -O2 -g3
COPTDEBUG	= $(CDEBUG)
FOPTDEBUG	= $(FDEBUG)
#COPTDEBUG	= $(COPT)
#FOPTDEBUG	= $(FOPT)
# 7/1/99 johnt - also need to modify SHLIB_OPTDEBUG and CLDOPTDEBUG in windows.mk

#SHLIB_STRIP	= -x

# C compiler defines
CC		= cc
CXX             = g++
#DEFINES	= -DUSE_SRS
DEFINES        += -DUSE_NON_CONST
NO_SRS=1

# Uncomment this to enable biolims support. Also see io_lib/options.mk.
#BIOLIMS=1
ifdef BIOLIMS
DEFINES		+= -DUSE_BIOLIMS
endif
CFLAGS		= $(COPTDEBUG) $(DEFINES) $(INCLUDES)
CXXFLAGS	= $(COPTDEBUG) $(DEFINES) $(INCLUDES)


# 7/1/99 johnt - Command to make Windows Def file for each object file - defaults to echo - required by Visual C++
MKDEFC = @\#

# 7/1/99 johnt - flag to specify output object file - defaults to "-o " - required by Visual C++
SPACE=
COBJFLAG = -o $(SPACE)
FOBJFLAG = -o $(SPACE)
LDEXEFLAG = -o $(SPACE)


# FORTRAN compiler defines
F77		= f77
FFLAGS		= $(FOPTDEBUG)

# Linker definitions
# 7/1/99 johnt - added Linker abstractions required for Visual C++ suppoort
CLDOPTDEBUG = $(COPTDEBUG)
LINK_LIBFLAG = -l
LINK_PATHFLAG = -L
LIB_EXT      =

CLD_PROG	= $(CC)
CXXLD_PROG	= $(CXX)
CLDFLAGS	= $(CLDFLAGS_S) $(CLDOPTDEBUG) $(LINK_PATHFLAG)$(L) $(CLDFLAGS_E)
CXXLDFLAGS      = $(CLDFLAGS_S) $(CLDOPTDEBUG) $(LINK_PATHFLAG)$(L) $(CLDFLAGS_E)
CLD		= $(CLD_PROG) $(CLDFLAGS)
CXXLD           = $(CXXLD_PROG) $(CXXLDFLAGS)
FLD_PROG	= $(F77)
FLDFLAGS	= $(FLDFLAGS_S) $(FOPTDEBUG) $(LINK_PATHFLAG)$(L) $(FLDFLAGS_E)
FLD		= $(FLD_PROG) $(FLDFLAGS)
GAP4SH_LD       = $(FLD) $(FLDFLAGS)
LINKF		= $(FLD) $(FLDFLAGS) $(FOBJFLAG) args $(LIBSF)
# At the end of every compile line. Useful in the system makefiles.
LIBSF		= $(LIBSF_S) $(EXTRA_LIBS) $(LIBSF_E)
LIBSC		= $(LIBSC_S) $(EXTRA_LIBS) $(LIBSC_E)
EXTRA_LIBS	=

RANLIB		= ranlib


# Default includes
VPATH           = $(SRC)
INCLUDES	= $(INCLUDES_S) -I$(SRC) $(INCLUDES_E) -I$(BUILD)

STADLIB		= ../../lib
TCLBIN		= $(L)
TKBIN		= $(L)
ITCLBIN		= $(L)
ITKBIN		= $(L)
TTBIN		=
TCLVERS		= 8.4
TKVERS		= 8.4
TCLSRC		= $(SRCROOT)/tcl8.4.6/generic
TKSRC		= $(SRCROOT)/tk8.4.6/generic
ITCLSRC		= $(SRCROOT)/incrTcl-3.3cvs/itcl/generic
ITKSRC		= $(SRCROOT)/incrTcl-3.3cvs/itk/generic
TTSRC		=
ITCLVERS	= 3.3
ITKVERS		= 3.3

#26/3/99 johnt - added corba support
ifdef CORBA
CORBA_OBJS      = corba.o trace.o basicServer.o

CORBA_INCDIR    = /usr/local/mico/include
CORBA_LIBDIR    = /usr/local/mico/lib
COSS_VER        = micocoss225
CORBA_VER       = mico225
endif

# 7/1/99 johnt - added abstractions to support Visual C++
# Standard chunks to add to the link line
ifdef CORBA
CORBA_LIB    = $(LINK_PATHFLAG)$(CORBA_LIBDIR) $(LINK_LIBFLAG)$(COSS_VER)$(LIB_EXT) \
		$(LINK_LIBFLAG)$(CORBA_VER)$(LIB_EXT)
endif
MATH_LIB     = $(MATH_LIB_S) -lm $(MATH_LIB_E)
MISC_LIB     = $(MISC_LIB_S) $(LINK_LIBFLAG)misc$(LIB_EXT) $(MISC_LIB_E)
TCL_LIB	     = $(TCL_LIB_S) $(LINK_PATHFLAG)$(TCLBIN) $(LINK_LIBFLAG)tcl$(TCLVERS)$(LIB_EXT) $(MATH_LIB) \
	           $(TCL_LIB_E)
TK_LIB	     = $(TK_LIB_S) $(LINK_PATHFLAG)$(TKBIN) $(LINK_LIBFLAG)tk$(TKVERS)$(LIB_EXT) $(TCL_LIB) $(TK_LIB_E)
ITCL_LIB	     = $(ITCL_LIB_S) $(LINK_PATHFLAG)$(ITCLBIN) $(LINK_LIBFLAG)itcl$(ITCLVERS)$(LIB_EXT) $(MATH_LIB) \
	           $(ITCL_LIB_E)
ITK_LIB	     = $(ITK_LIB_S) $(LINK_PATHFLAG)$(ITKBIN) $(LINK_LIBFLAG)itk$(ITKVERS)$(LIB_EXT) $(ITCL_LIB) $(ITK_LIB_E)
TT_LIB	     = $(TT_LIB_S) $(TTBIN:%=-L%) $(TT_LIBRARY) $(TT_LIB_E)
# io-utils and read libraries have now been merged into one.
SCF_LIB	     = $(SCF_LIB_S) $(LINK_LIBFLAG)scf$(LIB_EXT) $(LINK_LIBFLAG)io-utils$(LIB_EXT) $(SCF_LIB_E)
EXP_LIB	     = $(EXP_LIB_S) $(LINK_LIBFLAG)exp$(LIB_EXT) $(LINK_LIBFLAG)io-utils$(LIB_EXT) $(EXP_LIB_E)
G_LIB	     = $(G_LIB_S) $(LINK_LIBFLAG)g$(LIB_EXT) $(G_LIB_E)
GAP_LIB	     = $(GAP_LIB_S) $(LINK_LIBFLAG)gap$(LIB_EXT) $(GAP_LIB_E)
TCLUTILS_LIB =
TKUTILS_LIB  = $(TKUTILS_LIB_S) $(LINK_LIBFLAG)tk_utils$(LIB_EXT) $(TKUTILS_LIB_E)
TEXTUTILS_LIB= $(TEXTUTILS_LIB_S) $(LINK_LIBFLAG)text_utils$(LIB_EXT) $(TK_LIB) $(TEXTUTILS_LIB_E)
SEQUTILS_LIB = $(SEQUTILS_LIB_S) $(LINK_LIBFLAG)seq_utils$(LIB_EXT) $(SEQUTILS_LIB_E)
SEQLIB_LIB   = $(SEQLIB_LIB_S) $(LINK_LIBFLAG)seqlib$(LIB_EXT) $(SEQ_LIB_E)
SPIN_LIB     = $(SPIN_LIB_S) $(LINK_LIBFLAG)spin$(LIB_EXT) $(SPIN_LIB_E)
COPYREADS_LIB = $(COPYREADS_LIB_S) $(LINK_LIBFLAG)copy_reads$(LIB_EXT) $(COPYREADS_LIB_E)
SCFEXPIO_LIB = $(SCFEXPIO_LIB_S) $(LINK_PATHFLAG)$(SRCROOT)/fakii/scf_exp_io/$(O) $(LINK_LIBFLAG)scf_exp_io$(LIB_EXT) $(SCFEXPIO_LIB_E)
ZLIB_LIB     = $(ZLIB_LIB_S) $(LINK_LIBFLAG)z$(LIB_EXT) $(ZLIB_LIB_E)
MUT_LIB      = $(MUT_LIB_S) $(LINK_LIBFLAG)mutlib$(LIB_EXT) $(MUT_LIB_E)
P3_LIB	     = $(P3_LIB_S)  $(LINK_LIBFLAG)primer3$(LIB_EXT) $(P3_LIB_E)
PNG_LIB	     = $(PNG_LIB_S) $(LINK_LIBFLAG)png12$(LIB_EXT) $(PNG_LIB_E)


# Standard chunks to add to the compile line
ifdef CORBA
CORBA_INC       = -I$(CORBA_INCDIR)
endif
MISC_INC	= -I$(MISCSRC)
TCL_INC		= -I$(TCLSRC)
TK_INC		= -I$(TKSRC) -I$(TCLSRC)
ITCL_INC	= -I$(ITCLSRC)
ITK_INC		= -I$(ITKSRC) -I$(ITCLSRC)
TT_INC		= $(TTSRC:%=-I%)
G_INC		= -I$(GSRC)
GAP4_INC	= -I$(GAP4SRC)
NIP4_INC	= -I$(NIP4SRC)
TCLUTILS_INC	=
TKUTILS_INC	= -I$(TKUTILSSRC)
TEXTUTILS_INC	= -I$(TEXTUTILSSRC)
SEQUTILS_INC	= -I$(SEQUTILSSRC)
SEQLIB_INC	= -I$(SEQLIBSRC)
SPIN_INC	= -I$(SPINSRC)
COPYREADS_INC	= -I$(COPYREADSSRC)
MUT_INC		= -I$(SRCROOT) -I$(MUTSRC)
P3_INC		= -I$(SRCROOT)/primer3/src
PNG_INC	        = -I$(SRCROOT)/libpng

#
# Where the objects are relative to this makefile and vice versa.
# Use this to shorten later references (makes things tidier).
#
# Src to Obj
#O = $(MACHINE)-binaries
O=.
# Obj to Src
S = $(SRCROOT)/$(SUBDIR)
# Lib
L = $(SRCROOT)/lib
LD_LIBRARY_PATH := $(L):$(LD_LIBRARY_PATH)

#
# Default Tool Talk object files, for gap only.
#
TT_OBJS=


#------------------------------------------------------------------------------
# MJ: Appends .exe extension for windows executables, $(EXE_SUFFIX) expands
#     to an empty string for UNIX.
#------------------------------------------------------------------------------
PROGSTOCLEAN=$(PROGS:=$(EXE_SUFFIX))



#------------------------------------------------------------------------------
# Default targets
#
all:	$(PROGS)

CLEANCOMMAND=-rm -f *.a *.o
clean:
	$(CLEANCOMMAND)

cleanprogs:
	-rm -f $(PROGSTOCLEAN)


CLEANLIBSCOMMAND=-rm -f $(PROGLIBS)

cleanlibs:
	$(CLEANLIBSCOMMAND)

relink: cleanprogs cleanlibs all

spotless: clean cleanprogs cleanlibs


#------------------------------------------------------------------------------
# Swap these lines around if you have a fortran compiler and wish to
# edit legacy.f. To keep dependencies down and simplify dynamic linking for
# now we just use the (edited) f2c derived version instead.
#GAP4_LEGACY	= legacy.o
GAP4_LEGACY	= legacy_f2c.o

#------------------------------------------------------------------------------
# Where to put things
#
ABISRC		= $(SRCROOT)/abi
ALFSRC		= $(SRCROOT)/alf
CLIPSRC		= $(SRCROOT)/clip
CONVERTSRC	= $(SRCROOT)/convert
COPSRC		= $(SRCROOT)/cop
EBASRC		= $(SRCROOT)/eba
EXPGETSEQSRC	= $(SRCROOT)/expGetSeq
FROGSRC		= $(SRCROOT)/frog
GSRC		= $(SRCROOT)/g
GAPSRC		= $(SRCROOT)/gap
GAP4SRC		= $(SRCROOT)/gap4
REPE_OBJ	= repe.o
NIP4SRC		= $(SRCROOT)/nip4
GETMCHSRC	= $(SRCROOT)/getMCH
GETSCFFIELD	= $(SRCROOT)/get_scf_field
INDEXSEQLIBSSRC	= $(SRCROOT)/indexseqlibs
INITEXPSRC	= $(SRCROOT)/init_exp
LOOKUPSRC	= $(SRCROOT)/lookup
MAKESCFSRC	= $(SRCROOT)/makeSCF
MISCSRC		= $(SRCROOT)/Misc
OSPSRC	  	= $(GAPSRC)/osp-bits
OSP4SRC	  	= $(GAP4SRC)/osp-bits
SEQUTILSSRC	= $(SRCROOT)/seq_utils
SEQLIBSRC	= $(SRCROOT)/seqlib
SIP4SRC		= $(SRCROOT)/sip4
SPINSRC		= $(SRCROOT)/spin
COPYREADSSRC	= $(SRCROOT)/copy_reads
STADENSRC	= $(SRCROOT)/staden
TCLUTILSSRC	=
TKUTILSSRC	= $(SRCROOT)/tk_utils
TEXTUTILSSRC	= $(SRCROOT)/text_utils
TEDSRC		= $(SRCROOT)/ted
TRACECLIPSRC	= $(SRCROOT)/trace_clip
TREVSRC		= $(SRCROOT)/trev
VEPESRC		= $(SRCROOT)/vepe
QCLIPSRC	= $(SRCROOT)/qclip
MUTSRC		= $(SRCROOT)/mutlib

ABIBIN		= $(ABISRC)/$(O)
ALFBIN		= $(ALFSRC)/$(O)
CLIPBIN		= $(CLIPSRC)/$(O)
CONVERTBIN	= $(CONVERTSRC)/$(O)
COPBIN		= $(COPSRC)/$(O)
EBABIN		= $(EBASRC)/$(O)
EXPGETSEQBIN	= $(EXPGETSEQSRC)/$(O)
FROGBIN		= $(FROGSRC)/$(O)
GBIN		= $(GSRC)/$(O)
GAPBIN		= $(GAPSRC)/$(O)
GAP4BIN		= $(GAP4SRC)/$(O)
NIP4BIN		= $(NIP4SRC)/$(O)
GETMCHBIN	= $(GETMCHSRC)/$(O)
GETSCFFIBIN	= $(GETSCFFIELD)/$(O)
INDEXSEQLIBSBIN	= $(INDEXSEQLIBSSRC)/$(O)
INITEXPBIN	= $(INITEXPSRC)/$(O)
LOOKUPBIN	= $(LOOKUPSRC)/$(O)
MAKESCFBIN	= $(MAKESCFSRC)/$(O)
MISCBIN		= $(MISCSRC)/$(O)
OSPBIN	  	= $(OSPSRC)/$(O)
OSP4BIN	  	= $(OSP4SRC)/$(O)
SEQUTILSBIN	= $(SEQUTILSSRC)/$(O)
SEQLIBBIN	= $(SEQLIBSRC)/$(O)
SIP4BIN		= $(SIP4SRC)/$(O)
SPINBIN		= $(SPINSRC)/$(O)
COPYREADSBIN	= $(COPYREADSSRC)/$(O)
STADENBIN	= $(STADENSRC)/$(O)
TCLUTILSBIN	= $(TKUTILSSRC)/$(O)
TKUTILSBIN	= $(TKUTILSSRC)/$(O)
TEXTUTILSBIN	= $(TEXTUTILSSRC)/$(O)
TEDBIN		= $(TEDSRC)/$(O)
TRACECLIPBIN	= $(TRACECLIPSRC)/$(O)
TREVBIN		= $(TREVSRC)/$(O)
VEPEBIN		= $(VEPESRC)/$(O)
QCLIPBIN	= $(QCLIPSRC)/$(O)
MUTBIN		= $(MUTSRC)/$(O)


#------------------------------------------------------------------------------
# For files in remote directories that require recompilation, we should invoke
# the corresponding makefile.
# The standard compilation rules below are called once the cd has occured.
# The $(@F) removes the pathname match and hence removes recursion.
#
$(GAPBIN)/%:
	cd $(GAPBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(GBIN)/%.a:
	cd $(GBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)	
$(OSPBIN)/%.o: $(OSPSRC)/%.c
	cd $(OSPBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(OSP4BIN)/%.o: $(OSP4SRC)/%.c
	cd $(OSP4BIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(STADENBIN)/%.o: $(STADENSRC)/%.c
	cd $(STADENBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(STADENBIN)/%.o: $(STADENSRC)/%.f
	cd $(STADENBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(STADENBIN)/%.o: $(TEDSRC)/%.c
	cd $(STADENBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(STADENBIN)/%.a: $(STADENSRC)/%.f
	cd $(STADENBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(TEDBIN)/%.o: $(TEDSRC)/%.c
	cd $(TEDBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(CONVERTBIN)/%.o: $(CONVERTSRC)/%.c
	cd $(CONVERTBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(MISCBIN)/%.a:
	cd $(MISCBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(SEQLIBBIN)/%.o: $(SEQLIBSRC)/%.c
	cd $(SEQLIBBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(SEQLIBBIN)/%.a:
	cd $(SEQLIBBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(SEQUTILSBIN)/%.o: $(SEQUTILSSRC)/%.c
	cd $(SEQUTILSBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(SEQUTILSBIN)/%.a:
	cd $(SEQUTILSBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)
$(MUTBIN)/%.a: $(MUTSRC)/%.cpp
	cd $(MUTBIN)/$(S);$(MAKE) $(MFLAGS) $(@F)


# Automatically make output directories
%.dir:
	-mkdir -p $(@D) 2>/dev/null
	touch $@


#
# Files requiring simple C and FORTRAN compilation ($(SRC)/thing.c -> thing.o)
# in the local directory.
# 
# 7/1/99 johnt - added abstractions to support Visual C++, and mkdef command
# to build Windows DEF files
%.o:	%.c 
	$(CC) $(CFLAGS) $(CDEFS) $(COBJFLAG)$@ -c $<
	$(MKDEFC) $(MKFLAGS) $@

%.o:	%.f
	$(F77) $(FFLAGS) $(FDEFS) $(FOBJFLAG)$@ -c $<
	$(MKDEFC) $(MKFLAGS) $@

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEFS) $(COBJFLAG)$@ -c $<


#
# ``Interesting'' suffices that make needs to know about
#
# 7/1/99 johnt - added Windows suffixes
.SUFFIXES: .a .def .dll .cpp $(SUFFIXES)


#
# The UNIX library (archive) mechanism is fairly weak so the following
# procedure is used to generate libraries from the large fortran files.
#   * split the source file up into files which each contain
#     a single function or subroutine
#   * compile each of those files
#   * insert the object files individually into the library
# Each library is built by using a temporary directory.
# The fortran compilations must be done sequentially in order to
# avoid filling up the /tmp directory with compiler debugging information.
#
%.a: %.f
	SRCDIR=`pwd`; \
	TMPDIR=/tmp/staden$$$$; \
	if test ! -d $$TMPDIR; \
	then mkdir $$TMPDIR; \
	fi; \
	cd $$TMPDIR; \
	rm -f *.f *.o; \
	fsplit $$SRCDIR/`expr $< : '.*/\(.*\)' \| $<`; \
	$(F77) $(FFLAGS) -c *.f; \
	rm -f *.f; \
	rm -f $$SRCDIR/$@; \
	ar rcv $$SRCDIR/$@ *.o; \
	$(RANLIB) $$SRCDIR/$@; \
	rm -f *.o; \
	cd $$SRCDIR; \
	rm -rf $$TMPDIR

#------------------------------------------------------------------------------
# Default "distsrc" rule. This copies files from the source tree into a
# source distribution directory defined in $(DISTSRC). The defaults are
# sensible, but it is expected that this will be overridden. Relies on
# GNU make.
.PHONY: distsrc_dirs
distsrc_dirs:
	-mkdir -p $(DISTSRC)/$(SUBDIR)

distsrc: DIRNAME=$(DISTSRC)/$(SUBDIR)

#------------------------------------------------------------------------------
# Depend rule should work on most systems except Solaris.
# We assume here that our object files are .c files. This is an ok assumption
# as makedepend does not quit when it cannot read a file (just issues warnings)
# and our fortran files do not use #include anyway (as that's non ANSI).

depend:
	@# Run makedepend on our sources
	-DEPEND_SRC=`echo $(DEPEND_OBJ:.o=.c) $(DEPEND_OBJ:.o=.cpp) \
	| sed 's/\([^ ]*\) */ $(subst /,\/,$(VPATH))\/\1/g'`; \
	touch ./dependencies.tmp; \
	makedepend -f ./dependencies.tmp -- $(CFLAGS) -- $$DEPEND_SRC 2>&-

	@# Remove system paths and strip out local paths if they exist in
	@# one of our -I<dir> locations
	sort < ./dependencies.tmp | uniq | sed -e 's; /usr/[^ ]*;;g' | \
	  sed -e 's;.*/\([^:]*\):;\1:;'  | \
	  egrep -v '^[^:]*:[     ]*$$' | \
	  sed -e 's#$(subst .,\.,$(SRCROOT))#$$(SRCROOT)#g' \
	      -e 's#: .*/staden_config.h#: $$(PWD)/staden_config.h#g' | \
	  egrep -v ': /' > dependencies.tmp2

	@# Copy the dependencies into the Makefile
	l=`egrep -n 'DO NOT DELETE' Makefile | head -1 | sed 's/:.*//'`; \
	([ "x$$l" != "x" ] && \
	( mv Makefile Makefile.bak; \
	  ( head -$$l Makefile.bak; \
	    egrep -v '^#' dependencies.tmp2 ) > Makefile; ) \
	) || echo 'No "# DO NOT DELETE" line found in Makefile'

	@# tidy up
	rm ./dependencies.tmp*
