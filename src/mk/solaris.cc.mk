COPT		= -O -g
FOPT		= -O2 -g
#SHLIB_STRIP	= -s

#
# Solaris 2.5
#
#DEFINES        += -DBIG_ENDIAN -DSYSV_REGEX -DIMAGEDISPLAY -DNOSTRDUP

# __EXTENSIONS__ gets things like popen and pclose. These may also be picked
# up using -D_XOPEN_SOURCE, but that blocks things like timeval. Why?
DEFINES		+= -D__EXTENSIONS__

# Used by gap. Only define this when we also define -DIMAGEDISPLAY
TT_OBJS		= $(GAPBIN)/tooltalklib.o

# The X libraries are housed in an abnormal place
XBIN		= /usr/openwin/lib
TT_BIN		= /usr/openwin/lib
XSRC		= /usr/openwin/include
TTSRC		= /usr/openwin/include
X_LIB_S		= -R/usr/openwin/lib
XAW_LIB_S	= -R/usr/openwin/lib
TT_LIB_S	= -R/usr/openwin/lib
TT_LIBRARY	= -ltt

# Ranlib functionality is folded into ar
RANLIB		= /bin/true

#------------------------------------------------------------------------------
#24/5/00 johnt - added Biolims support whe BIOLIMS defined
# BIOLIMS=1

#------------------------------------------------------------------------------
# Biolims stuff. Uncommenting this will enable biolims support. You'll also
# need to edit global.mk to uncomment the "BIOLIMS=1" line. (For technical
# reasons it doesn't work if defined here.)
# Also see io_lib/options.mk
BIOLIMS_LIB  = $(LINK_LIBFLAG)biolimsio$(LIB_EXT)
BIOGAP_LIB   = $(LINK_LIBFLAG)biolimsgap$(LIB_EXT)


#------------------------------------------------------------------------------
# The sysv regex code seems to be buggy and cause crashes. The BSD re_exec code
# doesn't work too well on Solaris 2.4, so in that case we have our own library
# with regex.o copied out of the Solaris 2.5 libc.a.
# libour_regex.a needs to be created in lib/solaris-binaries, or just remove
# this line if using Solaris 2.5 onwards.
# EXTRA_LIBS=-lour_regex

#------------------------------------------------------------------------------
# Local stuff. We use gcc v2.3.3 and SparcWorks Fortran v2.0.1.
# Points to note:
# 1. When linking with fortran, we need to link the gcc library.
# 2. We have no copy of the dynamic gcc library so no need for -Bstatic
# 3. Fortran has dynamic libraries, but we must link static to distribute
# 4. Fortran ALWAYS adds a -lF77 after the command line that we specify.
#    Hence we must have -Bstatic at the end of the Fortran link line, which
#    in turn means we have to link dynamically all the system stuff first.
CC	    = cc
CXX	    = CC
#CC	    = gcc
#GCCWARNINGS = -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat
#GCCWARNINGS = \
#	-W \
#	-Wimplicit-function-declaration \
#	-Wmain \
#	-Wswitch \
#	-Wcomment \
#	-Wtrigraphs \
#	-Wformat \
#	-Wchar-subscripts \
#	-Wpointer-arith
#
#	-Wuninitialized

#CFLAGS     += -ansi -pedantic $(GCCWARNINGS)
CFLAGS      += -Xa
#GCCLIB	    = /usr/local/lib/gcc-lib/sparc-sun-solaris2.6/2.95.2/
#EXTRA_LIBS += -L$(GCCLIB) -lgcc -ldl -lnsl -lsocket
EXTRA_LIBS += -ldl -lnsl -lsocket
LIBSF_E    += -lm -lc -ldl -Bstatic
#LIBSF_E    += -lm -lc -ldl


# Linking a dynamic library.
#SHLIB_CFLAGS  		= -fPIC
SHLIB_CFLAGS		=
SHLIB_LD      		= cc
SHLIB_LDFLAGS 		= -L$(L) $(SHLIB_STRIP) $(LINK_PATHFLAG)$(L) -G -o
SHLIB_LDXX     		= CC
SHLIB_LDXXFLAGS 	= $(SHLIB_LDFLAGS)
SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .so
SHLIB_SONAME		=
SHLIB_DEP		= -Bdynamic -ldl -lresolv
CLDFLAGS_E	       += -R $(INSTALLLIB)
FLDFLAGS_E	       += -R $(INSTALLLIB)
CXX_DEP			= -L/opt/SUNWspro/WS6U1/lib/ -Bstatic -lCstd -Bdynamic -lCrun
F77_DEP			= -L/opt/SUNWspro/WS6U1/lib/ -Bstatic -lF77 -lM77 -lm -lsunmath -Bdynamic

# To link statically instead of dynamically uncomment these lines.
#L = $(SRCROOT)/lib-static/$(O)
#SHLIB_SUFFIX=.a
#SHLIB_LD=ar
#SHLIB_LDFLAGS=rv
