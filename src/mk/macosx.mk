#
# Makefile for MacOS X 10.1
#

# Complain about everything sensible, but not the stylistic complaints
# introduced by -Wall
COPT += -Wuninitialized 
GCCWARNINGS = -Wimplicit -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat
CFLAGS	+= $(GCCWARNINGS)
CXX=c++
#CXXFLAGS += -I$(SRCROOT)/stlport/macosx

# X is installed in /usr/X11R6
XBIN	= /usr/X11R6/lib
INCLUDES_E  += -I/usr/X11R6/include

# Tcl/Tk are a Frameworks...
#TCLSRC=/Library/Frameworks/Tcl.framework/Headers
#TKSRC=/Library/Frameworks/Tk.framework/Headers


#TKUTILS_EXTRAS = $(TKUTILSBIN)/tkMacX.o

# The MacOSX-aqua version is based around tcl8.4. We also assume that this
# has been installed on the system in the normal place (/Library/Frameworks).
#TCLSRC=/Users/jkb/latest_tcl_source/tcl8.4b2/generic
#TKSRC=/Users/jkb/latest_tcl_source/tk8.4b2/generic
#TCLVERS = 8.4
#TKVERS = 8.4
#TCL_LIB=-F$(SRCROOT)/ftp/tcl-cvs/macosx/build -framework Tcl
#TK_LIB=$(TCL_LIB) -F$(SRCROOT)/ftp/tk-cvs/macosx/build -framework Tk
#TCL_LIB=-ltcl8.4
#X_LIB_S=-L/usr/X11R6/lib
#TK_LIB=$(TCL_LIB) -ltk8.4 $(X_LIB)

#TCL_DEP=$(TCL_LIB)
#TK_DEP=$(TK_LIB)
IOLIB_DEP=$(IOLIB_LIB)
MISC_DEP=$(MISC_LIB)
TKUTILS_DEP=$(TKUTILS_LIB)

# Hack for now as things are wrongly linking against TK_UTILS to work around
# limitations in windows linking. This causes the Mac problems as
# linking against TK_UTILS then requires linking against Tk.
EXTRA_LIBS=$(TK_LIB)

# Similarly for png. I do not understand why linking tk_utils againgst
# png means that all apps linked against tk_utils then need to
# explicitly link against png, but it appears they do.
TKUTILS_LIB_E=$(PNG_LIB)
EXTRA_LIBS += $(PNG_LIB)

FLD=$(CLD)

# Fortran code only really still exists in Gap4's legacy and repe.
# We use f2c to convert and rename. The renaming is so that the old %.f rules
# still apply on the other systems.
%_f2c.c: %.f
	../f2c/f2c $*.f
	mv $*.c $*_f2c.c
F77_DEP=-L../../lib.third/macosx-binaries -lF
F77_INC = -I../f2c

GAP4_LEGACY	= legacy_f2c.o
REPE_OBJ	= repe_f2c.o

# Dynamic linking
SHLIB_CFLAGS		= -fPIC
SHLIB_LD		= cc
SHLIB_LDFLAGS		= $(SHLIB_STRIP) -L$(L) -dynamiclib -bind_at_load -o
SHLIB_LDXX		= c++
SHLIB_LDXXFLAGS		= $(SHLIB_LDFLAGS)
SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .dylib
SHLIB_SEARCH_FLAGS	=
SHLIB_SONAME		=
