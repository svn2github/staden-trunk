# SunOS 4.1 system definitions
# DEFINES	+= -DBIG_ENDIAN -DNOMEMMOVE -DNOSTRERROR

# Fortran requires an increase in the symbol table space
FFLAGS	+= -Nx500

# We have X11R6 installed, but some systems may have openwindows instead,
# in which case set XBIN to /usr/openwin/lib.
# Due to this we link X statically. Under X11R6 this also means that we
# need to link the SM and ICE libraries (part of X).
XBIN		= /usr/X11R6/lib
X_LIB_S		= -Bstatic
X_LIB_E		= -lSM -lICE -Bdynamic
XAW_LIB_S	= -Bstatic
XAW_LIB_E	= -lSM -lICE -Bdynamic

# Gcc cannot support linking partially static and partially dynamic.
# Additionally it doesn't like the standard -Bstatic notation, but rather uses
# -static. We need to link X programs statically (see above) and so must use
# fortran to perform the link phase. Hopefully this isn't needed if you're
# recompiling with the static X bits removed. A consequence of using F77 for
# linking is that it then need explicit linking of the gcc library.
# If you're using SPARCworks cc compiler then I expect this is so much easier.
#
# Also note that Fortran has dynamic libraries, but we wish distribute with
# these linked static. We do this by specifying the system stuff as dynamic
# and then let f77 link the rest as static.
CLD_PROG 	= $(F77)
EXTRA_LIBS     += -lgcc
LIBSC_E	       += -lm -lc -ldl -Bstatic
LIBSF_E	       += -lm -lc -ldl -Bstatic

# Linking a dynamic library.
SHLIB_CFLAGS  		= -fPIC
SHLIB_FFLAGS  		= -PIC
SHLIB_LD      		= ld
SHLIB_LDFLAGS 		= -assert pure-text -o
SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .so.1.0
SHLIB_SEARCH_FLAGS	=
SHLIB_SONAME		=

# To link statically instead of dynamically uncomment these lines.
#L = $(SRCROOT)/lib-static/$(O)
#SHLIB_SUFFIX=.a
#SHLIB_LD=ar
#SHLIB_LDFLAGS=rv

#-----------------------------------------------------------------------------
# Local definitions. We use SPARCworks fortran v1.0 with gcc v2.7.0.
CC		= gcc
GCCWARNINGS	= -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat
CFLAGS	       += -ansi -pedantic $(GCCWARNINGS)

# Sun version numbers are handled differently when compiling tcl/tk.
TCLVERS		= 80
TKVERS		= 80
