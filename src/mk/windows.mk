#
# Defines for a windows machine running MinGW
# GNU make 3.64 or newer must be used on this system.
#

/* DEFINES += -DLITTLE_ENDIAN */

# Use g77 for compiling our fortran
F77	 = g77
# Ewe:
F77_DEP  = -lg2c
LIBSF	= $(F77_DEP)
# Jura
# F77_DEP	 = /usr/lib/gcc-lib/i386-glibc21-linux/egcs-2.91.66/libg2c.a

# Complain about everything sensible, but not the stylistic complaints
# introduced by -Wall
CC = gcc
CXX = g++
COPT += -Wuninitialized 
GCCWARNINGS = -Wimplicit -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat \
	      -Wstrict-prototypes
CFLAGS	+= $(GCCWARNINGS)

CFLAGS += -I$(SRCROOT)/windows/include -I$(TKSRC)/../win

#CXXFLAGS += -I$(SRCROOT)/stlport/linux

XSRC=$(TKSRC)/../xlib

ZLIB_LIB     = $(ZLIB_LIB_S) $(LINK_LIBFLAG)zlib1$(LIB_EXT) $(ZLIB_LIB_E)
TCLVERS		= 84
TKVERS		= 84
X_LIB=		-lgdi32

#CLDFLAGS_E	+= -shared
#FLDFLAGS_E	+= -shared

TKUTILS_EXTRAS=$(TKUTILSBIN)/tkWinX.o

# Dynamic linking
#SHLIB_CFLAGS		= -fpic
# pic with g77 generates opcodes unknown by the 386 assembler. Odd
# It doesn't appear to be needed (but at a speed cost?).
#SHLIB_FFLAGS		= -fpic
SHLIB_LD		= $(CC)
SHLIB_LDFLAGS		= --enable-auto-import -L$(L) $(SHLIB_STRIP) -shared -lcomdlg32 -o
SHLIB_PREFIX		= 
SHLIB_SUFFIX		= .dll
SHLIB_SEARCH_FLAGS	=
SHLIB_SONAME		=
#EXTRA_LIBS		+= -ldl

SHLIB_LDXX		= $(CXX)
SHLIB_LDXXFLAGS		= $(SHLIB_LDFLAGS)
