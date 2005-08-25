#
# Defines for a RedHat Linux system.
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
GCCWARNINGS = -Wimplicit -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat \
	      -Wstrict-prototypes
CFLAGS	+= $(GCCWARNINGS)

# For large file support (off_t by default is only 32-bit)
CFLAGS += -D_FILE_OFFSET_BITS=64
CXXFLAGS += -D_FILE_OFFSET_BITS=64

#CXXFLAGS += -I$(SRCROOT)/stlport/linux

# X is installed in /usr/X11R6
XBIN	= /usr/X11R6/lib
# For AMD64 based systems use XBIN=/usr/X11R6/lib64

#CLDFLAGS_E	+= -shared
#FLDFLAGS_E	+= -shared

# Dynamic linking
# For AMD64 based systems try -fPIC instead.
SHLIB_CFLAGS		= -fpic

# pic with g77 generates opcodes unknown by the 386 assembler. Odd
# It doesn't appear to be needed (but at a speed cost?).
# For AMD64 based systems try -fPIC instead.
#SHLIB_FFLAGS		= -fpic
SHLIB_LD		= cc
SHLIB_LDFLAGS		= -L$(L) $(SHLIB_STRIP) -shared -o
SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .so
SHLIB_SEARCH_FLAGS	=
SHLIB_SONAME		=
EXTRA_LIBS		+= -ldl

SHLIB_LDXX		= $(CXX)
SHLIB_LDXXFLAGS		= $(SHLIB_LDFLAGS)
