#
# Defines for a RedHat Linux system.
# GNU make 3.64 or newer must be used on this system.
#

#PNG_LIB	     = /usr/lib/libpng12.so.0
#ZLIB_LIB     = /usr/lib/libz.so.1

COPT		= -O2 -DNDEBUG
FOPT		= -O2 -DNDEBUG

DEFINES += -DSP_LITTLE_ENDIAN


# Use g77 for compiling our fortran
F77	 = icc
# Ewe:
GAP4_LEGACY	= legacy_f2c.o
#F77_DEP  = -lg2c
LIBSF	= $(F77_DEP)
# Jura
# F77_DEP	 = /usr/lib/gcc-lib/i386-glibc21-linux/egcs-2.91.66/libg2c.a

# Complain about everything sensible, but not the stylistic complaints
# introduced by -Wall
CC = icc
CXX = icpc
#COPT += -Wuninitialized 
#GCCWARNINGS = -Wimplicit -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat \
	      -Wstrict-prototypes
#CXXFLAGS += -I$(SRCROOT)/stlport/linux

# X is installed in /usr/X11R6
XBIN	= /usr/X11R6/lib

#CLDFLAGS_E	+= -shared
#FLDFLAGS_E	+= -shared

# Dynamic linking
SHLIB_CFLAGS		= -fpic
SHLIB_CXXFLAGS		= -fpic
# pic with g77 generates opcodes unknown by the 386 assembler. Odd
# It doesn't appear to be needed (but at a speed cost?).
#SHLIB_FFLAGS		= -fpic
SHLIB_LD		= $(CC)
SHLIB_LDFLAGS		= -L$(L) $(SHLIB_STRIP) -shared -o
SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .so
SHLIB_SEARCH_FLAGS	=
SHLIB_SONAME		=
EXTRA_LIBS		+= -ldl

SHLIB_LDXX		= $(CXX)
SHLIB_LDXXFLAGS		= $(SHLIB_LDFLAGS)
