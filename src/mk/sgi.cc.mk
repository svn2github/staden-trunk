#
# SGI system Makefile using SGI's cc. Rename this file to sgi.mk.
# GNU make 3.64 or newer must be used on this system.
#

/* DEFINES += -DBIG_ENDIAN -DNOSTRDUP */
COPT	+= -Olimit 1200

# ABI - new 32-bit calling convention.
# Note that the optimiser (-O2) for -o32 has a bug in it which
# causes io_lib to compile wrongly.
ABI=-n32

NO_SRS=1

# You may also wish to add -mips2 on later SGI systems to make user of the
# newer processor specific instructions. Similarly for fortran.
CFLAGS	+= -ansi $(ABI)
CLDFLAGS+= $(ABI)

FOPT	+= -Olimit 2500
FFLAGS	+= -old_rl $(ABI)
FLDFLAGS+= $(ABI)

# Using fortran to link appears to cause dynamic linking of the fortran
# libraries, even if we specify them explicitly. We'll use cc instead.
#FLD_PROG = $(CC)
#LIBSF_S	 = -lF77 -lU77 -lI77 -lm -lisam
FLD_PROG = $(F77)
LIBSF_S	 = -lm

# Ranlib functionality has been folded into the ar utility.
RANLIB	 = /bin/true

# Linking a dynamic library
# The SGI ld command treats -o in a very odd way when creating dynamic
# libraries.

# The argument after the -o is treated exactly as the library SONAME, which
# means that subsequent links with this library use this too. SONAMES
# containing a slash are considered as pathnames, which changes the rld
# behaviour to ignore rpath, LD_LIBRARY_PATH, _RLD_ROOT and the default
# library paths when searching for the library.

# The solution is to use -soname to reset the soname. For instance:
# ld -shared -o sgi-binaries/libfred.so -soname libfred.so (... etc)
# Hence the provision os SHLIB_SONAME in the Makefiles. The assumption made
# here is that all dynamic libraries are placed in $(L).
SHLIB_CFLAGS		=
SHLIB_LD		= ld
SHLIB_LDFLAGS		= $(SHLIB_STRIP) $(ABI) \
			  -update_registry $(L)/so_locations \
			  -shared -rdata_shared -o
SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .so
SHLIB_SONAME		= $(@:$(L)/%=-soname %)
CLDFLAGS_E	       += -Wl,-rpath,$(INSTALLLIB)
FLDFLAGS_E	       += -Wl,-rpath,$(INSTALLLIB)
