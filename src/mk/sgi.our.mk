#
# SGI system Makefile using gcc, but SGI's f77.
# GNU make 3.64 or newer must be used on this system.
#

# You may also wish to add -mips2 on later SGI systems to make user of the
# newer processor specific instructions. Similarly for fortran.
#CC	= gcc
#CFLAGS	+= -ansi -mabi=n32 
#FOPT	+= -Olimit 2500
#FFLAGS	+= -old_rl -n32 #-mips2
CC = cc
CFLAGS = -n32 -mips3

# Using fortran to link appears to cause dynamic linking of the fortran
# libraries, even if we specify them explicitly. We'll use cc instead.
FLD_PROG = $(CC)
#LIBSF_S	 = -L/usr/lib32 -lF77 -lU77 -lI77 -lm -lisam
LIBSF_S	 = -lftn -lm

#F77	= g77
#FLD_PROG  = g77
#FFLAGS	+= -old_rl -mabi=n32

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
SHLIB_LDFLAGS		= $(SHLIB_STRIP) -n32 \
			  -update_registry $(L)/so_locations \
			  -shared -rdata_shared -o
SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .so
SHLIB_SONAME		= $(@:$(L)/%=-soname %)
CLDFLAGS_E	       += -Wl,-rpath,$(INSTALLLIB)
FLDFLAGS_E	       += -Wl,-rpath,$(INSTALLLIB)
