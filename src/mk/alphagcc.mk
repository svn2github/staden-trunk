#
# Defines for a DEC Alpha running Digital Unix
# GNU make 3.64 or newer must be used on this system.
#

# DEFINES	+= -DLITTLE_ENDIAN

#COPT	+= -Olimit 1000
#CFLAGS	+= -std1 -verbose -trapuv #-ieee_with_inexact
CC=gcc
CFLAGS += -fbounds-checking
CLDFLAGS += -fbounds-checking

# Link the fortran programs using cc. We now have to specify the fortran
# libraries to do this, but it solves problems of fortran complaining about
# our use of a C main.
# Also note that the Digital Unix linker seems incapable of linking some
# libraries static and others static. For this reason we specify the .a files
# explicitly to force only these to be static.
# The -qlc_r is needed (-ql => non-erroring) when linking against libfor.a
# from Digital Unix 4.0, but running on Digital Unix 3.0.
# 
T = /usr/lib/lib
#LIBSF_E	+= $(T)for.a $(T)util.a $(T)Ufor.a $(T)ots.a $(T)Futil.a -lm\
#	   $(T)for.a $(T)util.a $(T)Ufor.a $(T)ots.a $(T)Futil.a -lm
LIBSF_E	+= $(T)for.a $(T)util.a $(T)Ufor.a $(T)Futil.a -lm\
	   $(T)for.a $(T)util.a $(T)Ufor.a $(T)Futil.a -lm -lc_r
FLD_PROG = $(CC)

EXTRA_LIBS += -ldnet_stub -lots -lmld

# Linking a dynamic library.
SHLIB_CFLAGS  		=
SHLIB_LD      		= ld
SHLIB_LDFLAGS 		= $(SHLIB_STRIP) \
			  -check_registry /usr/shlib/so_locations \
			  -check_registry $(L)/so_locations \
			  -update_registry $(L)/so_locations \
			  -shared -expect_unresolved "*" -lots -o
SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .so
SHLIB_SONAME		=
CLDFLAGS_E	       += -Wl,-rpath,$(INSTALLLIB)
FLDFLAGS_E	       += -Wl,-rpath,$(INSTALLLIB)

# To link statically instead of dynamically uncomment these lines.
#L = $(SRCROOT)/lib-static/$(O)
#SHLIB_SUFFIX=.a
#SHLIB_LD=ar
#SHLIB_LDFLAGS=rv

# Tmp stuff
#TCLBIN	 = /home6/jkb/ftp/tcl7.6/unix
#TCLSRC	 = /home6/jkb/ftp/tcl7.6/generic
#TKBIN	 = /home6/jkb/ftp/tk4.2/unix
#TKSRC	 = /home6/jkb/ftp/tk4.2/generic

#COPTDEBUG	= $(COPT)
#FOPTDEBUG	= $(FOPT)
