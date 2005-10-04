#
# Defines for a DEC Alpha running Digital Unix
# GNU make 3.64 or newer must be used on this system.
#

PNG_INC=-I/usr/local/include

# DEFINES	+= -DLITTLE_ENDIAN

#COPT	+= -Olimit 1000
#CFLAGS	+= -std1 -verbose -trapuv #-ieee_with_inexact
#CFLAGS  += -msg_enable level3
#CC = /nfs/team71/psg/dgm/sys/alpha/bin/gcc-3.2
#GCCWARNINGS = -Wimplicit -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat \
#	      -Wstrict-prototypes
#CFLAGS	+= $(GCCWARNINGS) -gcoff

# For DEC cc
CFLAGS += -w0 -warnprotos


#CXX= /nfs/team71/psg/dgm/sys/alpha/bin/g++-3.2
#CXXFLAGS += -I$(SRCROOT)/stlport/alpha
CXX=cxx
CXXFLAGS += -D__USE_STD_IOSTREAM
#CXXFLAGS += -std strict_ansi

#CXX = cxx
#CXXFLAGS += -std strict_ansi

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
#F77=g77
#FLD_PROG=g77
#F77_DEP = /usr/local/lib/gcc-lib/alpha-dec-osf5.1/3.0.3/libg2c.a
#F77_DEP	 = $(T)for.a $(T)util.a $(T)Ufor.a $(T)Futil.a -lm\
#	   $(T)for.a $(T)util.a $(T)Ufor.a $(T)Futil.a -lm -lc_r
#FLD_PROG = $(CC)

F77=f77
FLD_PROG=$(CC)
F77_DEP=

EXTRA_LIBS += -ldnet_stub -lots -lmld

# Linking a dynamic library.
SHLIB_CFLAGS  		=
SHLIB_LD      		= ld

#ifeq ($(USER),pubseq)
#SHLIB_LDFLAGS 		= $(SHLIB_STRIP) \
#			  -check_registry /usr/shlib/so_locations \
#			  -check_registry $(L)/so_locations \
#			  -update_registry $(L)/so_locations \
#			  -shared -expect_unresolved "*" -lots -o
#else
#SHLIB_LDFLAGS 		= $(SHLIB_STRIP) \
#			  -check_registry /usr/shlib/so_locations \
#			  -check_registry $(L)/so_locations \
#			  -shared -expect_unresolved "*" -lots -o
#endif

# Don't bother with maintaining our own so_locations file. It does reduce
# library relocation, but the overhead is minimal compared to the typical
# length of time that the apps run for.
SHLIB_LDFLAGS 		= -L -L$(L) -L/usr/shlib $(SHLIB_STRIP) \
			  -shared -expect_unresolved "*" -lots -o

SHLIB_PREFIX		= lib
SHLIB_SUFFIX		= .so
SHLIB_SONAME		=
#CLDFLAGS_E	       += -rpath '$${STADLIBBIN}'
#FLDFLAGS_E	       += -rpath '$${STADLIBBIN}'

#SHLIB_LDXX		= g++
SHLIB_LDXX		= $(CXX)
SHLIB_LDXXFLAGS		= -L$(L) $(SHLIB_STRIP) \
			  -shared -lots -o

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
