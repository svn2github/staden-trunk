#-----------------------------------------------------------------------------
# NB: This Makefile is way out of date - it won't work with the dynamic
# library setup in use now. However it may be instructive on compilation
# for older Solaris machines.
#-----------------------------------------------------------------------------

#
# For Solaris 2.[234]
#
# We link as much as possible dynamic so that our binaries will attempt to run
# on as many solaris systems as possible. However the Athena Widget library is
# an optionally installed package, so we should not require the user to have it
# installed. Hence we need to statically link with the X libraries. This in
# turn implies dynamic linking with a LOT of other libraries. Also, for Fortran
# we do not wish to have the bother of distributing the run time libraries, so
# we end with a -Bstatic. This then implies we need to have linked various C
# and system libraries (-lc -lsys) explicitly dynamic before hand. All rather
# messy!
#

DEFINES	  = -DSYSV_REGEX -DIMAGEDISPLAY -DNOSTRDUP -DBIG_ENDIAN
# Only define this is we also define -DIMAGEDISPLAY
TT_OBJS	  = $(BAPBIN)/tooltalklib.o

# Add a few extra places onto our library search path
SPRO	  = /opt/SUNWspro/SC2.0.1
GCCLIB	  = /opt/cygnus-sol2-1.1/lib/gcc-lib/sparc-sun-solaris2/cygnus-2.0.2
LDFLAGS_X = -L$(SPRO) -L/usr/ccs/lib -L$(GCCLIB)
CLDFLAGS_X= -L$(SPRO) -L/usr/ccs/lib -L$(GCCLIB)

# Link the X libraries - they're housed in an abnormal place too
OPENWIN   = /usr/openwin
XLIBS_X	  = -L$(OPENWIN)/lib -R$(OPENWIN)/lib
CLIBS_X	  = -L/usr/local/lib 
# -lXext only needed on Solaris 2.4
XLIBS	 += -lXext -lX11 -lce -ltt -Bdynamic
INCLUDES_X= -I$(OPENWIN)/include

# Gcc library
EXTRA_LIBS += -lform -lgcc

# Statically link the rest that Fortran adds.
#FLIBS_X   = -Bstatic
FLIBS	 += -Bstatic -lV77 -lF77 -lM77 -Bdynamic

CC	= gcc
CLD	= $(F77)
GCCWARNINGS = -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat
CFLAGS += -ansi -pedantic $(GCCWARNINGS)

FFLAGS += -C

RANLIB=/bin/true

#----------------------------------------------------------------------------
#
# For Solaris 2.1
#
# Not the optimal solution - and doesn't deal with the issues described above,
# but this is known to have worked on a Solaris 2.1 system.
#
#---
#
#DEFINES	  = -DSYSV_REGEX -DBIG_ENDIAN
#
#SPRO	  = /opt/SUNWspro/SC2.0.1
#LDFLAGS_X = -L$(SPRO) -L/usr/ccs/lib -L/usr/ucblib -Bstatic
#
#OPENWIN   = /usr/openwin
#XLIBS_X	  = -L$(OPENWIN)/lib
#INCLUDES += -I$(OPENWIN)/include
#
#CLIBS	 += -lm -lucb -lelf -Bdynamic -ldl -Bstatic -lsocket -lnsl -lintl -lgcc
#
#FLIBS	 += -lV77 -lF77 -lM77
#
#CC	= gcc
#GCCWARNINGS = -Wreturn-type -Wunused -Wswitch -Wcomment -Wformat
#CFLAGS += -ansi -pedantic $(GCCWARNINGS)
#
#FFLAGS += -C
#
#RANLIB=/bin/true
