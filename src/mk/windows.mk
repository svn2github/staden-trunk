# Makefile for compiling with GNU make under Windows using Microsoft Visual C++

#
# Debugging On/Off
#

#_DEBUG=1


#
# Override Inherited flags
#
DEFINES += /D_WINDOWS
TCLVERS  = 84
TKVERS   = 84
ITCLVERS = 33
ITKVERS  = 33


ifdef _DEBUG
DEBUG_SUFFIX   = d
SHLIB_OPTDEBUG = $(SHLIB_DEBUG)
CLDOPTDEBUG    = $(CLDDEBUG)
COPTDEBUG      = $(CDEBUG)
FOPTDEBUG      = $(FDEBUG)
else
DEBUG_SUFFIX   =
SHLIB_OPTDEBUG = $(SHLIB_OPT)
CLDOPTDEBUG    = $(CLDOPT)
COPTDEBUG      = $(COPT)
FOPTDEBUG      = $(FOPT)
DEFINES       += /DNDEBUG=1
endif


#
# Directories
#
TCLBIN    = $(SRCROOT)/src/windows/windows-binaries
TKBIN     = $(SRCROOT)/src/windows/windows-binaries
MATH_LIB =
X_LIB    =
XSRC     = $(SRCROOT)/tk8.4.0/xlib
XAW_LIB  =
TKWININC = $(SRCROOT)/tk8.4.0/win


#
# Compiling
#
MAKE = make
ABSOFT_PATH = e:/absoft


#TKINCDIR must come before ...windows\include so we get the TK version of some X include files !
INCLUDES_E  += $(TK_INC) -I$(TKWININC) -I$(SRCROOT)/windows/include

RC        = rc.exe
CC        = cl.exe /nologo
CXX       = $(CC)
CLD_PROG  = link.exe /MAP
CXXLD_PROG=$(CLD_PROG)
FLD       = $(CLD)
CDEBUG    = /D_DEBUG /MD$(DEBUG_SUFFIX) /Zi /Od /W3
COPT      = /MD$(DEBUG_SUFFIX) /Ox /O2
CXXFLAGS += /GX /TP
FDEBUG    =
FFLAGS    = -f -N15 -K
MKDEFC    = $(SRCROOT)/windows/mkdef/release/mkdef -c
RM        = rm
NO_SRS    = 1
COBJFLAG  = /Fo
LDEXEFLAG = /OUT:


#
# Library Exports
#
# Unfortunately, we need one definition for each DLL since during a build
# they have dependencies on each other.
#
G_DLL         = /DBUILDING_G_DLL
MISC_DLL      = /DBUILDING_MISC_DLL
TK_UTILS_DLL  = /DBUILDING_TK_UTILS_DLL
SEQ_UTILS_DLL = /DBUILDING_SEQ_UTILS_DLL


#
# Linking, always used dynamic CRT. /NODEFAULTLIB forces symbols in static
# libraries to be searched for in the dynamic CRT instead of libc.lib which
# is the default.
#
CRTLIBS       = /NODEFAULTLIB:libc.lib msvcrt$(DEBUG_SUFFIX).lib
LINK_PATHFLAG = /LIBPATH:


#
# Linking a dynamic library
#
SHLIB_LD        = link.exe
SHLIB_LDXX      = $(SHLIB_LD)
SHLIB_LDFLAGS   = /MAP /DLL $(SHLIB_OPTDEBUG) $(CRTLIBS) $(LINK_PATHFLAG)$(SRCROOT)/windows/windows-binaries $(LINK_PATHFLAG)$(L) /DEF:$*.def gdi32.lib user32.lib comdlg32.lib advapi32.lib wsock32.lib
SHLIB_LDXXFLAGS = $(SHLIB_LDFLAGS)
SHLIB_PREFIX    =
SHLIB_SUFFIX    = .dll
SHLIB_SONAME    =
SHLIB_OUTFLAG   = /OUT:
SHLIB_OPT       =
SHLIB_DEBUG     = /DEBUG
SHLIB_CFLAGS    =
SHLIB_CXXFLAGS  =


#
# Linking an executable
#

LINK_LIBFLAG  =
LIB_EXT       = .lib
REPESTACK     = /STACK:0x200000
CLDDEBUG      = /DEBUG
CLDOPT        = $(CRTLIBS) /MACHINE:IX86
CLDFLAGS_S   += $(LINK_PATHFLAG)$(SRCROOT)/windows/windows-binaries
EXE_SUFFIX    = .exe
EXTRA_LIBS    = user32.lib gdi32.lib comdlg32.lib advapi32.lib wsock32.lib
F77_DEP       = ntf77.lib


# WinStash Icon
#RESFILE      = $(O)/resource.res
#SUBSYSTEMWIN = /SUBSYSTEM:windows $(RESFILE)
SUBSYSTEMWIN = /SUBSYSTEM:windows


#Options for Console EXE mode - use setargv.obj to get wildcard expansion
SUBSYSTEMCONSOLE = /SUBSYSTEM:console setargv.obj



# Supplies missing XWindows calls needed by windows
TKUTILS_EXTRAS = $(TKUTILSBIN)/tkWinX.o


#$(RESFILE): resource.rc
#   $(RC) /fo $(RESFILE) resource.rc


CLEANCOMMAND=-$(RM) $(O)/*.o $(O)/*.a $(O)/*.map $(O)/*.def $(O)/*.exp $(O)/*.lib $(O)/*.ilk $(O)/*.pdb $(O)/*.exe
CLEANLIBSCOMMAND=-rm -f $(PROGLIBS:.dll=.*)

#
# Macros to require, and build windows module definition file
#
DEF_FILE = $(L)/$(SHLIB_PREFIX)$(LIBS).def
MKDEFL    = $(SRCROOT)/windows/mkdef/release/mkdef -l
