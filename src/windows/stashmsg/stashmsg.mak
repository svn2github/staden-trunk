# Microsoft Developer Studio Generated NMAKE File, Based on stashmsg.dsp
!IF "$(CFG)" == ""
CFG=stashmsg - Win32 Debug
!MESSAGE No configuration specified. Defaulting to stashmsg - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "stashmsg - Win32 Release" && "$(CFG)" !=\
 "stashmsg - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "stashmsg.mak" CFG="stashmsg - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "stashmsg - Win32 Release" (based on\
 "Win32 (x86) Dynamic-Link Library")
!MESSAGE "stashmsg - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "stashmsg - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "messages.rc" "$(OUTDIR)\stashmsg.dll"

!ELSE 

ALL : "messages.rc" "$(OUTDIR)\stashmsg.dll"

!ENDIF 

CLEAN :
	-@erase "$(INTDIR)\stadenmsg.res"
	-@erase "$(OUTDIR)\stashmsg.dll"
	-@erase "$(OUTDIR)\stashmsg.exp"
	-@erase "$(OUTDIR)\stashmsg.lib"
	-@erase "messages.rc"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)\stashmsg.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
CPP_OBJS=.\Release/
CPP_SBRS=.

.c{$(CPP_OBJS)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(CPP_OBJS)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(CPP_OBJS)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(CPP_SBRS)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(CPP_SBRS)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(CPP_SBRS)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

MTL=midl.exe
MTL_PROJ=/nologo /D "NDEBUG" /mktyplib203 /o NUL /win32 
RSC=rc.exe
RSC_PROJ=/l 0x409 /fo"$(INTDIR)\stadenmsg.res" /d "NDEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\stashmsg.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:windows /dll /incremental:no\
 /pdb:"$(OUTDIR)\stashmsg.pdb" /machine:I386 /out:"$(OUTDIR)\stashmsg.dll"\
 /implib:"$(OUTDIR)\stashmsg.lib" /noentry 
LINK32_OBJS= \
	"$(INTDIR)\stadenmsg.res"

"$(OUTDIR)\stashmsg.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "stashmsg - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "messages.rc" "$(OUTDIR)\stashmsg.dll"

!ELSE 

ALL : "messages.rc" "$(OUTDIR)\stashmsg.dll"

!ENDIF 

CLEAN :
	-@erase "$(INTDIR)\stadenmsg.res"
	-@erase "$(OUTDIR)\stashmsg.dll"
	-@erase "$(OUTDIR)\stashmsg.exp"
	-@erase "$(OUTDIR)\stashmsg.ilk"
	-@erase "$(OUTDIR)\stashmsg.lib"
	-@erase "$(OUTDIR)\stashmsg.pdb"
	-@erase "messages.rc"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MTd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)\stashmsg.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
CPP_OBJS=.\Debug/
CPP_SBRS=.

.c{$(CPP_OBJS)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(CPP_OBJS)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(CPP_OBJS)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(CPP_SBRS)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(CPP_SBRS)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(CPP_SBRS)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

MTL=midl.exe
MTL_PROJ=/nologo /D "_DEBUG" /mktyplib203 /o NUL /win32 
RSC=rc.exe
RSC_PROJ=/l 0x409 /fo"$(INTDIR)\stadenmsg.res" /d "_DEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\stashmsg.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:windows /dll /incremental:yes\
 /pdb:"$(OUTDIR)\stashmsg.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)\stashmsg.dll" /implib:"$(OUTDIR)\stashmsg.lib" /pdbtype:sept\
 /noentry 
LINK32_OBJS= \
	"$(INTDIR)\stadenmsg.res"

"$(OUTDIR)\stashmsg.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(CFG)" == "stashmsg - Win32 Release" || "$(CFG)" ==\
 "stashmsg - Win32 Debug"
SOURCE=.\messages.mc

!IF  "$(CFG)" == "stashmsg - Win32 Release"

InputPath=.\messages.mc

"messages.rc"	 : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	mc messages.mc

!ELSEIF  "$(CFG)" == "stashmsg - Win32 Debug"

InputPath=.\messages.mc

"messages.rc"	 : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	mc messages.mc

!ENDIF 

SOURCE=.\stadenmsg.rc
DEP_RSC_STADE=\
	".\icon1.ico"\
	

"$(INTDIR)\stadenmsg.res" : $(SOURCE) $(DEP_RSC_STADE) "$(INTDIR)"
	$(RSC) $(RSC_PROJ) $(SOURCE)



!ENDIF 

