# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

!IF "$(CFG)" == ""
CFG=gpoly - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to gpoly - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "gpoly - Win32 Release" && "$(CFG)" != "gpoly - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "gpoly.mak" CFG="gpoly - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "gpoly - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "gpoly - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "gpoly - Win32 Debug"
F90=fl32.exe

!IF  "$(CFG)" == "gpoly - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "."
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\.
INTDIR=.\Release

ALL : "$(OUTDIR)\rfunc_7_9.lib"

CLEAN : 
	-@erase "..\rfunc_7_9.lib"
	-@erase ".\Release\fgpoly.obj"
	-@erase ".\Release\gpoly_gradients.obj"
	-@erase ".\Release\gpoly_surfaces.obj"
	-@erase ".\Release\rfunc_7_9.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /4R8 /I "Release/" /c /nologo
F90_PROJ=/Ox /4R8 /I "Release/" /c /nologo /Fo"Release/" 
F90_OBJS=.\Release/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/gpoly.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\rfunc_7_9.lib"
LIB32_FLAGS=/nologo /out:"..\rfunc_7_9.lib" 
LIB32_OBJS= \
	"$(INTDIR)/fgpoly.obj" \
	"$(INTDIR)/gpoly_gradients.obj" \
	"$(INTDIR)/gpoly_surfaces.obj" \
	"$(INTDIR)/rfunc_7_9.obj" \
	"..\..\..\..\..\software\msdev\LIB\MATHD.LIB"

"$(OUTDIR)\rfunc_7_9.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "gpoly - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "."
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\.
INTDIR=.\Debug

ALL : "$(OUTDIR)\rfunc_7_9.lib"

CLEAN : 
	-@erase "..\rfunc_7_9.lib"
	-@erase ".\Debug\gpoly_surfaces.obj"
	-@erase ".\Debug\fgpoly.obj"
	-@erase ".\Debug\gpoly_gradients.obj"
	-@erase ".\Debug\rfunc_7_9.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "Debug/" /c /nologo
# ADD F90 /4R8 /Z7 /I "Debug/" /c /nologo
F90_PROJ=/4R8 /Z7 /I "Debug/" /c /nologo /Fo"Debug/" 
F90_OBJS=.\Debug/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/gpoly.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\rfunc_7_9.lib"
LIB32_FLAGS=/nologo /out:"..\rfunc_7_9.lib" 
LIB32_OBJS= \
	"$(INTDIR)/gpoly_surfaces.obj" \
	"$(INTDIR)/fgpoly.obj" \
	"$(INTDIR)/gpoly_gradients.obj" \
	"$(INTDIR)/rfunc_7_9.obj" \
	"..\..\..\..\..\software\msdev\LIB\MATHD.LIB"

"$(OUTDIR)\rfunc_7_9.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "gpoly - Win32 Release"
# Name "gpoly - Win32 Debug"

!IF  "$(CFG)" == "gpoly - Win32 Release"

!ELSEIF  "$(CFG)" == "gpoly - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\gpoly_surfaces.for

"$(INTDIR)\gpoly_surfaces.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\fgpoly.for

"$(INTDIR)\fgpoly.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\gpoly_gradients.f90

"$(INTDIR)\gpoly_gradients.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=\software\msdev\LIB\MATHD.LIB

!IF  "$(CFG)" == "gpoly - Win32 Release"

!ELSEIF  "$(CFG)" == "gpoly - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\rfunc_7_9.f90

"$(INTDIR)\rfunc_7_9.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
