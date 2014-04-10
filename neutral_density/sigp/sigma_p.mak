# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

!IF "$(CFG)" == ""
CFG=sigma_p - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to sigma_p - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "sigma_p - Win32 Release" && "$(CFG)" !=\
 "sigma_p - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "sigma_p.mak" CFG="sigma_p - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "sigma_p - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "sigma_p - Win32 Debug" (based on "Win32 (x86) Static Library")
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
# PROP Target_Last_Scanned "sigma_p - Win32 Debug"
F90=fl32.exe

!IF  "$(CFG)" == "sigma_p - Win32 Release"

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

ALL : "$(OUTDIR)\sigma_p.lib"

CLEAN : 
	-@erase ".\sigma_p.lib"
	-@erase ".\Release\sigp_gradients.obj"
	-@erase ".\Release\sigp_surface.obj"
	-@erase ".\Release\sigp_surfaces.obj"
	-@erase ".\Release\fsigp.obj"
	-@erase ".\Release\sctp_interp.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /I "Release/" /c /nologo
F90_PROJ=/Ox /I "Release/" /c /nologo /Fo"Release/" 
F90_OBJS=.\Release/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/sigma_p.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/sigma_p.lib" 
LIB32_OBJS= \
	"$(INTDIR)/sigp_gradients.obj" \
	"$(INTDIR)/sigp_surface.obj" \
	"$(INTDIR)/sigp_surfaces.obj" \
	"$(INTDIR)/fsigp.obj" \
	"$(INTDIR)/sctp_interp.obj" \
	"..\..\..\..\software\msdev\LIB\MATHD.LIB"

"$(OUTDIR)\sigma_p.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "sigma_p - Win32 Debug"

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

ALL : "$(OUTDIR)\sigma_p.lib"

CLEAN : 
	-@erase ".\sigma_p.lib"
	-@erase ".\Debug\sigp_gradients.obj"
	-@erase ".\Debug\sigp_surface.obj"
	-@erase ".\Debug\sigp_surfaces.obj"
	-@erase ".\Debug\fsigp.obj"
	-@erase ".\Debug\sctp_interp.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "Debug/" /c /nologo
# ADD F90 /Z7 /I "Debug/" /c /nologo
F90_PROJ=/Z7 /I "Debug/" /c /nologo /Fo"Debug/" 
F90_OBJS=.\Debug/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/sigma_p.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/sigma_p.lib" 
LIB32_OBJS= \
	"$(INTDIR)/sigp_gradients.obj" \
	"$(INTDIR)/sigp_surface.obj" \
	"$(INTDIR)/sigp_surfaces.obj" \
	"$(INTDIR)/fsigp.obj" \
	"$(INTDIR)/sctp_interp.obj" \
	"..\..\..\..\software\msdev\LIB\MATHD.LIB"

"$(OUTDIR)\sigma_p.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
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

# Name "sigma_p - Win32 Release"
# Name "sigma_p - Win32 Debug"

!IF  "$(CFG)" == "sigma_p - Win32 Release"

!ELSEIF  "$(CFG)" == "sigma_p - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\sigp_surfaces.f90

"$(INTDIR)\sigp_surfaces.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\sigp_surface.f90

"$(INTDIR)\sigp_surface.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\sigp_gradients.f90

"$(INTDIR)\sigp_gradients.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\fsigp.f90

"$(INTDIR)\fsigp.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=\software\msdev\LIB\MATHD.LIB

!IF  "$(CFG)" == "sigma_p - Win32 Release"

!ELSEIF  "$(CFG)" == "sigma_p - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\sctp_interp.f90

"$(INTDIR)\sctp_interp.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
