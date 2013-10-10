FISOC
=====

Framework for ice sheet ocean coupling

The ESMF coupling framework is used.  
Hierarchy is caller -> parent -> everything else.
FISOC_caller.f90       The calling program, registers routines with ESMF
FISOC_parent.f90       ESMF gridded component, parent to ice and ocean components
FISOC_coupler.f90      ESMF coupler component, handles regridding
FISOC_proc.f90         ESMF gridded component, processes temporal operations on the ice grid
FISOC_Elmer.f90        ESMF gridded component, ice sheet model
FISOC_ROMS.f90         ESMF gridded component, ocean model
FISOC_Elmer_dummy.f90  ESMF gridded component, for testing in absence of Elmer
FISOC_ROMS_dummy.f90   ESMF gridded component, for testing in absence of ROMS
