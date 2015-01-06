FISOC
=====

Framework for Ice Sheet Ocean Coupling

The ESMF coupling framework is used.  
Hierarchy is caller -> parent -> everything else.

FISOC_caller.f90       The calling program, registers routines with ESMF
FISOC_parent.f90       ESMF gridded component, parent to ice and ocean components

Note: OM is for Ocean Model, and ISM is for Ice Sheet Model.

FISOC_coupler.f90        ESMF coupler component, handles regridding
FISOC_proc.f90           ESMF gridded component, processes temporal operations on the ice grid
FISOC_ISM_Elmer.f90      ESMF gridded component, ice sheet model
FISOC_OM_ROMS.f90        ESMF gridded component, ocean model
FISOC_ISM_ElmerDummy.f90 ESMF gridded component, for testing in absence of Elmer
FISOC_OM_ROMSdummy.f90   ESMF gridded component, for testing in absence of ROMS


FISOC_config is a Configuration file for FISOC.  Each component may also use its own 
configuration file(s).  ts_ocn_sec is the timestep of the ocean component in seconds.  ts_ratio 
is the ratio of ocean timesteps to ice sheet timesteps.  Example:
&config
     ts_ocn_sec = 2000,   ts_ratio = 4, 
     start_year = 2010,  start_month = 1,
     end_year = 2010,    end_month = 2, 
     tight_coupling = .TRUE./
(here the ice sheet timesteps are 8000 seconds)

Notes on extending the coupling to use different OM and ISM components:
use the dummy fortran modules (or template modules?) as a starting point.  You can add code to 
these without needing to change the soubroutine names.  Follow the ESMF guidance about making 
your component compatible with ESMF.
***more specifics about variables to exchange and about regridding