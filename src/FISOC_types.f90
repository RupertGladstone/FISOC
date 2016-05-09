
! this module is intended to contain any FISOC-specific derived types and global parameters
MODULE FISOC_types_MOD
  
  USE ESMF

  IMPLICIT NONE

  REAL(ESMF_KIND_R8),PARAMETER :: FISOC_secPerYear = 365.2422 * 24.0 * 60.0 * 60.0
  REAL(ESMF_KIND_R8),PARAMETER :: FISOC_missingData = -99999999999.9
  INTEGER,PARAMETER            :: FISOC_mpic_missing = -99
  INTEGER,PARAMETER            :: OM_outputUnit = 31, ISM_outputUnit = 32

END MODULE FISOC_types_MOD
