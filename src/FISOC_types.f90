
! this module is intended to contain any FISOC-specific derived types and global parameters
MODULE FISOC_types_MOD
  
  USE ESMF

  IMPLICIT NONE

  REAL(ESMF_KIND_R8),PARAMETER :: FISOC_missingData = -99999999999.9
  INTEGER,PARAMETER            :: FISOC_mpic_missing = -99

END MODULE FISOC_types_MOD
