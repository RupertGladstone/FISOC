
! this module is intended to contain any FISOC-specific f90 derived types and global parameters
MODULE FISOC_types_MOD
  
  USE ESMF

  IMPLICIT NONE

!  REAL(ESMF_KIND_R8),PARAMETER :: FISOC_secPerYear = 365.2422 * 24.0 * 60.0 * 60.0
  REAL(ESMF_KIND_R8),PARAMETER :: FISOC_secPerYear = 360.0 * 24.0 * 60.0 * 60.0
  REAL(ESMF_KIND_R8),PARAMETER :: FISOC_missingData = -99999999999.9, FISOC_missing_R8 = -9999.0000000000000000000000000000000000000000000000000
  INTEGER,PARAMETER            :: FISOC_mpic_missing = -99, FISOC_missing = -9999
  INTEGER,PARAMETER            :: OM_outputUnit = 37, ISM_outputUnit = 32
  INTEGER,PARAMETER            :: CLOCKWISE=1, ANTI_CLOCKWISE=2

  ! This time information is set at run time.
  TYPE(ESMF_Time)              :: FISOC_time, FISOC_startTime, FISOC_endTime
!  TYPE(ESMF_TimeInterval)      :: FISOC_OM_TI, FISOC_ISM_TI
  ! TODO: intialise these to zero?

END MODULE FISOC_types_MOD
