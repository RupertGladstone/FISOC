MODULE FISOC_ISM
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_ISM_register
  
!  TYPE(ESMF_GridComp), SAVE :: FISOC_ice, FISOC_ocean, FISOC_proc
!  TYPE(ESMF_CplComp),  SAVE :: FISOC_coupler
!  TYPE(ESMF_State),    SAVE :: FISOC_ice_import, FISOC_ice_export
  
CONTAINS
  
  SUBROUTINE FISOC_ISM_register(FISOC_ISM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
!    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_RUN, &
!         userRoutine=FISOC_ISM_run, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) RETURN
    
!    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_FINALIZE, &
!         userRoutine=FISOC_ISM_finalize, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) RETURN

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_ISM_register

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_init(FISOC_ISM, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount, localrc
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_SUCCESS

print *,"************************************"

  END SUBROUTINE FISOC_ISM_init

END MODULE FISOC_ISM
