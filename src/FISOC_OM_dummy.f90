MODULE FISOC_OM
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_OM_register
    
  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS
  
  SUBROUTINE FISOC_OM_register(FISOC_OM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_OM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_OM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_OM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_OM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_OM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_OM_register


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_init_phase1(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    msg = "OM initialise phase 1 is empty in dummy case"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_init_phase1
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_init_phase2(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_init_phase2
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_run(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_run
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_finalise(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_finalise  

END MODULE FISOC_OM
