MODULE FISOC_ISM
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_ISM_register
  
CONTAINS
  
  SUBROUTINE FISOC_ISM_register(FISOC_ISM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_ISM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_ISM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_ISM_register

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_init(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)     :: FISOC_clock
    TYPE(ESMF_config)    :: config
    TYPE(ESMF_mesh)      :: ISM_mesh
    INTEGER              :: localrc
    INTEGER, INTENT(OUT) :: rc
    CHARACTER(len=ESMF_MAXSTR) :: ISM_meshFile 
   
    rc = ESMF_SUCCESS

    ! Get mesh file name from the config file
    CALL ESMF_GridCompGet(FISOC_ISM, config=config, rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    CALL ESMF_ConfigGetAttribute(config, ISM_meshFile, label='ISM_meshFile:', rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! create ISM mesh and use it to create zeroed fields for the ISM import and export states
    ISM_mesh = ESMF_MeshCreate(filename=ISM_meshFile, &
            filetypeflag=ESMF_FILEFORMAT_ESMFMESH, &
            rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

print *,'!create mesh here and use it to create fields#'
print *,'!put the fields in the state objects, set to zero for now probably'
print *,'need some log file writes here and in parent'

  END SUBROUTINE FISOC_ISM_init
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_run(FISOC_ISM, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount, localrc
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_SUCCESS

    print *," FISOC_ISM_run needs writing"
    
  END SUBROUTINE FISOC_ISM_run
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_finalise(FISOC_ISM, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount, localrc
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_SUCCESS

    print *," FISOC_ISM_finalise needs writing"
    
  END SUBROUTINE FISOC_ISM_finalise
  
END MODULE FISOC_ISM
