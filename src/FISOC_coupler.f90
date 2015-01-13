MODULE FISOC_coupler
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_coupler_register
    
  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_register(FISOC_coupler, rc)
    
    TYPE(ESMF_CplComp)  :: FISOC_coupler
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_coupler_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_coupler_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_RUN, &
         userRoutine=FISOC_coupler_run_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_RUN, &
         userRoutine=FISOC_coupler_run_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_coupler_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_coupler_register

  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two phases.  The first phase is to convert 
  ! the ISM state to the OM grid, and the second phase is to convert the OM state 
  ! to the ISM mesh.
  SUBROUTINE FISOC_coupler_init_phase1(FISOC_coupler, ISM_ExpSt, OM_ImpSt, FISOC_clock, rc)

    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: ISM_ExpSt, OM_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_fieldBundle)       :: ISM_ExpFB, OM_ImpFB
    TYPE(ESMF_field)             :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_grid)              :: OM_grid
    TYPE(ESMF_mesh)              :: ISM_mesh
    TYPE(ESMF_config)            :: config
    CHARACTER(len=ESMF_MAXSTR)   :: ISM_name, OM_name
    INTEGER                      :: fieldCount, ii
    TYPE(ESMF_Field),ALLOCATABLE :: fieldList(:)

    rc = ESMF_FAILURE

    ! Establish which ISM and OM components we are using (though we aim to remove dependency 
    ! on this in the coupler if possible)
    CALL ESMF_cplCompGet(FISOC_coupler, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, ISM_name, label='ISM_name:', rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, OM_name, label='OM_name:', rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)


    ! Extract ISM field bundle for regridding...
    CALL ESMF_StateGet(ISM_ExpSt, "ISM export fields", ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=fieldCount, fieldList=fieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    print *,"coupler needs log writes"
    print *,"get OM_grid from OM imp state... via variable?"
    print *,"get ISM_mesh from first ISM field"
    print *,"create a route handle to add to the state objects"


    CALL ESMF_FieldGet(fieldList(1), mesh=ISM_mesh, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    loop_over_fields: DO ii = 1,fieldCount 

       print *,"now regrid the current field"
       
    END DO loop_over_fields

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_coupler_init_phase1


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_init_phase2(FISOC_coupler, OM_ExpSt, ISM_ImpSt, FISOC_clock, rc)
    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: OM_ExpSt, ISM_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_coupler_init_phase2
  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_run_phase1(FISOC_coupler, ISM_ExpSt, OM_ImpSt, FISOC_clock, rc)
    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: ISM_ExpSt, OM_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    print*,"is run phase 1 the same as init phase 1? if so can we call the same routine in both cases?"

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_coupler_run_phase1


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_run_phase2(FISOC_coupler, OM_ExpSt, ISM_ImpSt, FISOC_clock, rc)
    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: OM_ExpSt, ISM_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_coupler_run_phase2


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_finalise(FISOC_coupler, dummy_ExpSt, dummy_ImpSt, FISOC_clock, rc)
    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: dummy_ExpSt, dummy_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_coupler_finalise


  !------------------------------------------------------------------------------
  TYPE(ESMF_grid) FUNCTION dummyCreateGrid(mesh, rc)

    TYPE(ESMF_mesh),INTENT(IN) :: mesh 
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_grid)        :: grid 

    rc = ESMF_FAILURE

    dummyCreateGrid = grid

    rc = ESMF_SUCCESS

  END FUNCTION  dummyCreateGrid

END MODULE FISOC_coupler
