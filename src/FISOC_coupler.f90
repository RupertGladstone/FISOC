MODULE FISOC_coupler
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_coupler_register
    
  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS
  
  SUBROUTINE FISOC_coupler_register(FISOC_coupler, rc)
    
    TYPE(ESMF_CplComp)  :: FISOC_coupler
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_coupler_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_coupler_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_RUN, &
         userRoutine=FISOC_coupler_run_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_RUN, &
         userRoutine=FISOC_coupler_run_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_coupler_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

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

    TYPE(ESMF_fieldBundle) :: ISM_ExpFB, OM_ImpFB
    TYPE(ESMF_field)       :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_grid)        :: OM_grid
    TYPE(ESMF_config)      :: config
    CHARACTER(len=ESMF_MAXSTR) :: ISM_name, OM_name


     integer                    :: fieldCount
     type(ESMF_Field),ALLOCATABLE           :: fieldList(:)
     character(len=ESMF_MAXSTR),ALLOCATABLE :: fieldNameList(:)
!     character(len=*),          :: name
    rc = ESMF_FAILURE



!     subroutine ESMF_FieldBundleGetListAll(fieldbundle, &
!       itemorderflag, geomtype, grid, locstream, mesh, xgrid, &
!       fieldCount, fieldList, fieldNameList, name, rc)

! make array of fields so we dont have to know all their names here in the coupler? 
! instead we just regrid them all?

    ! Establish which ISM and OM components we are using (though we aim to remove dependency 
    ! on this in the coupler if possible)
    CALL ESMF_cplCompGet(FISOC_coupler, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    CALL ESMF_ConfigGetAttribute(config, ISM_name, label='ISM_name:', rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, OM_name, label='OM_name:', rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)


    ! Extract ISM field bundle for regridding...
    CALL ESMF_StateGet(ISM_ExpSt, "ISM export fields", ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=fieldCount, fieldList=fieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return


    ! the current dummy coupling creates an OM grid on the fly from the ISM mesh
    IF ((ISM_name.EQ."dummy").AND.(OM_name.EQ."dummy")) THEN
print*,"get mesh from field"
print*,"use mesh to make grid"
!fieldList(1)
!       OM_grid = function***(ISM_mesh)
    ELSE
       print*,"get OM_grid from OM state"
    END IF

!***coupler needs both grids
!***therefore both OM and ISM must have initialisation phase 1 and 2 (OM1 > ISM1 > cpl1> OM2 > cpl2 > ISM2
!print*,fieldNameList(:)
!print*,size(fieldNameList)

!    CALL ESMF_StateGet(state, "ISM export fields", ISM_temperature_l0, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__, rcToReturn=rc)) return




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
