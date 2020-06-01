MODULE FISOC_ISM_MOD
  
  USE ESMF
  USE FISOC_ISM_Wrapper
  USE FISOC_utils_MOD
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_ISM_register
  
  CHARACTER(len=ESMF_MAXSTR),SAVE :: ISM_impFBname = "ISM import fields"
  CHARACTER(len=ESMF_MAXSTR),SAVE :: ISM_expFBname = "ISM export fields"

CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_register(FISOC_ISM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_ISM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_ISM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_ISM_register

  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two stages.  The first stage is for independent 
  ! initialisation of the ISM and the second stage occurrs after the OM has been 
  ! initialised, to allow inter-component consistency checks or use of the OM 
  ! state to complete ISM initalisation.
  SUBROUTINE FISOC_ISM_init_phase1(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)        :: FISOC_ISM
    TYPE(ESMF_State)           :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)           :: FISOC_clock
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_config)          :: FISOC_config
    TYPE(ESMF_mesh)            :: ISM_mesh
    TYPE(ESMF_grid)            :: ISM_grid, OM_grid
    TYPE(ESMF_fieldBundle)     :: ISM_ExpFB,OM_ExpFB
    TYPE(ESMF_VM)              :: vm
    TYPE(ESMF_field),ALLOCATABLE :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR) :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: ISM_DerVarList(:)
    INTEGER                    :: fieldCount
    LOGICAL                    :: ISM_UseOMGrid, OM_UseISMGrid

    rc = ESMF_FAILURE

    msg = "ISM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompGet(FISOC_ISM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! if the ISM and OM are on the same grid then the ISM export field bundle is 
    ! used as the OM import field bundle
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_UseOMGrid, 'ISM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_UseISMGrid, 'OM_UseISMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ISM_UseOMGrid.OR.OM_UseISMGrid) THEN
      ISM_impFBname = "OM export fields"
    END IF

    label = 'FISOC_ISM_DerVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, ISM_DerVarList,rc=rc)
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
      ALLOCATE(ISM_DerVarList(0))
    ELSE IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_UseOMGrid, 'ISM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! create empty field bundle
    ISM_ExpFB = ESMF_FieldBundleCreate(name=ISM_expFBname, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! The model-specific initialisation adds the ISM vars to the field bundle.
    ! This is followed by adding non-model-specific derived variables. 
# if defined(FISOC_ISM_GRID)
    IF (ISM_UseOMGrid) THEN
      ! Get grid from field from field list from field bundle from state.
      CALL ESMF_StateGet(ISM_ImpSt, ISM_impFBname, OM_ExpFB, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=fieldCount, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      ALLOCATE(fieldList(fieldCount))
      CALL ESMF_FieldBundleGet(OM_ExpFB, fieldList=fieldList,rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)        
      CALL ESMF_FieldGet(FieldList(1), grid=OM_grid, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      ! We assign the ISM_grid to the OM_grid.  This is not taking a copy, 
      ! but using the same object.
      ISM_grid = OM_grid
      
    END IF
    
    CALL FISOC_ISM_Wrapper_Init_Phase1(FISOC_config,vm,ISM_ExpFB,ISM_grid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_populateFieldBundle(ISM_DerVarList,ISM_ExpFB,ISM_grid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
# elif defined(FISOC_ISM_MESH)
    CALL FISOC_ISM_Wrapper_Init_Phase1(FISOC_config,vm,ISM_ExpFB,ISM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_populateFieldBundle(ISM_DerVarList,ISM_ExpFB,ISM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

# else
    msg = "ERROR: FISOC does not recognise ISM geom type."
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif


    ! Calculate values for derived variables from the model-specific ISM vars.
    CALL FISOC_ISM_calcDerivedFields(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! we only add the fields to the import state as a way of letting the coupler get hold of the 
    ! grid.  There must be a better way to do this. 
    ! [edit: there is now, see email from Gerhard Theurich]
    CALL ESMF_StateAdd(ISM_ImpSt, (/ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(ISM_ExpSt, (/ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "ISM initialise phase 1 complete"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_init_phase1


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_init_phase2(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)

    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_config)      :: FISOC_config
    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_VM)          :: vm

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_ISM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(ISM_ImpSt, ISM_impFBname, ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(ISM_ExpSt, ISM_expFBname, ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_ISM_Wrapper_Init_Phase2(FISOC_config,vm,ISM_ImpFB,ISM_ExpFB,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ISM_maskOMfields(FISOC_config,ISM_ImpFB,ISM_ExpFB,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "ISM initialise phase 2 complete (allows the ISM access to the OM initial state) "
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_init_phase2
  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_run(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)

    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc
    
    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    LOGICAL                :: verbose_coupling
    TYPE(ESMF_config)      :: FISOC_config
!    INTEGER                :: localPet
    TYPE(ESMF_VM)          :: vm

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_ISM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! extract field bundles from import and export states to send to model-specific wrapper

    CALL ESMF_StateGet(ISM_ImpSt, ISM_impFBname, ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(ISM_ExpSt, ISM_expFBname, ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    CALL FISOC_ISM_calcDerivedFields_pre(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_ISM_maskOMfields(FISOC_config,ISM_ImpFB,ISM_ExpFB,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ISM_Wrapper_Run(FISOC_config,vm,ISM_ImpFB=ISM_ImpFB, &
         ISM_ExpFB=ISM_ExpFB,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    CALL FISOC_ISM_calcDerivedFields_post(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    msg = "ISM run complete for current timestep"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_run


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_finalise(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc

    INTEGER                      :: localPet
    TYPE(ESMF_VM)                :: vm
    TYPE(ESMF_config)            :: FISOC_config
    TYPE(ESMF_fieldbundle)       :: ISM_ImpFB, ISM_ExpFB
    INTEGER                      :: FieldCount,ii
    TYPE(ESMF_field),ALLOCATABLE :: ImpFieldList(:)
    TYPE(ESMF_field),ALLOCATABLE :: ExpFieldList(:)
    TYPE(ESMF_mesh)              :: ISM_mesh
    LOGICAL                      :: ISM_UseOMGrid, OM_UseISMGrid

    CHARACTER(len=ESMF_MAXSTR) :: name
    
    rc = ESMF_FAILURE

    msg = "ISM finalise: destroy fields, bundles and mesh"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompGet(FISOC_ISM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_UseOMGrid, 'ISM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_UseISMGrid, 'OM_UseISMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ISM_Wrapper_Finalize(FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    ! destroy ISM export fields
    CALL ESMF_StateGet(ISM_ExpSt, ISM_expFBname, ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       
    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! ... get list of fields from bundle.
    ALLOCATE(ExpFieldList(FieldCount))
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldList=ExpFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,FieldCount
       
       CALL ESMF_FieldGet(field=ExpFieldList(ii),name=name, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldDestroy(ExpFieldList(ii), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO
    
    DEALLOCATE(ExpFieldList)
    
    CALL ESMF_FieldBundleDestroy(ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    
    IF ( (.NOT.ISM_UseOMGrid).AND.(.NOT.OM_UseISMGrid) ) THEN
      ! destroy ISM import fields    
      CALL ESMF_StateGet(ISM_ImpSt, ISM_impFBname, ISM_ImpFB, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
      
      ! ...how many fields?...
      CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldCount=FieldCount, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      ! ... get list of fields from bundle.
      ALLOCATE(ImpFieldList(FieldCount))
      CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldList=ImpFieldList,rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
!    CALL ESMF_FieldGet(ImpFieldList(1), mesh=ISM_mesh, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
      DO ii = 1,FieldCount
        CALL ESMF_FieldDestroy(ImpFieldList(ii), rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      END DO
      
!    CALL ESMF_MeshDestroy(ISM_mesh,rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
      DEALLOCATE(ImpFieldList)
    
      CALL ESMF_FieldBundleDestroy(ISM_ImpFB, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    END IF
      
    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_ISM_finalise
  


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_maskOMfields(FISOC_config,ISM_ImpFB,ISM_ExpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB, ISM_ImpFB
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    LOGICAL                         :: ISM_maskOMvars
    TYPE(ESMF_field)                :: ISM_gmask
    REAL(ESMF_KIND_R8),POINTER      :: mask_ptr(:),ptr(:)
    INTEGER                         :: fieldCount, ii
    TYPE(ESMF_Field),ALLOCATABLE    :: fieldList(:)
!    CHARACTER(len=ESMF_MAXSTR)      :: fieldName

    rc = ESMF_FAILURE

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_maskOMvars, 'ISM_maskOMvars',rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF (.NOT. ISM_maskOMvars) THEN
      rc = ESMF_SUCCESS
      RETURN
    END IF

    ! Get grounded mask from ISM export fields
    CALL ESMF_FieldBundleGet(ISM_ExpFB, "ISM_gmask", field=ISM_gmask, rc=rc)
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
      msg = "ERROR: ISM_maskOMvars set to true but ISM_gmask not available"
    ELSE IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    CALL ESMF_FieldGet(ISM_gmask, farrayPtr=mask_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! Loop through ISM import fields, mask them all
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,fieldCount
      CALL ESMF_FieldGet(fieldList(ii), farrayPtr=ptr, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!       CALL ESMF_FieldGet(fieldList(ii), name=fieldName, rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

      ptr = ptr * mask_ptr
      
      IF (ASSOCIATED(ptr)) THEN
        NULLIFY(ptr)
      END IF
      
    END DO

    IF (ASSOCIATED(mask_ptr)) THEN
      NULLIFY(mask_ptr)
    END IF
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_maskOMfields

  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_calcDerivedFields(ISM_ExpFB,FISOC_config,FISOC_clock,rc)

    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    CALL FISOC_ISM_calcDerivedFields_pre(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    CALL FISOC_ISM_calcDerivedFields_post(ISM_ExpFB,FISOC_config,FISOC_clock,rc)

  END SUBROUTINE FISOC_ISM_calcDerivedFields

  

  !------------------------------------------------------------------------------
  ! Fields that need to be derived before the call to the ISM run method.
  !
  SUBROUTINE FISOC_ISM_calcDerivedFields_pre(ISM_ExpFB,FISOC_config,FISOC_clock,rc)

    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: FISOC_ISM_DerVarList(:),FISOC_ISM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)             :: label
    INTEGER                                :: ii, numDerVars

    rc = ESMF_FAILURE

    ! extract a list of derived ISM variables from the configuration object
    label = 'FISOC_ISM_DerVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_ISM_DerVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    label = 'FISOC_ISM_ReqVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    numDerVars = SIZE(FISOC_ISM_DerVarList)

    DO ii = 1,numDerVars

       SELECT CASE(FISOC_ISM_DerVarList(ii))

       CASE ("ISM_z_l0_previous")
!          CALL ISM_derive_z_l0_previous(ISM_ExpFB,rc)
          CALL ISM_derive_previous("ISM_z_l0", "ISM_z_l0_previous",ISM_ExpFB,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE ("ISM_z_lts_previous")
!          CALL ISM_derive_z_lts_previous(ISM_ExpFB,rc)
          CALL ISM_derive_previous("ISM_z_lts", "ISM_z_lts_previous",ISM_ExpFB,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE ("ISM_dddt","ISM_dsdt","ISM_dTdz_l0","ISM_z_l0_linterp","ISM_z_lts")
          ! Allow vars that will be derived after the ISM run call

       CASE DEFAULT
          msg="ERROR: derived variable name not recognised: "//FISOC_ISM_DerVarList(ii)
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END SELECT

    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_calcDerivedFields_pre



  !------------------------------------------------------------------------------
  ! Fields that need to be derived after the call to the ISM run method.
  !
  SUBROUTINE FISOC_ISM_calcDerivedFields_post(ISM_ExpFB,FISOC_config,FISOC_clock,rc)

    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: FISOC_ISM_DerVarList(:),FISOC_ISM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)             :: label
    INTEGER                                :: ii, numDerVars

    rc = ESMF_FAILURE

    ! extract a list of derived ISM variables from the configuration object
    label = 'FISOC_ISM_DerVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_ISM_DerVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    label = 'FISOC_ISM_ReqVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    numDerVars = SIZE(FISOC_ISM_DerVarList)

    DO ii = 1,numDerVars

       SELECT CASE(FISOC_ISM_DerVarList(ii))

          ! These are handled elsewhere, hence the empty case
       CASE ("ISM_z_l0_previous","ISM_z_lts_previous","ISM_z_l0_linterp")
          ! ignore vars that should have been derived before the ISM run call

       CASE ("ISM_z_lts")
          CALL ISM_derive_z_lts(ISM_ExpFB,FISOC_config,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!       CASE ("ISM_z_lts")
          ! Needed if we have only ice thikness and lower surface height 
          ! from the ice sheet model.
!          CALL ISM_derive_z_lts(ISM_ExpFB,FISOC_config,rc)
!          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!               line=__LINE__, file=__FILE__)) &
!               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE ("ISM_dddt")
!          CALL ISM_derive_rate("ISM_z_l0","ISM_z_l0_previous",  &
!               "ISM_dddt",ISM_ExpFB,FISOC_config,FISOC_clock,rc)
          CALL ISM_derive_dddt(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE ("ISM_dsdt")
          CALL ISM_derive_dsdt(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE ("ISM_dTdz_l0")
          CALL ISM_derive_dTdz_l0(ISM_ExpFB,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE DEFAULT
          msg="ERROR: derived variable name not recognised: "//FISOC_ISM_DerVarList(ii)
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END SELECT

    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_calcDerivedFields_post



  !------------------------------------------------------------------------------
  ! If ice thickness and lower surface elevation are available we can calculate 
  ! the upper surface elevation.
  SUBROUTINE ISM_derive_z_lts(ISM_ExpFB,FISOC_config,rc)
    
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER, INTENT(OUT),OPTIONAL         :: rc
    
    TYPE(ESMF_field)           :: ISM_z_lts, ISM_thick, ISM_z_l0
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_lts_ptr(:), ISM_thick_ptr(:), &
         ISM_z_l0_ptr(:)


    RC = ESMF_FAILURE

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_lts", field=ISM_z_lts, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_lts, localDe=0, farrayPtr=ISM_z_lts_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_thick", field=ISM_thick, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_thick, localDe=0, farrayPtr=ISM_thick_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_z_lts_ptr = ISM_z_l0_ptr + ISM_thick_ptr

    NULLIFY(ISM_z_lts_ptr)
    NULLIFY(ISM_z_l0_ptr)
    NULLIFY(ISM_thick_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_z_lts


  !------------------------------------------------------------------------------
  ! Copy current data.  At next timestep this copy will contain "previous" data.
  !
  SUBROUTINE ISM_derive_previous(varName,varName_previous,ISM_ExpFB,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL         :: rc
    CHARACTER(len=*),INTENT(IN)           :: varName,varName_previous

    TYPE(ESMF_field)           :: var, var_previous

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName=TRIM(varName), field=var, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName=TRIM(varName_previous), field=var_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldCopy(var_previous, var, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_previous


  !------------------------------------------------------------------------------
  ! Copy current data.  At next timestep this copy will contain "previous" data.
  !
  ! TODO: remove this once the generic routine derive_previous is tested
  SUBROUTINE ISM_derive_z_l0_previous(ISM_ExpFB,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)            :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL                   :: rc

    TYPE(ESMF_field)           :: ISM_z_l0
    TYPE(ESMF_field)           :: ISM_z_l0_previous

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0_previous", field=ISM_z_l0_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldCopy(ISM_z_l0_previous, ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_z_l0_previous


  !------------------------------------------------------------------------------
  ! Copy current data.  At next timestep this copy will contain "previous" data.
  !
  SUBROUTINE ISM_derive_z_lts_previous(ISM_ExpFB,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)            :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL                   :: rc

    TYPE(ESMF_field)           :: ISM_z_lts
    TYPE(ESMF_field)           :: ISM_z_lts_previous
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_lts_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_lts_previous_ptr(:) 

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_lts", field=ISM_z_lts, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_lts_previous", field=ISM_z_lts_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldCopy(ISM_z_lts_previous, ISM_z_lts, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_z_lts_previous


  !------------------------------------------------------------------------------
  ! TODO: a lot of code duplication in the dddt and dsdt routines could do with 
  ! rationalising.

  !------------------------------------------------------------------------------
  ! dsdt, short for ds/dt, is the rate of change of upper ice surface with time, 
  ! for use by ocean model.
  !
  SUBROUTINE ISM_derive_dsdt(ISM_ExpFB,FISOC_config,FISOC_clock,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc

    INTEGER                :: rank
    
    rc = ESMF_FAILURE
    rank = 0

    CALL FISOC_getFirstFieldRank(ISM_ExpFB,rank,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    SELECT CASE (rank)

    CASE(1)
       CALL  ISM_derive_dsdt_1D(ISM_ExpFB,FISOC_config,FISOC_clock,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CASE(2)
       CALL  ISM_derive_dsdt_2D(ISM_ExpFB,FISOC_config,FISOC_clock,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CASE DEFAULT
       msg = "ERROR: dimension of array ptr NYI"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_dsdt

  
  !------------------------------------------------------------------------------
  SUBROUTINE ISM_derive_dsdt_1D(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    TYPE(ESMF_field)           :: ISM_z_lts
    TYPE(ESMF_field)           :: ISM_z_lts_previous
    TYPE(ESMF_field)           :: ISM_dsdt
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_lts_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_lts_previous_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dsdt_ptr(:) 
    REAL(ESMF_KIND_R8)         :: ISM_dt
    INTEGER                    :: ISM_dt_int
    INTEGER(ESMF_KIND_I8)      :: advancecount

    RC = ESMF_FAILURE

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_lts", field=ISM_z_lts, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_lts, localDe=0, farrayPtr=ISM_z_lts_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_lts_previous", field=ISM_z_lts_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_lts_previous, localDe=0, farrayPtr=ISM_z_lts_previous_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dsdt", field=ISM_dsdt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dsdt, localDe=0, farrayPtr=ISM_dsdt_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_int, 'ISM_dt_sec',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_dt = REAL(ISM_dt_int,ESMF_KIND_R8)
    ISM_dsdt_ptr = (ISM_z_lts_ptr - ISM_z_lts_previous_ptr) / ISM_dt

    ! first time we set it to zero
    CALL ESMF_ClockGet(FISOC_clock, advanceCount=advanceCount, rc=rc)
    IF (advanceCount.LE.1) THEN
       ISM_dsdt_ptr = 0.0
    END IF
!    WHERE (ISM_z_lts_previous_ptr .EQ. 0.0) ISM_dsdt_ptr = 0.0

    NULLIFY(ISM_z_lts_previous_ptr)
    NULLIFY(ISM_z_lts_ptr)
    NULLIFY(ISM_dsdt_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_dsdt_1D


  !------------------------------------------------------------------------------
  SUBROUTINE ISM_derive_dsdt_2D(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    TYPE(ESMF_field)           :: ISM_z_lts
    TYPE(ESMF_field)           :: ISM_z_lts_previous
    TYPE(ESMF_field)           :: ISM_dsdt
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_lts_ptr(:,:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_lts_previous_ptr(:,:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dsdt_ptr(:,:) 
    REAL(ESMF_KIND_R8)         :: ISM_dt
    INTEGER                    :: ISM_dt_int
    INTEGER(ESMF_KIND_I8)      :: advancecount

    RC = ESMF_FAILURE

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_lts", field=ISM_z_lts, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_lts, localDe=0, farrayPtr=ISM_z_lts_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_lts_previous", field=ISM_z_lts_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_lts_previous, localDe=0, farrayPtr=ISM_z_lts_previous_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dsdt", field=ISM_dsdt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dsdt, localDe=0, farrayPtr=ISM_dsdt_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_int, 'ISM_dt_sec',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_dt = REAL(ISM_dt_int,ESMF_KIND_R8)
    ISM_dsdt_ptr = (ISM_z_lts_ptr - ISM_z_lts_previous_ptr) / ISM_dt

    ! first time we set it to zero
    CALL ESMF_ClockGet(FISOC_clock, advanceCount=advanceCount, rc=rc)
    IF (advanceCount.LE.1) THEN
       ISM_dsdt_ptr = 0.0
    END IF

    NULLIFY(ISM_z_lts_previous_ptr)
    NULLIFY(ISM_z_lts_ptr)
    NULLIFY(ISM_dsdt_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_dsdt_2D


  !------------------------------------------------------------------------------
  !
  SUBROUTINE ISM_derive_rate(varName,varName_previous,varName_rate, &
       ISM_ExpFB,FISOC_config,FISOC_clock,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    CHARACTER(len=*),INTENT(IN)           :: varName,varName_previous,varName_rate
    TYPE(ESMF_clock),INTENT(IN)           :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL         :: rc

    INTEGER                :: rank
    
    rc = ESMF_FAILURE
    rank = 0

    CALL FISOC_getFirstFieldRank(ISM_ExpFB,rank,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    SELECT CASE (rank)

    CASE(1)
       CALL  ISM_derive_rate_1D(varName,varName_previous,varName_rate, &
            ISM_ExpFB,FISOC_config,FISOC_clock,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CASE(2)
       CALL  ISM_derive_rate_2D(varName,varName_previous,varName_rate, &
            ISM_ExpFB,FISOC_config,FISOC_clock,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CASE DEFAULT
       msg = "ERROR: dimension of array ptr NYI"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_rate


  !------------------------------------------------------------------------------
  ! dddt, short for dd/dt, is the rate of change of depth (or draft) with time, 
  ! for use by ocean model.
  !
  SUBROUTINE ISM_derive_dddt(ISM_ExpFB,FISOC_config,FISOC_clock,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc

    INTEGER                :: rank
    
    rc = ESMF_FAILURE
    rank = 0

    CALL FISOC_getFirstFieldRank(ISM_ExpFB,rank,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    SELECT CASE (rank)

    CASE(1)
       CALL  ISM_derive_dddt_1D(ISM_ExpFB,FISOC_config,FISOC_clock,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CASE(2)
       CALL  ISM_derive_dddt_2D(ISM_ExpFB,FISOC_config,FISOC_clock,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CASE DEFAULT
       msg = "ERROR: dimension of array ptr NYI"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_dddt

  
  !------------------------------------------------------------------------------
  SUBROUTINE ISM_derive_rate_1D(varName,varName_previous,varName_rate, &
       ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    

    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    CHARACTER(len=*),INTENT(IN)            :: varName,varName_previous
    CHARACTER(len=*),INTENT(IN)            :: varName_rate
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    TYPE(ESMF_field)           :: field
    TYPE(ESMF_field)           :: field_previous
    TYPE(ESMF_field)           :: rate
    REAL(ESMF_KIND_R8),POINTER :: field_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: field_previous_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: rate_ptr(:) 
    REAL(ESMF_KIND_R8)         :: ISM_dt
    INTEGER                    :: ISM_dt_int
    INTEGER(ESMF_KIND_I8)      :: advancecount

    RC = ESMF_FAILURE

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName=varName, field=field, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=field, localDe=0, farrayPtr=field_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName=varName_previous, field=field_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=field_previous, localDe=0, farrayPtr=field_previous_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName=varName_rate, field=rate, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=rate, localDe=0, farrayPtr=rate_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_int, 'ISM_dt_sec',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_dt = REAL(ISM_dt_int,ESMF_KIND_R8)
    rate_ptr = (field_ptr - field_previous_ptr) / ISM_dt

    ! first time we set it to zero
    CALL ESMF_ClockGet(FISOC_clock, advanceCount=advanceCount, rc=rc)
    IF (advanceCount.LE.1) THEN
       rate_ptr = 0.0
    END IF

    NULLIFY(field_previous_ptr)
    NULLIFY(field_ptr)
    NULLIFY(rate_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_rate_1D


  !------------------------------------------------------------------------------
  SUBROUTINE ISM_derive_rate_2D(varName,varName_previous,varName_rate, &
       ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    

    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    CHARACTER(len=*),INTENT(IN)            :: varName,varName_previous
    CHARACTER(len=*),INTENT(IN)            :: varName_rate
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    TYPE(ESMF_field)           :: field
    TYPE(ESMF_field)           :: field_previous
    TYPE(ESMF_field)           :: rate
    REAL(ESMF_KIND_R8),POINTER :: field_ptr(:,:) 
    REAL(ESMF_KIND_R8),POINTER :: field_previous_ptr(:,:) 
    REAL(ESMF_KIND_R8),POINTER :: rate_ptr(:,:) 
    REAL(ESMF_KIND_R8)         :: ISM_dt
    INTEGER                    :: ISM_dt_int
    INTEGER(ESMF_KIND_I8)      :: advancecount

    RC = ESMF_FAILURE

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName=varName, field=field, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=field, localDe=0, farrayPtr=field_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName=varName_previous, field=field_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=field_previous, localDe=0, farrayPtr=field_previous_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName=varName_rate, field=rate, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=rate, localDe=0, farrayPtr=rate_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_int, 'ISM_dt_sec',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_dt = REAL(ISM_dt_int,ESMF_KIND_R8)
    rate_ptr = (field_ptr - field_previous_ptr) / ISM_dt

    ! first time we set it to zero
    CALL ESMF_ClockGet(FISOC_clock, advanceCount=advanceCount, rc=rc)
    IF (advanceCount.LE.1) THEN
       rate_ptr = 0.0
    END IF

    NULLIFY(field_previous_ptr)
    NULLIFY(field_ptr)
    NULLIFY(rate_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_rate_2D


  !------------------------------------------------------------------------------
  SUBROUTINE ISM_derive_dddt_1D(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    TYPE(ESMF_field)           :: ISM_z_l0
    TYPE(ESMF_field)           :: ISM_z_l0_previous
    TYPE(ESMF_field)           :: ISM_dddt
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_previous_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dddt_ptr(:) 
    REAL(ESMF_KIND_R8)         :: ISM_dt
    INTEGER                    :: ISM_dt_int
    INTEGER(ESMF_KIND_I8)      :: advancecount

    RC = ESMF_FAILURE

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0_previous", field=ISM_z_l0_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0_previous, localDe=0, farrayPtr=ISM_z_l0_previous_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dddt", field=ISM_dddt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dddt, localDe=0, farrayPtr=ISM_dddt_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_int, 'ISM_dt_sec',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_dt = REAL(ISM_dt_int,ESMF_KIND_R8)
    ISM_dddt_ptr = (ISM_z_l0_ptr - ISM_z_l0_previous_ptr) / ISM_dt

    ! first time we set it to zero
    CALL ESMF_ClockGet(FISOC_clock, advanceCount=advanceCount, rc=rc)
    IF (advanceCount.LE.1) THEN
       ISM_dddt_ptr = 0.0
    END IF

    NULLIFY(ISM_z_l0_previous_ptr)
    NULLIFY(ISM_z_l0_ptr)
    NULLIFY(ISM_dddt_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_dddt_1D


  !------------------------------------------------------------------------------
  SUBROUTINE ISM_derive_dddt_2D(ISM_ExpFB,FISOC_config,FISOC_clock,rc)
    
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_clock),INTENT(IN)            :: FISOC_clock
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    TYPE(ESMF_field)           :: ISM_z_l0
    TYPE(ESMF_field)           :: ISM_z_l0_previous
    TYPE(ESMF_field)           :: ISM_dddt
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:,:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_previous_ptr(:,:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dddt_ptr(:,:) 
    REAL(ESMF_KIND_R8)         :: ISM_dt
    INTEGER                    :: ISM_dt_int
    INTEGER(ESMF_KIND_I8)      :: advancecount

    RC = ESMF_FAILURE

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0_previous", field=ISM_z_l0_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0_previous, localDe=0, farrayPtr=ISM_z_l0_previous_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dddt", field=ISM_dddt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dddt, localDe=0, farrayPtr=ISM_dddt_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_int, 'ISM_dt_sec',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_dt = REAL(ISM_dt_int,ESMF_KIND_R8)
    ISM_dddt_ptr = (ISM_z_l0_ptr - ISM_z_l0_previous_ptr) / ISM_dt

    ! first time we set it to zero
    CALL ESMF_ClockGet(FISOC_clock, advanceCount=advanceCount, rc=rc)
    IF (advanceCount.LE.1) THEN
       ISM_dddt_ptr = 0.0
    END IF

    NULLIFY(ISM_z_l0_previous_ptr)
    NULLIFY(ISM_z_l0_ptr)
    NULLIFY(ISM_dddt_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_dddt_2D



  !------------------------------------------------------------------------------
  SUBROUTINE ISM_derive_dTdz_l0(ISM_ExpFB,rc)

!    CHARACTER(len=ESMF_MAXSTR),INTENT(IN)   :: FISOC_ISM_ReqVarList(:)
    TYPE(ESMF_fieldBundle),INTENT(INOUT)    :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL           :: rc

    TYPE(ESMF_field)           :: ISM_z_l0,ISM_z_l1
    TYPE(ESMF_field)           :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_field)           :: ISM_dTdz_l0
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:), ISM_z_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_temperature_l0_ptr(:),ISM_temperature_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:)

    INTEGER :: localDeCount

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_temperature_l0", field=ISM_temperature_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l0, localDe=0, farrayPtr=ISM_temperature_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_temperature_l1", field=ISM_temperature_l1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l1, localDe=0, farrayPtr=ISM_temperature_l1_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldGet(field=ISM_z_l0, localDeCount=localDeCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l1", field=ISM_z_l1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l1, localDe=0, farrayPtr=ISM_z_l1_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dTdz_l0", field=ISM_dTdz_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dTdz_l0, localDe=0, farrayPtr=ISM_dTdz_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! the simplest approximation for the vertical temperature gradient at the ice base is to 
    ! use the lowest two levels and assume the gradient doesn't change between these two 
    ! levels.
    ISM_dTdz_l0_ptr(:) = (ISM_temperature_l1_ptr-ISM_temperature_l0_ptr) / (ISM_z_l1_ptr-ISM_z_l0_ptr)

    NULLIFY(ISM_dTdz_l0_ptr)
    NULLIFY(ISM_temperature_l1_ptr)
    NULLIFY(ISM_temperature_l0_ptr)
    NULLIFY(ISM_z_l1_ptr)
    NULLIFY(ISM_z_l0_ptr)

    rc = ESMF_SUCCESS
    
  END SUBROUTINE ISM_derive_dTdz_l0
  
END MODULE FISOC_ISM_MOD
